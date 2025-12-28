// gpu_field.rs - GPU-accelerated vortex field simulation
//
// The O(n²) pairwise interactions are perfect for GPU parallelization.
// Each GPU thread computes forces on one vortex from all others.
//
// For 50k vortices: 2.5 billion interactions per tick
// GPU can do this in milliseconds instead of minutes.

use std::collections::HashMap;
use bytemuck::{Pod, Zeroable};
use wgpu::util::DeviceExt;

use crate::vortex::{Vortex, OrbitalPair, Interaction, determine_interaction, merge_vortices};

/// GPU-compatible vortex data (must be Pod for GPU transfer)
#[repr(C)]
#[derive(Copy, Clone, Debug, Pod, Zeroable)]
pub struct GpuVortex {
    pub x: f32,
    pub y: f32,
    pub z: f32,
    pub vx: f32,
    pub vy: f32,
    pub vz: f32,
    pub energy: f32,
    pub frequency: f32,
    pub phase: f32,
    pub spin: f32,  // +1.0 or -1.0
    pub _pad1: f32,
    pub _pad2: f32,
}

/// GPU-computed acceleration output (3D)
#[repr(C)]
#[derive(Copy, Clone, Debug, Pod, Zeroable)]
pub struct GpuAccel {
    pub ax: f32,
    pub ay: f32,
    pub az: f32,
    pub _pad: f32,
}

/// Simulation parameters for GPU
#[repr(C)]
#[derive(Copy, Clone, Debug, Pod, Zeroable)]
pub struct SimParams {
    pub dt: f32,
    pub num_vortices: u32,
    pub field_size: f32,
    pub _pad: f32,
}

/// GPU-accelerated field
pub struct GpuField {
    // GPU resources
    device: wgpu::Device,
    queue: wgpu::Queue,
    compute_pipeline: wgpu::ComputePipeline,
    bind_group_layout: wgpu::BindGroupLayout,
    
    // Simulation state (CPU side, synced with GPU)
    pub vortices: HashMap<u64, Vortex>,
    pub orbital_pairs: Vec<OrbitalPair>,
    pub next_id: u64,
    
    // Void positions (absolute nothing - blocks bonds)
    pub voids: Vec<(f64, f64, f64, f64)>,  // (x, y, z, radius)
    
    // Field properties
    pub size: f64,
    pub orbit_threshold: f64,
    
    // Statistics
    pub tick_count: u64,
    pub total_merges: u64,
    pub total_orbits_formed: u64,
    pub bonds_blocked_by_void: u64,  // Track LOS failures
}

/// Maximum bonds per vortex - the kissing number for 3D spheres
/// Prevents singularity-like over-connection
pub const MAX_BONDS_PER_VORTEX: usize = 13;

impl GpuField {
    /// Create a new GPU-accelerated field
    pub async fn new_async(size: f64) -> Self {
        // Initialize GPU
        let instance = wgpu::Instance::new(wgpu::InstanceDescriptor {
            backends: wgpu::Backends::all(),
            ..Default::default()
        });
        
        let adapter = instance
            .request_adapter(&wgpu::RequestAdapterOptions {
                power_preference: wgpu::PowerPreference::HighPerformance,
                compatible_surface: None,
                force_fallback_adapter: false,
            })
            .await
            .expect("Failed to find GPU adapter");
        
        println!("🎮 GPU: {}", adapter.get_info().name);
        println!("   Backend: {:?}", adapter.get_info().backend);
        
        let (device, queue) = adapter
            .request_device(
                &wgpu::DeviceDescriptor {
                    label: Some("Vortex GPU"),
                    required_features: wgpu::Features::empty(),
                    required_limits: wgpu::Limits::default(),
                    memory_hints: wgpu::MemoryHints::Performance,
                },
                None,
            )
            .await
            .expect("Failed to create GPU device");
        
        // Create compute shader
        let shader = device.create_shader_module(wgpu::ShaderModuleDescriptor {
            label: Some("Vortex Compute Shader"),
            source: wgpu::ShaderSource::Wgsl(COMPUTE_SHADER.into()),
        });
        
        // Bind group layout
        let bind_group_layout = device.create_bind_group_layout(&wgpu::BindGroupLayoutDescriptor {
            label: Some("Vortex Bind Group Layout"),
            entries: &[
                // Params
                wgpu::BindGroupLayoutEntry {
                    binding: 0,
                    visibility: wgpu::ShaderStages::COMPUTE,
                    ty: wgpu::BindingType::Buffer {
                        ty: wgpu::BufferBindingType::Uniform,
                        has_dynamic_offset: false,
                        min_binding_size: None,
                    },
                    count: None,
                },
                // Vortices (input)
                wgpu::BindGroupLayoutEntry {
                    binding: 1,
                    visibility: wgpu::ShaderStages::COMPUTE,
                    ty: wgpu::BindingType::Buffer {
                        ty: wgpu::BufferBindingType::Storage { read_only: true },
                        has_dynamic_offset: false,
                        min_binding_size: None,
                    },
                    count: None,
                },
                // Accelerations (output)
                wgpu::BindGroupLayoutEntry {
                    binding: 2,
                    visibility: wgpu::ShaderStages::COMPUTE,
                    ty: wgpu::BindingType::Buffer {
                        ty: wgpu::BufferBindingType::Storage { read_only: false },
                        has_dynamic_offset: false,
                        min_binding_size: None,
                    },
                    count: None,
                },
            ],
        });
        
        let pipeline_layout = device.create_pipeline_layout(&wgpu::PipelineLayoutDescriptor {
            label: Some("Vortex Pipeline Layout"),
            bind_group_layouts: &[&bind_group_layout],
            push_constant_ranges: &[],
        });
        
        let compute_pipeline = device.create_compute_pipeline(&wgpu::ComputePipelineDescriptor {
            label: Some("Vortex Compute Pipeline"),
            layout: Some(&pipeline_layout),
            module: &shader,
            entry_point: "main",
            compilation_options: Default::default(),
            cache: None,
        });
        
        Self {
            device,
            queue,
            compute_pipeline,
            bind_group_layout,
            vortices: HashMap::new(),
            orbital_pairs: Vec::new(),
            next_id: 1,
            voids: Vec::new(),
            size,
            orbit_threshold: size / 20.0,
            tick_count: 0,
            total_merges: 0,
            total_orbits_formed: 0,
            bonds_blocked_by_void: 0,
        }
    }
    
    /// Blocking version of new
    pub fn new(size: f64) -> Self {
        pollster::block_on(Self::new_async(size))
    }
    
    /// Add a void (absolute nothing) at position with radius
    pub fn add_void(&mut self, x: f64, y: f64, z: f64, radius: f64) {
        self.voids.push((x, y, z, radius));
    }
    
    /// Check if line-of-sight between two points is blocked by any void
    pub fn los_blocked(&self, p1: (f64, f64, f64), p2: (f64, f64, f64)) -> bool {
        // Check each void for intersection with line segment p1-p2
        for &(vx, vy, vz, radius) in &self.voids {
            // Ray-sphere intersection test
            // Line from p1 to p2: P = p1 + t*(p2-p1), t in [0,1]
            let dx = p2.0 - p1.0;
            let dy = p2.1 - p1.1;
            let dz = p2.2 - p1.2;
            
            let fx = p1.0 - vx;
            let fy = p1.1 - vy;
            let fz = p1.2 - vz;
            
            let a = dx*dx + dy*dy + dz*dz;
            let b = 2.0 * (fx*dx + fy*dy + fz*dz);
            let c = fx*fx + fy*fy + fz*fz - radius*radius;
            
            let discriminant = b*b - 4.0*a*c;
            
            if discriminant >= 0.0 && a > 0.0 {
                let sqrt_d = discriminant.sqrt();
                let t1 = (-b - sqrt_d) / (2.0 * a);
                let t2 = (-b + sqrt_d) / (2.0 * a);
                
                // Check if intersection is within the segment [0,1]
                if (t1 >= 0.0 && t1 <= 1.0) || (t2 >= 0.0 && t2 <= 1.0) || (t1 < 0.0 && t2 > 1.0) {
                    return true;  // Bond blocked by void!
                }
            }
        }
        false
    }
    
    /// Spawn a vortex
    pub fn spawn_vortex(&mut self, mut vortex: Vortex) -> u64 {
        let id = self.next_id;
        self.next_id += 1;
        vortex.id = id;
        self.vortices.insert(id, vortex);
        id
    }
    
    /// Run one tick using GPU for force calculations
    pub fn tick(&mut self, dt: f64) {
        self.tick_count += 1;
        
        let n = self.vortices.len();
        if n < 2 {
            return;
        }
        
        // Phase 1: Advance phases (CPU - simple)
        for vortex in self.vortices.values_mut() {
            vortex.tick_phase(dt);
        }
        
        // Phase 2: GPU-accelerated force calculation (3D)
        let accelerations = self.compute_forces_gpu(dt as f32);
        
        // Apply 3D accelerations
        let ids: Vec<u64> = self.vortices.keys().cloned().collect();
        for (i, id) in ids.iter().enumerate() {
            if let Some(vortex) = self.vortices.get_mut(id) {
                vortex.vx += accelerations[i].0 * dt;
                vortex.vy += accelerations[i].1 * dt;
                vortex.vz += accelerations[i].2 * dt;
            }
        }
        
        // Phase 3: Update positions
        for vortex in self.vortices.values_mut() {
            vortex.tick_position(dt);
        }
        
        // Phase 4: Process interactions (CPU - needs complex logic)
        self.process_interactions();
        
        // Phase 5: Boundary handling (3D sphere)
        for vortex in self.vortices.values_mut() {
            let dist = (vortex.x * vortex.x + vortex.y * vortex.y + vortex.z * vortex.z).sqrt();
            if dist > self.size * 0.9 {
                let nx = -vortex.x / dist;
                let ny = -vortex.y / dist;
                let nz = -vortex.z / dist;
                vortex.vx += nx * 0.1;
                vortex.vy += ny * 0.1;
                vortex.vz += nz * 0.1;
                vortex.vx *= 0.9;
                vortex.vy *= 0.9;
                vortex.vz *= 0.9;
            }
        }
    }
    
    /// Compute forces on GPU - returns (ax, ay, az) for each vortex
    fn compute_forces_gpu(&self, dt: f32) -> Vec<(f64, f64, f64)> {
        let n = self.vortices.len();
        
        // Convert to GPU format (3D)
        let ids: Vec<u64> = self.vortices.keys().cloned().collect();
        let gpu_vortices: Vec<GpuVortex> = ids.iter()
            .map(|id| {
                let v = &self.vortices[id];
                GpuVortex {
                    x: v.x as f32,
                    y: v.y as f32,
                    z: v.z as f32,
                    vx: v.vx as f32,
                    vy: v.vy as f32,
                    vz: v.vz as f32,
                    energy: v.energy as f32,
                    frequency: v.frequency as f32,
                    phase: v.phase as f32,
                    spin: v.spin as f32,
                    _pad1: 0.0,
                    _pad2: 0.0,
                }
            })
            .collect();
        
        // Create GPU buffers
        let params = SimParams {
            dt,
            num_vortices: n as u32,
            field_size: self.size as f32,
            _pad: 0.0,
        };
        
        let params_buffer = self.device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
            label: Some("Params Buffer"),
            contents: bytemuck::cast_slice(&[params]),
            usage: wgpu::BufferUsages::UNIFORM,
        });
        
        let vortex_buffer = self.device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
            label: Some("Vortex Buffer"),
            contents: bytemuck::cast_slice(&gpu_vortices),
            usage: wgpu::BufferUsages::STORAGE,
        });
        
        let accel_size = (n * std::mem::size_of::<GpuAccel>()) as u64;
        let accel_buffer = self.device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("Accel Buffer"),
            size: accel_size,
            usage: wgpu::BufferUsages::STORAGE | wgpu::BufferUsages::COPY_SRC,
            mapped_at_creation: false,
        });
        
        let staging_buffer = self.device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("Staging Buffer"),
            size: accel_size,
            usage: wgpu::BufferUsages::MAP_READ | wgpu::BufferUsages::COPY_DST,
            mapped_at_creation: false,
        });
        
        // Create bind group
        let bind_group = self.device.create_bind_group(&wgpu::BindGroupDescriptor {
            label: Some("Vortex Bind Group"),
            layout: &self.bind_group_layout,
            entries: &[
                wgpu::BindGroupEntry {
                    binding: 0,
                    resource: params_buffer.as_entire_binding(),
                },
                wgpu::BindGroupEntry {
                    binding: 1,
                    resource: vortex_buffer.as_entire_binding(),
                },
                wgpu::BindGroupEntry {
                    binding: 2,
                    resource: accel_buffer.as_entire_binding(),
                },
            ],
        });
        
        // Dispatch compute shader
        let mut encoder = self.device.create_command_encoder(&wgpu::CommandEncoderDescriptor {
            label: Some("Compute Encoder"),
        });
        
        {
            let mut pass = encoder.begin_compute_pass(&wgpu::ComputePassDescriptor {
                label: Some("Force Compute Pass"),
                timestamp_writes: None,
            });
            pass.set_pipeline(&self.compute_pipeline);
            pass.set_bind_group(0, &bind_group, &[]);
            pass.dispatch_workgroups((n as u32 + 255) / 256, 1, 1);
        }
        
        encoder.copy_buffer_to_buffer(&accel_buffer, 0, &staging_buffer, 0, accel_size);
        
        self.queue.submit(std::iter::once(encoder.finish()));
        
        // Read back results
        let buffer_slice = staging_buffer.slice(..);
        let (sender, receiver) = std::sync::mpsc::channel();
        buffer_slice.map_async(wgpu::MapMode::Read, move |result| {
            sender.send(result).unwrap();
        });
        
        self.device.poll(wgpu::Maintain::Wait);
        receiver.recv().unwrap().unwrap();
        
        let data = buffer_slice.get_mapped_range();
        let accels: &[GpuAccel] = bytemuck::cast_slice(&data);
        
        let result: Vec<(f64, f64, f64)> = accels.iter()
            .map(|a| (a.ax as f64, a.ay as f64, a.az as f64))
            .collect();
        
        drop(data);
        staging_buffer.unmap();
        
        result
    }
    
    /// Process interactions (merges and orbits) - CPU side
    fn process_interactions(&mut self) {
        use rayon::prelude::*;
        use std::collections::HashSet;
        
        let ids: Vec<u64> = self.vortices.keys().cloned().collect();
        let n = ids.len();
        
        let vortex_data: Vec<(u64, Vortex)> = self.vortices.iter()
            .map(|(id, v)| (*id, v.clone()))
            .collect();
        
        let orbit_threshold = self.orbit_threshold;
        
        let existing_orbits: HashSet<(u64, u64)> = self.orbital_pairs.iter()
            .map(|p| if p.a < p.b { (p.a, p.b) } else { (p.b, p.a) })
            .collect();
        
        // Count current bonds per vortex for max-13 limit
        let mut bond_count: HashMap<u64, usize> = HashMap::new();
        for pair in &self.orbital_pairs {
            *bond_count.entry(pair.a).or_insert(0) += 1;
            *bond_count.entry(pair.b).or_insert(0) += 1;
        }
        
        // Get void data for LOS checking (clone for parallel access)
        let voids = self.voids.clone();
        
        // Parallel interaction detection
        let interactions: Vec<(u64, u64, bool)> = (0..n).into_par_iter()
            .flat_map(|i| {
                let mut local = Vec::new();
                let (id_a, ref va) = vortex_data[i];
                
                for j in (i + 1)..n {
                    let (id_b, ref vb) = vortex_data[j];
                    
                    match determine_interaction(va, vb, orbit_threshold) {
                        Interaction::Merge => {
                            local.push((id_a, id_b, true));
                        }
                        Interaction::Orbit => {
                            let key = if id_a < id_b { (id_a, id_b) } else { (id_b, id_a) };
                            if !existing_orbits.contains(&key) {
                                // Check LOS - is there a void between these vortices?
                                let p1 = (va.x, va.y, va.z);
                                let p2 = (vb.x, vb.y, vb.z);
                                let mut blocked = false;
                                
                                for &(vx, vy, vz, radius) in &voids {
                                    // Ray-sphere intersection
                                    let dx = p2.0 - p1.0;
                                    let dy = p2.1 - p1.1;
                                    let dz = p2.2 - p1.2;
                                    let fx = p1.0 - vx;
                                    let fy = p1.1 - vy;
                                    let fz = p1.2 - vz;
                                    let a = dx*dx + dy*dy + dz*dz;
                                    let b = 2.0 * (fx*dx + fy*dy + fz*dz);
                                    let c = fx*fx + fy*fy + fz*fz - radius*radius;
                                    let discriminant = b*b - 4.0*a*c;
                                    
                                    if discriminant >= 0.0 && a > 0.0 {
                                        let sqrt_d = discriminant.sqrt();
                                        let t1 = (-b - sqrt_d) / (2.0 * a);
                                        let t2 = (-b + sqrt_d) / (2.0 * a);
                                        if (t1 >= 0.0 && t1 <= 1.0) || (t2 >= 0.0 && t2 <= 1.0) || (t1 < 0.0 && t2 > 1.0) {
                                            blocked = true;
                                            break;
                                        }
                                    }
                                }
                                
                                if !blocked {
                                    local.push((id_a, id_b, false));
                                }
                            }
                        }
                        _ => {}
                    }
                }
                local
            })
            .collect();
        
        // Process results
        let mut merges = Vec::new();
        let mut new_orbits = Vec::new();
        let mut merged_ids = HashSet::new();
        let mut void_blocked = 0u64;
        
        for (id_a, id_b, is_merge) in interactions {
            if is_merge {
                if !merged_ids.contains(&id_a) && !merged_ids.contains(&id_b) {
                    merges.push((id_a, id_b));
                    merged_ids.insert(id_a);
                    merged_ids.insert(id_b);
                }
            } else {
                // Check max 13 bonds (kissing number) - no singularities!
                let count_a = bond_count.get(&id_a).copied().unwrap_or(0);
                let count_b = bond_count.get(&id_b).copied().unwrap_or(0);
                
                if count_a < MAX_BONDS_PER_VORTEX && count_b < MAX_BONDS_PER_VORTEX {
                    new_orbits.push((id_a, id_b));
                    // Update counts for subsequent checks this tick
                    *bond_count.entry(id_a).or_insert(0) += 1;
                    *bond_count.entry(id_b).or_insert(0) += 1;
                }
            }
        }
        
        self.bonds_blocked_by_void += void_blocked;
        
        // Execute merges
        for (id_a, id_b) in merges {
            if let (Some(a), Some(b)) = (self.vortices.get(&id_a), self.vortices.get(&id_b)) {
                let merged = merge_vortices(a, b, self.next_id);
                self.next_id += 1;
                self.vortices.remove(&id_a);
                self.vortices.remove(&id_b);
                self.vortices.insert(merged.id, merged);
                self.total_merges += 1;
                self.orbital_pairs.retain(|p| p.a != id_a && p.a != id_b && p.b != id_a && p.b != id_b);
            }
        }
        
        // Create new orbits
        for (id_a, id_b) in new_orbits {
            let pair = OrbitalPair {
                id: self.next_id,
                a: id_a,
                b: id_b,
                orbital_phase: 0.0,
                stability: 0.5,
            };
            self.next_id += 1;
            self.orbital_pairs.push(pair);
            self.total_orbits_formed += 1;
        }
    }
    
    pub fn total_energy(&self) -> f64 {
        self.vortices.values().map(|v| v.energy).sum()
    }
    
    pub fn vortex_count(&self) -> usize {
        self.vortices.len()
    }
    
    pub fn report(&self) -> String {
        let mut s = format!(
            "GPU Field State (tick {})\n  Vortices: {}\n  Orbital pairs: {}\n  Total energy: {:.2}\n  Merges: {}\n  Orbits formed: {}",
            self.tick_count,
            self.vortices.len(),
            self.orbital_pairs.len(),
            self.total_energy(),
            self.total_merges,
            self.total_orbits_formed
        );
        
        if !self.voids.is_empty() {
            s.push_str(&format!("\n  Voids: {} (bonds blocked: {})", self.voids.len(), self.bonds_blocked_by_void));
        }
        
        // Report bond distribution
        let mut bond_counts: HashMap<usize, usize> = HashMap::new();
        let mut per_vortex: HashMap<u64, usize> = HashMap::new();
        for pair in &self.orbital_pairs {
            *per_vortex.entry(pair.a).or_insert(0) += 1;
            *per_vortex.entry(pair.b).or_insert(0) += 1;
        }
        for &count in per_vortex.values() {
            *bond_counts.entry(count).or_insert(0) += 1;
        }
        
        if let Some(&max_bonds) = per_vortex.values().max() {
            let at_limit = per_vortex.values().filter(|&&c| c >= 13).count();
            s.push_str(&format!("\n  Max bonds/vortex: {} (at limit 13: {})", max_bonds, at_limit));
        }
        
        s
    }
}

/// WGSL compute shader for 3D force calculation
const COMPUTE_SHADER: &str = r#"
struct Params {
    dt: f32,
    num_vortices: u32,
    field_size: f32,
    _pad: f32,
}

struct Vortex {
    x: f32,
    y: f32,
    z: f32,
    vx: f32,
    vy: f32,
    vz: f32,
    energy: f32,
    frequency: f32,
    phase: f32,
    spin: f32,
    _pad1: f32,
    _pad2: f32,
}

struct Accel {
    ax: f32,
    ay: f32,
    az: f32,
    _pad: f32,
}

@group(0) @binding(0) var<uniform> params: Params;
@group(0) @binding(1) var<storage, read> vortices: array<Vortex>;
@group(0) @binding(2) var<storage, read_write> accels: array<Accel>;

@compute @workgroup_size(256)
fn main(@builtin(global_invocation_id) global_id: vec3<u32>) {
    let idx = global_id.x;
    if (idx >= params.num_vortices) {
        return;
    }
    
    let me = vortices[idx];
    var ax: f32 = 0.0;
    var ay: f32 = 0.0;
    var az: f32 = 0.0;
    
    // Sum forces from all other vortices (3D)
    for (var j: u32 = 0u; j < params.num_vortices; j = j + 1u) {
        if (j == idx) {
            continue;
        }
        
        let other = vortices[j];
        let dx = other.x - me.x;
        let dy = other.y - me.y;
        let dz = other.z - me.z;
        let dist_sq = dx * dx + dy * dy + dz * dz;
        
        if (dist_sq < 0.000001) {
            continue;
        }
        
        let dist = sqrt(dist_sq);
        
        // F = E1 * E2 / r² (inverse square law works in any dimension)
        let force = (me.energy * other.energy) / dist_sq;
        let accel = force / me.energy;
        
        // Accumulate 3D acceleration (normalized direction)
        ax = ax + (dx / dist) * accel;
        ay = ay + (dy / dist) * accel;
        az = az + (dz / dist) * accel;
    }
    
    accels[idx].ax = ax;
    accels[idx].ay = ay;
    accels[idx].az = az;
}
"#;
