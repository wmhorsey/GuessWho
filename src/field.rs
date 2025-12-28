// field.rs - The space where vortices dance
//
// The field is pure potential. Empty space that energy can move through.
// But λ=0 (true nothing) is unstable - energy avoids the void.
//
// The field doesn't push. It only allows attraction to happen.

use std::collections::HashMap;
use rayon::prelude::*;
use crate::vortex::{Vortex, Interaction, OrbitalPair, determine_interaction, merge_vortices};

/// The field - space where vortices exist and interact
pub struct Field {
    pub vortices: HashMap<u64, Vortex>,
    pub orbital_pairs: Vec<OrbitalPair>,
    pub next_id: u64,
    
    // Field properties
    pub size: f64,              // radius of the field
    pub orbit_threshold: f64,   // distance at which orbits can form
    
    // Statistics
    pub tick_count: u64,
    pub total_merges: u64,
    pub total_orbits_formed: u64,
}

impl Field {
    pub fn new(size: f64) -> Self {
        Self {
            vortices: HashMap::new(),
            orbital_pairs: Vec::new(),
            next_id: 1,
            size,
            orbit_threshold: size / 20.0,  // relative to field size
            tick_count: 0,
            total_merges: 0,
            total_orbits_formed: 0,
        }
    }
    
    /// Spawn a new vortex in the field
    pub fn spawn(&mut self, x: f64, y: f64, energy: f64) -> u64 {
        let id = self.next_id;
        self.next_id += 1;
        
        let vortex = Vortex::new(id, x, y, energy);
        self.vortices.insert(id, vortex);
        id
    }
    
    /// Spawn with full configuration
    pub fn spawn_vortex(&mut self, mut vortex: Vortex) -> u64 {
        let id = self.next_id;
        self.next_id += 1;
        vortex.id = id;
        self.vortices.insert(id, vortex);
        id
    }
    
    /// Run one tick of the simulation
    pub fn tick(&mut self, dt: f64) {
        self.tick_count += 1;
        
        // Phase 1: All vortices advance their phase (they're always spinning)
        for vortex in self.vortices.values_mut() {
            vortex.tick_phase(dt);
        }
        
        // Phase 2: Compute attractions in PARALLEL
        // Gather all vortex data into a Vec for parallel processing
        let vortex_data: Vec<(u64, f64, f64, f64)> = self.vortices.values()
            .map(|v| (v.id, v.x, v.y, v.energy))
            .collect();
        
        // Parallel computation of accelerations
        let accelerations: Vec<(u64, f64, f64)> = vortex_data.par_iter()
            .map(|(id, x, y, energy)| {
                let mut accel_x = 0.0;
                let mut accel_y = 0.0;
                
                for (other_id, ox, oy, oe) in &vortex_data {
                    if id == other_id {
                        continue;
                    }
                    
                    let dx = ox - x;
                    let dy = oy - y;
                    let dist_sq = dx * dx + dy * dy;
                    
                    if dist_sq < 0.000001 {
                        continue;
                    }
                    
                    let dist = dist_sq.sqrt();
                    
                    // F = E1 * E2 / r²
                    let force = (energy * oe) / dist_sq;
                    let accel = force / energy;
                    
                    // Direction (normalized)
                    accel_x += (dx / dist) * accel;
                    accel_y += (dy / dist) * accel;
                }
                
                (*id, accel_x, accel_y)
            })
            .collect();
        
        // Apply accelerations (sequential, but fast)
        for (id, ax, ay) in accelerations {
            if let Some(vortex) = self.vortices.get_mut(&id) {
                vortex.vx += ax * dt;
                vortex.vy += ay * dt;
            }
        }
        
        // Phase 3: Update positions
        for vortex in self.vortices.values_mut() {
            vortex.tick_position(dt);
        }
        
        // Phase 4: Check for interactions (merges, new orbits)
        self.process_interactions();
        
        // Phase 5: Handle boundary (soft boundary - things slow down at edge)
        for vortex in self.vortices.values_mut() {
            let dist_from_center = (vortex.x * vortex.x + vortex.y * vortex.y).sqrt();
            if dist_from_center > self.size * 0.9 {
                // Soft reflection - reduce velocity toward center
                let factor = 0.9;
                let nx = -vortex.x / dist_from_center;
                let ny = -vortex.y / dist_from_center;
                
                // Add inward velocity component
                vortex.vx += nx * 0.1;
                vortex.vy += ny * 0.1;
                
                // Dampen outward velocity
                vortex.vx *= factor;
                vortex.vy *= factor;
            }
        }
    }
    
    /// Process merges and orbit formation - PARALLEL VERSION
    fn process_interactions(&mut self) {
        let ids: Vec<u64> = self.vortices.keys().cloned().collect();
        let n = ids.len();
        
        // Collect vortex data for parallel processing
        let vortex_data: Vec<(u64, Vortex)> = self.vortices.iter()
            .map(|(id, v)| (*id, v.clone()))
            .collect();
        
        let orbit_threshold = self.orbit_threshold;
        
        // Existing orbital pairs as a set for fast lookup
        let existing_orbits: std::collections::HashSet<(u64, u64)> = self.orbital_pairs.iter()
            .map(|p| if p.a < p.b { (p.a, p.b) } else { (p.b, p.a) })
            .collect();
        
        // Parallel search for interactions
        // We split the work into chunks of pairs
        let interactions: Vec<(u64, u64, bool)> = (0..n).into_par_iter()
            .flat_map(|i| {
                let mut local_interactions = Vec::new();
                let (id_a, ref vortex_a) = vortex_data[i];
                
                for j in (i + 1)..n {
                    let (id_b, ref vortex_b) = vortex_data[j];
                    
                    match determine_interaction(vortex_a, vortex_b, orbit_threshold) {
                        Interaction::Merge => {
                            local_interactions.push((id_a, id_b, true)); // true = merge
                        }
                        Interaction::Orbit => {
                            let key = if id_a < id_b { (id_a, id_b) } else { (id_b, id_a) };
                            if !existing_orbits.contains(&key) {
                                local_interactions.push((id_a, id_b, false)); // false = orbit
                            }
                        }
                        _ => {}
                    }
                }
                
                local_interactions
            })
            .collect();
        
        // Separate merges and orbits (sequential, but small)
        let mut merges: Vec<(u64, u64)> = Vec::new();
        let mut new_orbits: Vec<(u64, u64)> = Vec::new();
        let mut merged_ids: std::collections::HashSet<u64> = std::collections::HashSet::new();
        
        for (id_a, id_b, is_merge) in interactions {
            if is_merge {
                // Skip if either already merged
                if !merged_ids.contains(&id_a) && !merged_ids.contains(&id_b) {
                    merges.push((id_a, id_b));
                    merged_ids.insert(id_a);
                    merged_ids.insert(id_b);
                }
            } else {
                new_orbits.push((id_a, id_b));
            }
        }
        
        // Execute merges
        for (id_a, id_b) in merges {
            if let (Some(a), Some(b)) = (self.vortices.get(&id_a), self.vortices.get(&id_b)) {
                let merged = merge_vortices(a, b, self.next_id);
                self.next_id += 1;
                
                self.vortices.remove(&id_a);
                self.vortices.remove(&id_b);
                self.vortices.insert(merged.id, merged);
                
                self.total_merges += 1;
                
                // Remove any orbital pairs that involved merged vortices
                self.orbital_pairs.retain(|p| p.a != id_a && p.a != id_b && p.b != id_a && p.b != id_b);
            }
        }
        
        // Create new orbital pairs
        for (id_a, id_b) in new_orbits {
            let pair = OrbitalPair {
                id: self.next_id,
                a: id_a,
                b: id_b,
                orbital_phase: 0.0,
                stability: 0.5,  // start at medium stability
            };
            self.next_id += 1;
            self.orbital_pairs.push(pair);
            self.total_orbits_formed += 1;
        }
    }
    
    /// Get total energy in the field (should be conserved)
    pub fn total_energy(&self) -> f64 {
        self.vortices.values().map(|v| v.energy).sum()
    }
    
    /// Get count of vortices
    pub fn vortex_count(&self) -> usize {
        self.vortices.len()
    }
    
    /// Print current state
    pub fn report(&self) -> String {
        let mut s = String::new();
        s.push_str(&format!("Field State (tick {})\n", self.tick_count));
        s.push_str(&format!("  Vortices: {}\n", self.vortices.len()));
        s.push_str(&format!("  Orbital pairs: {}\n", self.orbital_pairs.len()));
        s.push_str(&format!("  Total energy: {:.2}\n", self.total_energy()));
        s.push_str(&format!("  Merges: {}\n", self.total_merges));
        s.push_str(&format!("  Orbits formed: {}\n", self.total_orbits_formed));
        s
    }
    
    /// ASCII visualization of the field
    pub fn visualize(&self, width: usize, height: usize) -> String {
        let mut grid = vec![vec![' '; width]; height];
        
        // Draw boundary
        for i in 0..width {
            grid[0][i] = '.';
            grid[height - 1][i] = '.';
        }
        for j in 0..height {
            grid[j][0] = '.';
            grid[j][width - 1] = '.';
        }
        
        // Draw vortices
        for vortex in self.vortices.values() {
            // Map position to grid
            let gx = ((vortex.x / self.size + 1.0) * 0.5 * (width - 2) as f64) as usize + 1;
            let gy = ((vortex.y / self.size + 1.0) * 0.5 * (height - 2) as f64) as usize + 1;
            
            if gx > 0 && gx < width - 1 && gy > 0 && gy < height - 1 {
                // Symbol based on energy
                let symbol = if vortex.energy > 10.0 {
                    '@'
                } else if vortex.energy > 5.0 {
                    'O'
                } else if vortex.energy > 2.0 {
                    'o'
                } else if vortex.energy > 1.0 {
                    '*'
                } else {
                    '·'
                };
                
                // Spin indicator
                let symbol = if vortex.spin > 0 { symbol } else { symbol.to_ascii_lowercase() };
                
                grid[gy][gx] = symbol;
            }
        }
        
        // Draw center
        let cx = width / 2;
        let cy = height / 2;
        if grid[cy][cx] == ' ' {
            grid[cy][cx] = '+';
        }
        
        // Convert to string
        grid.iter()
            .map(|row| row.iter().collect::<String>())
            .collect::<Vec<_>>()
            .join("\n")
    }
}
