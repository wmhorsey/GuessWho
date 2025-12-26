//! Packet — Light/Energy units carrying information
//!
//! At creation each packet is assigned a type and a memory reference,
//! then freed into the system. The packet doesn't know where it's going —
//! it just resonates, and the field sorts it.
//!
//! Packets travel through 3D space with their resonance signature
//! determining which nodes attract or repel them.
//!
//! KEY INSIGHT: Resonance is derived FROM content at WORD level.
//! This enables LANGUAGE-AGNOSTIC clustering — "hello", "hola", "bonjour"
//! cluster by meaning, not by character spelling patterns.

use crate::resonance::{Resonance, Velocity};
use crate::semantic::ContentClass;
use serde::{Serialize, Deserialize};

/// Unique identifier for a packet in the system
pub type PacketId = u64;

/// Memory reference — where this packet "snaps back" to
pub type MemoryRef = u64;

/// The classification of packet at birth
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum PacketType {
    /// Active signal — something happening now
    Signal,
    /// Memory trace — something that happened
    Memory,
    /// Query — something seeking information
    Query,
    /// Echo — reflection of another packet
    Echo,
    /// Anomaly — defies categorization (the 0.01%)
    Anomaly,
    /// Cosmic ray — hit the CMB, came back changed
    CosmicRay,
}

impl PacketType {
    /// Get the base resonance wavelength for this packet type
    pub fn base_wavelength(&self) -> f64 {
        match self {
            PacketType::Signal => 0.1,
            PacketType::Memory => 0.3,
            PacketType::Query => 0.5,
            PacketType::Echo => 0.7,
            PacketType::Anomaly => rand_wavelength(), // Anomalies are unpredictable
            PacketType::CosmicRay => rand_wavelength(), // Changed by the CMB
        }
    }

    /// Get base amplitude for this type
    pub fn base_amplitude(&self) -> f64 {
        match self {
            PacketType::Signal => 0.9,   // Strong and clear
            PacketType::Memory => 0.6,   // Fading over time
            PacketType::Query => 0.8,    // Active search
            PacketType::Echo => 0.4,     // Weakened reflection
            PacketType::Anomaly => rand_amplitude(),
            PacketType::CosmicRay => 1.0, // Energized by the CMB
        }
    }

    pub fn name(&self) -> &'static str {
        match self {
            PacketType::Signal => "Signal",
            PacketType::Memory => "Memory",
            PacketType::Query => "Query",
            PacketType::Echo => "Echo",
            PacketType::Anomaly => "Anomaly",
            PacketType::CosmicRay => "CosmicRay",
        }
    }

    pub fn symbol(&self) -> char {
        match self {
            PacketType::Signal => '!',
            PacketType::Memory => '#',
            PacketType::Query => '?',
            PacketType::Echo => '~',
            PacketType::Anomaly => '*',
            PacketType::CosmicRay => '@',
        }
    }
}

/// State of a packet in its lifecycle
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub enum PacketState {
    /// Floating in the field, seeking a node
    InFlight,
    /// Captured by a node, being held
    Captured,
    /// Consumed — snapped back to switchboard for lookup
    Consumed,
    /// Marked for release (garbage collection candidate)
    Released,
    /// Hit the CMB boundary
    AtBoundary,
}

/// 3D position in the universe
#[derive(Clone, Debug, Default, Serialize, Deserialize)]
pub struct Position {
    pub x: f64,
    pub y: f64,
    pub z: f64,
}

impl Position {
    pub fn new(x: f64, y: f64, z: f64) -> Self {
        Self { x, y, z }
    }

    pub fn random_in_sphere(radius: f64) -> Self {
        use rand::Rng;
        let mut rng = rand::thread_rng();
        // Random point inside a sphere
        let r = radius * rng.gen::<f64>().cbrt(); // Cube root for uniform distribution
        let theta = rng.gen_range(0.0..2.0 * std::f64::consts::PI);
        let phi = rng.gen_range(0.0..std::f64::consts::PI);
        Self {
            x: r * phi.sin() * theta.cos(),
            y: r * phi.sin() * theta.sin(),
            z: r * phi.cos(),
        }
    }

    pub fn distance_from_origin(&self) -> f64 {
        (self.x * self.x + self.y * self.y + self.z * self.z).sqrt()
    }

    pub fn distance_to(&self, other: &Position) -> f64 {
        let dx = self.x - other.x;
        let dy = self.y - other.y;
        let dz = self.z - other.z;
        (dx * dx + dy * dy + dz * dz).sqrt()
    }

    /// Get 2D projection for visualization (looking down Z axis)
    pub fn project_2d(&self) -> (f64, f64) {
        (self.x, self.y)
    }
}

/// A packet of light/energy carrying information through resonance
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct Packet {
    /// Unique identifier
    pub id: PacketId,
    /// The packet's resonance signature (wavelength, amplitude, velocity)
    pub resonance: Resonance,
    /// Classification at birth
    pub packet_type: PacketType,
    /// Memory reference for snap-back
    pub memory_ref: MemoryRef,
    /// Current lifecycle state
    pub state: PacketState,
    /// 3D position in the universe
    pub position: Position,
    /// Number of times this packet has been reused
    pub reuse_count: u32,
    /// Tick when this packet was created
    pub born_at: u64,
    /// Distance traveled (odometer)
    pub distance_traveled: f64,
    /// Times this packet has hit the CMB
    pub cmb_hits: u32,
    /// Optional content that generated this packet's resonance
    /// When present, resonance is derived from content, not random
    #[serde(default)]
    pub content: Option<String>,
}

impl Packet {
    /// Create a new packet with the given type
    pub fn new(id: PacketId, packet_type: PacketType, memory_ref: MemoryRef, universe_radius: f64) -> Self {
        use rand::Rng;
        let mut rng = rand::thread_rng();
        
        let base_wl = packet_type.base_wavelength();
        let base_amp = packet_type.base_amplitude();
        
        // Add some variation so not all packets of same type are identical
        let wl_variation = rng.gen_range(-0.025..0.025);
        let amp_variation = rng.gen_range(-0.1..0.1);
        
        let velocity = Velocity::random();
        let resonance = Resonance::full(
            base_wl + wl_variation,
            (base_amp + amp_variation).clamp(0.1, 1.0),
            velocity,
        );
        
        Self {
            id,
            resonance,
            packet_type,
            memory_ref,
            state: PacketState::InFlight,
            position: Position::random_in_sphere(universe_radius * 0.8),
            reuse_count: 0,
            born_at: 0,
            distance_traveled: 0.0,
            cmb_hits: 0,
            content: None,
        }
    }

    /// Create a packet from actual content — THE KEY INSIGHT
    /// 
    /// The resonance is DERIVED from the content, not assigned randomly.
    /// This means:
    /// - Similar content → similar resonance → clusters at same nodes
    /// - The "meaning" is encoded in the resonance signature
    /// - Sorting by resonance IS sorting by semantic similarity
    /// 
    /// This is the quantum tunneling mechanism: content "tunnels" to
    /// the right location based on its inherent resonance.
    pub fn from_content(
        id: PacketId,
        content: &str,
        memory_ref: MemoryRef,
        universe_radius: f64,
    ) -> Self {
        // Use WORD-LEVEL resonance for language-agnostic clustering
        // This means "hello", "hola", "bonjour" cluster by meaning, not spelling
        let resonance = crate::semantic::words_to_resonance(content);
        
        // Classify content to determine packet type
        let packet_type = match ContentClass::classify(content) {
            ContentClass::Query => PacketType::Query,
            ContentClass::Command => PacketType::Signal,
            ContentClass::Numeric => PacketType::Memory,
            ContentClass::Code => PacketType::Memory,
            ContentClass::Data => PacketType::Memory,
            ContentClass::Textual => PacketType::Signal,
        };
        
        Self {
            id,
            resonance,
            packet_type,
            memory_ref,
            state: PacketState::InFlight,
            position: Position::random_in_sphere(universe_radius * 0.8),
            reuse_count: 0,
            born_at: 0,
            distance_traveled: 0.0,
            cmb_hits: 0,
            content: Some(content.to_string()),
        }
    }
    
    /// Set resonance for this packet (builder pattern)
    pub fn with_resonance(mut self, resonance: Resonance) -> Self {
        self.resonance = resonance;
        self
    }

    /// Update position based on resonance velocity
    pub fn update_position(&mut self, dt: f64) {
        let old_pos = self.position.clone();
        
        self.position.x += self.resonance.velocity.x * dt;
        self.position.y += self.resonance.velocity.y * dt;
        self.position.z += self.resonance.velocity.z * dt;
        
        self.distance_traveled += old_pos.distance_to(&self.position);
        
        // Drift the resonance slightly over time
        self.resonance.drift(dt * 0.01);
    }

    /// Apply gravitational-like force toward a point
    pub fn apply_force(&mut self, target: &Position, force_magnitude: f64) {
        let dx = target.x - self.position.x;
        let dy = target.y - self.position.y;
        let dz = target.z - self.position.z;
        let dist = (dx * dx + dy * dy + dz * dz).sqrt().max(0.1);
        
        // Inverse square law
        let force = force_magnitude / (dist * dist);
        let norm = dist.max(0.001);
        
        self.resonance.velocity.x += (dx / norm) * force;
        self.resonance.velocity.y += (dy / norm) * force;
        self.resonance.velocity.z += (dz / norm) * force;
    }

    /// Apply velocity damping (space has slight resistance)
    pub fn apply_damping(&mut self, factor: f64) {
        self.resonance.velocity.x *= factor;
        self.resonance.velocity.y *= factor;
        self.resonance.velocity.z *= factor;
    }

    /// Check if packet has reached the CMB boundary
    pub fn check_cmb(&mut self, universe_radius: f64) -> bool {
        if self.position.distance_from_origin() >= universe_radius {
            self.state = PacketState::AtBoundary;
            self.cmb_hits += 1;
            true
        } else {
            false
        }
    }

    /// Bounce off the CMB (reflect velocity, possibly change)
    pub fn bounce_from_cmb(&mut self, universe_radius: f64) {
        // Move back inside
        let dist = self.position.distance_from_origin();
        if dist > 0.0 {
            let scale = (universe_radius * 0.99) / dist;
            self.position.x *= scale;
            self.position.y *= scale;
            self.position.z *= scale;
        }
        
        // Reflect velocity (bounce off the CMB)
        // Calculate normal at boundary (points toward center)
        let nx = -self.position.x / dist;
        let ny = -self.position.y / dist;
        let nz = -self.position.z / dist;
        
        // Reflect: v' = v - 2(v·n)n
        let dot = self.resonance.velocity.x * nx 
                + self.resonance.velocity.y * ny 
                + self.resonance.velocity.z * nz;
        
        self.resonance.velocity.x -= 2.0 * dot * nx;
        self.resonance.velocity.y -= 2.0 * dot * ny;
        self.resonance.velocity.z -= 2.0 * dot * nz;
        
        // CMB interaction can boost amplitude (cosmic rays are energized)
        self.resonance.boost(0.2);
        
        self.state = PacketState::InFlight;
    }

    /// Mark as captured by a node
    pub fn capture(&mut self) {
        self.state = PacketState::Captured;
    }

    /// Consume the packet — triggers snap-back to switchboard
    pub fn consume(&mut self) {
        self.state = PacketState::Consumed;
    }

    /// Mark for release (garbage collection)
    pub fn release(&mut self) {
        self.state = PacketState::Released;
    }

    /// Reincarnate the packet with a new memory reference
    pub fn reincarnate(&mut self, new_memory_ref: MemoryRef, universe_radius: f64) {
        self.memory_ref = new_memory_ref;
        self.state = PacketState::InFlight;
        self.position = Position::random_in_sphere(universe_radius * 0.8);
        self.resonance.velocity = Velocity::random();
        self.resonance.boost(0.3); // Fresh energy
        self.reuse_count += 1;
    }

    /// Get a display signature
    pub fn signature(&self) -> String {
        format!(
            "P{}[{}:{}]",
            self.id,
            self.packet_type.name(),
            self.resonance.signature()
        )
    }

    /// Get detailed signature
    pub fn signature_full(&self) -> String {
        format!(
            "P{}[{}] {} pos({:.1},{:.1},{:.1}) dist:{:.1} cmb:{}",
            self.id,
            self.packet_type.name(),
            self.resonance.signature_full(),
            self.position.x,
            self.position.y,
            self.position.z,
            self.distance_traveled,
            self.cmb_hits
        )
    }
}

// Helpers for randomness
fn rand_wavelength() -> f64 {
    use rand::Rng;
    rand::thread_rng().gen::<f64>()
}

fn rand_amplitude() -> f64 {
    use rand::Rng;
    rand::thread_rng().gen_range(0.3..1.0)
}
