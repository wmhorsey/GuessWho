// vortex.rs - The primitives of energy
// 
// What can energy do?
//   - Exist (not nothing)
//   - Attract (toward other energy)
//   - Flow (has direction)
//   - Spin (flow in a loop = stable pattern)
//   - Phase (where in the cycle)
//
// What emerges?
//   - Orbit (missed merge, now bound)
//   - Merge (same phase, direct approach)
//   - Shed (three-body instability ejects one)

use std::f64::consts::PI;

/// A vortex is localized spinning energy.
/// Not "stuff" - a pattern of flow.
#[derive(Clone, Debug)]
pub struct Vortex {
    pub id: u64,
    
    // Position in the field (3D)
    pub x: f64,
    pub y: f64,
    pub z: f64,
    
    // Motion (3D)
    pub vx: f64,
    pub vy: f64,
    pub vz: f64,
    
    // The resonance signature
    pub energy: f64,      // amplitude - how much (determines attraction strength)
    pub frequency: f64,   // how fast the cycle (wavelength = 1/frequency)
    pub phase: f64,       // where in cycle (0 to 2π)
    pub spin: i8,         // +1 clockwise, -1 counter-clockwise
}

/// What can happen when two vortices interact
#[derive(Debug, Clone, PartialEq)]
pub enum Interaction {
    None,           // Too far, no significant interaction
    Attract,        // Moving toward each other
    Orbit,          // Stable binding, circling
    Merge,          // Becoming one larger vortex
}

/// An orbital pair - two vortices bound together
#[derive(Clone, Debug)]
pub struct OrbitalPair {
    pub id: u64,
    pub a: u64,  // vortex id
    pub b: u64,  // vortex id
    pub orbital_phase: f64,  // where in the orbit
    pub stability: f64,      // how stable (1.0 = locked, 0.0 = about to break)
}

/// A stable cluster - multiple vortices in stable configuration
#[derive(Clone, Debug)]
pub struct Cluster {
    pub id: u64,
    pub members: Vec<u64>,
    pub center_x: f64,
    pub center_y: f64,
    pub total_energy: f64,
    pub resonance_signature: ResonanceSignature,
}

/// The combined resonance of a cluster
#[derive(Clone, Debug)]
pub struct ResonanceSignature {
    pub fundamental: f64,     // the dominant frequency
    pub harmonics: Vec<f64>,  // overtones present
    pub net_spin: i8,         // overall rotation direction
}

impl Vortex {
    pub fn new(id: u64, x: f64, y: f64, energy: f64) -> Self {
        Self {
            id,
            x,
            y,
            z: 0.0,
            vx: 0.0,
            vy: 0.0,
            vz: 0.0,
            energy,
            frequency: 1.0,
            phase: 0.0,
            spin: 1,
        }
    }
    
    pub fn new_3d(id: u64, x: f64, y: f64, z: f64, energy: f64) -> Self {
        Self {
            id,
            x,
            y,
            z,
            vx: 0.0,
            vy: 0.0,
            vz: 0.0,
            energy,
            frequency: 1.0,
            phase: 0.0,
            spin: 1,
        }
    }
    
    pub fn with_frequency(mut self, freq: f64) -> Self {
        self.frequency = freq;
        self
    }
    
    pub fn with_phase(mut self, phase: f64) -> Self {
        self.phase = phase % (2.0 * PI);
        self
    }
    
    pub fn with_spin(mut self, spin: i8) -> Self {
        self.spin = if spin >= 0 { 1 } else { -1 };
        self
    }
    
    pub fn with_velocity(mut self, vx: f64, vy: f64) -> Self {
        self.vx = vx;
        self.vy = vy;
        self
    }
    
    pub fn with_velocity_3d(mut self, vx: f64, vy: f64, vz: f64) -> Self {
        self.vx = vx;
        self.vy = vy;
        self.vz = vz;
        self
    }
    
    pub fn with_z(mut self, z: f64) -> Self {
        self.z = z;
        self
    }
    
    /// Distance to another vortex (3D)
    pub fn distance_to(&self, other: &Vortex) -> f64 {
        let dx = other.x - self.x;
        let dy = other.y - self.y;
        let dz = other.z - self.z;
        (dx * dx + dy * dy + dz * dz).sqrt()
    }
    
    /// Attraction force magnitude (energy attracts energy)
    /// Falls off with distance squared (like gravity, like electromagnetism)
    pub fn attraction_to(&self, other: &Vortex) -> f64 {
        let dist = self.distance_to(other);
        if dist < 0.001 {
            return 0.0; // Avoid singularity
        }
        // F = E1 * E2 / r²
        (self.energy * other.energy) / (dist * dist)
    }
    
    /// Phase difference (0 = in sync, π = opposite)
    pub fn phase_difference(&self, other: &Vortex) -> f64 {
        let diff = (self.phase - other.phase).abs();
        if diff > PI { 2.0 * PI - diff } else { diff }
    }
    
    /// Frequency ratio (for detecting harmonics)
    pub fn frequency_ratio(&self, other: &Vortex) -> f64 {
        if other.frequency > self.frequency {
            other.frequency / self.frequency
        } else {
            self.frequency / other.frequency
        }
    }
    
    /// Energy ratio (larger / smaller)
    pub fn energy_ratio(&self, other: &Vortex) -> f64 {
        if other.energy > self.energy {
            other.energy / self.energy
        } else {
            self.energy / other.energy
        }
    }
    
    /// Advance phase by one time step
    pub fn tick_phase(&mut self, dt: f64) {
        self.phase = (self.phase + self.frequency * dt * self.spin as f64) % (2.0 * PI);
        if self.phase < 0.0 {
            self.phase += 2.0 * PI;
        }
    }
    
    /// Update position based on velocity (3D)
    pub fn tick_position(&mut self, dt: f64) {
        self.x += self.vx * dt;
        self.y += self.vy * dt;
        self.z += self.vz * dt;
    }
    
    /// Apply attraction from another vortex (3D)
    pub fn attract_toward(&mut self, other: &Vortex, dt: f64) {
        let dist = self.distance_to(other);
        if dist < 0.001 {
            return;
        }
        
        let force = self.attraction_to(other);
        
        // Direction toward other (3D)
        let dx = (other.x - self.x) / dist;
        let dy = (other.y - self.y) / dist;
        let dz = (other.z - self.z) / dist;
        
        // Acceleration = force / energy (energy acts as inertia)
        let accel = force / self.energy;
        
        self.vx += dx * accel * dt;
        self.vy += dy * accel * dt;
        self.vz += dz * accel * dt;
    }
}

/// Check if a frequency ratio is harmonic (simple integer ratio)
pub fn is_harmonic(ratio: f64) -> bool {
    // Check common harmonic ratios
    let harmonics = [
        1.0,        // unison
        2.0,        // octave
        1.5,        // fifth (3:2)
        1.333,      // fourth (4:3)
        1.25,       // major third (5:4)
        1.2,        // minor third (6:5)
        1.618,      // golden ratio
    ];
    
    for h in harmonics {
        if (ratio - h).abs() < 0.05 {
            return true;
        }
    }
    false
}

/// Determine what interaction will occur between two vortices
pub fn determine_interaction(a: &Vortex, b: &Vortex, orbit_threshold: f64) -> Interaction {
    let dist = a.distance_to(b);
    let phase_diff = a.phase_difference(b);
    let energy_ratio = a.energy_ratio(b);
    let freq_ratio = a.frequency_ratio(b);
    
    // Too far for interaction
    if dist > orbit_threshold * 5.0 {
        return Interaction::None;
    }
    
    // Very close + similar phase + large energy difference = MERGE
    // The larger absorbs the smaller
    if dist < orbit_threshold * 0.1 && phase_diff < PI / 8.0 && energy_ratio > 5.0 {
        return Interaction::Merge;
    }
    
    // Close + harmonic frequencies + phase offset = ORBIT
    // They missed each other, now bound
    if dist < orbit_threshold && is_harmonic(freq_ratio) && phase_diff > PI / 6.0 {
        return Interaction::Orbit;
    }
    
    // Close enough to feel attraction
    if dist < orbit_threshold * 3.0 {
        return Interaction::Attract;
    }
    
    Interaction::None
}

/// Merge two vortices into one (3D)
pub fn merge_vortices(a: &Vortex, b: &Vortex, new_id: u64) -> Vortex {
    let total_energy = a.energy + b.energy;
    
    // Center of energy (weighted by energy) - 3D
    let x = (a.x * a.energy + b.x * b.energy) / total_energy;
    let y = (a.y * a.energy + b.y * b.energy) / total_energy;
    let z = (a.z * a.energy + b.z * b.energy) / total_energy;
    
    // Momentum preserved (weighted by energy since energy = inertia here) - 3D
    let vx = (a.vx * a.energy + b.vx * b.energy) / total_energy;
    let vy = (a.vy * a.energy + b.vy * b.energy) / total_energy;
    let vz = (a.vz * a.energy + b.vz * b.energy) / total_energy;
    
    // Frequency of larger dominates (or average?)
    let frequency = if a.energy > b.energy { a.frequency } else { b.frequency };
    
    // Phase... complex. For now, take from larger.
    let phase = if a.energy > b.energy { a.phase } else { b.phase };
    
    // Spin: if same, keep. If opposite... angular momentum determines.
    let spin = if a.spin == b.spin { 
        a.spin 
    } else {
        // Whichever has more energy wins
        if a.energy > b.energy { a.spin } else { b.spin }
    };
    
    Vortex {
        id: new_id,
        x,
        y,
        z,
        vx,
        vy,
        vz,
        energy: total_energy,
        frequency,
        phase,
        spin,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    
    #[test]
    fn test_attraction() {
        let a = Vortex::new(1, 0.0, 0.0, 1.0);
        let b = Vortex::new(2, 1.0, 0.0, 1.0);
        
        // At distance 1, with energy 1 each: F = 1*1/1² = 1
        assert!((a.attraction_to(&b) - 1.0).abs() < 0.001);
    }
    
    #[test]
    fn test_harmonic_detection() {
        assert!(is_harmonic(2.0));  // octave
        assert!(is_harmonic(1.5));  // fifth
        assert!(!is_harmonic(1.7)); // not harmonic
    }
}
