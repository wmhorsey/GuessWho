//! Resonance — The fundamental frequency signature
//!
//! Resonance is the sorting key of the universe. It's not about content,
//! it's about vibration. Similar resonances attract, dissimilar repel.
//!
//! Light has three fundamental properties:
//! 1. **Direction/Velocity** - Where it's going and how fast (momentum vector)
//! 2. **Wavelength/Frequency** - The color, the pitch, the λ
//! 3. **Amplitude** - The intensity, the loudness, the brightness
//!
//! All three must align for true resonance.

use std::f64::consts::PI;
use serde::{Serialize, Deserialize};

/// A 3D velocity/direction vector
#[derive(Clone, Debug, Default, Serialize, Deserialize)]
pub struct Velocity {
    pub x: f64,
    pub y: f64,
    pub z: f64,
}

impl Velocity {
    pub fn new(x: f64, y: f64, z: f64) -> Self {
        Self { x, y, z }
    }

    pub fn zero() -> Self {
        Self { x: 0.0, y: 0.0, z: 0.0 }
    }

    pub fn random() -> Self {
        use rand::Rng;
        let mut rng = rand::thread_rng();
        // Random direction on unit sphere, random speed 0.1-1.0
        let theta = rng.gen_range(0.0..2.0 * PI);
        let phi = rng.gen_range(0.0..PI);
        let speed = rng.gen_range(0.1..1.0);
        Self {
            x: speed * phi.sin() * theta.cos(),
            y: speed * phi.sin() * theta.sin(),
            z: speed * phi.cos(),
        }
    }

    pub fn magnitude(&self) -> f64 {
        (self.x * self.x + self.y * self.y + self.z * self.z).sqrt()
    }

    pub fn normalized(&self) -> Self {
        let mag = self.magnitude().max(0.001);
        Self {
            x: self.x / mag,
            y: self.y / mag,
            z: self.z / mag,
        }
    }

    /// Dot product — how aligned are two velocities?
    pub fn dot(&self, other: &Velocity) -> f64 {
        self.x * other.x + self.y * other.y + self.z * other.z
    }

    /// Alignment factor: 1.0 = same direction, -1.0 = opposite, 0.0 = perpendicular
    pub fn alignment(&self, other: &Velocity) -> f64 {
        let n1 = self.normalized();
        let n2 = other.normalized();
        n1.dot(&n2)
    }
}

/// A resonance signature — the three properties of light
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct Resonance {
    // === Vector 1: Direction/Velocity (momentum) ===
    /// The direction and speed of this resonance
    pub velocity: Velocity,

    // === Vector 2: Wavelength/Frequency ===
    /// Primary frequency (0.0 to 1.0, wrapping like a phase)
    pub wavelength: f64,
    /// Harmonic overtones that add character
    pub harmonics: [f64; 3],
    /// Phase offset — when in the cycle this resonance sits
    pub phase: f64,

    // === Vector 3: Amplitude (intensity) ===
    /// The intensity/brightness of this resonance (0.0 to 1.0)
    pub amplitude: f64,
    /// Amplitude decay rate (how quickly intensity fades)
    pub decay: f64,
}

impl Resonance {
    /// Create a new resonance with the given base wavelength
    pub fn new(wavelength: f64) -> Self {
        let wl = wavelength.rem_euclid(1.0);
        Self {
            velocity: Velocity::default(),
            wavelength: wl,
            harmonics: [
                (wl * 2.0).rem_euclid(1.0),
                (wl * 3.0).rem_euclid(1.0),
                (wl * 5.0).rem_euclid(1.0), // Prime harmonics
            ],
            phase: 0.0,
            amplitude: 1.0,
            decay: 0.001,
        }
    }

    /// Create a full resonance with all three vectors
    pub fn full(wavelength: f64, amplitude: f64, velocity: Velocity) -> Self {
        Self {
            velocity,
            amplitude: amplitude.clamp(0.0, 1.0),
            ..Self::new(wavelength)
        }
    }

    /// Create a resonance with random characteristics
    pub fn random() -> Self {
        use rand::Rng;
        let mut rng = rand::thread_rng();
        let wavelength = rng.gen::<f64>();
        let amplitude = rng.gen_range(0.3..1.0);
        let velocity = Velocity::random();
        Self::full(wavelength, amplitude, velocity)
    }

    /// Create with a specific phase
    pub fn with_phase(mut self, phase: f64) -> Self {
        self.phase = phase.rem_euclid(1.0);
        self
    }

    /// Create with a specific velocity
    pub fn with_velocity(mut self, velocity: Velocity) -> Self {
        self.velocity = velocity;
        self
    }

    /// Create with a specific amplitude
    pub fn with_amplitude(mut self, amplitude: f64) -> Self {
        self.amplitude = amplitude.clamp(0.0, 1.0);
        self
    }

    /// Calculate affinity between two resonances using all three vectors
    /// Returns a value from -1.0 (repulsion) to 1.0 (attraction)
    pub fn affinity(&self, other: &Resonance) -> f64 {
        // === Component 1: Wavelength match (most important) ===
        let wl_diff = (self.wavelength - other.wavelength).abs();
        let wl_affinity = 1.0 - (wl_diff.min(1.0 - wl_diff) * 2.0);

        // Harmonic alignment
        let harmonic_affinity: f64 = self.harmonics.iter()
            .zip(other.harmonics.iter())
            .map(|(a, b)| {
                let diff = (a - b).abs();
                1.0 - (diff.min(1.0 - diff) * 2.0)
            })
            .sum::<f64>() / 3.0;

        // Phase coherence
        let phase_diff = (self.phase - other.phase).abs();
        let phase_coherence = (phase_diff * 2.0 * PI).cos();

        // === Component 2: Velocity/Direction alignment ===
        let velocity_alignment = self.velocity.alignment(&other.velocity);
        // Speed similarity (magnitude ratio)
        let speed_ratio = {
            let s1 = self.velocity.magnitude().max(0.01);
            let s2 = other.velocity.magnitude().max(0.01);
            1.0 - ((s1 - s2).abs() / s1.max(s2)).min(1.0)
        };
        let velocity_affinity = velocity_alignment * 0.7 + speed_ratio * 0.3;

        // === Component 3: Amplitude resonance ===
        // Similar amplitudes resonate; very different ones don't
        let amp_diff = (self.amplitude - other.amplitude).abs();
        let amplitude_affinity = 1.0 - amp_diff;

        // === Combined affinity ===
        // Wavelength dominates, velocity matters, amplitude fine-tunes
        let raw_affinity = 
            wl_affinity * 0.40 +
            harmonic_affinity * 0.15 +
            phase_coherence * 0.05 +
            velocity_affinity * 0.25 +
            amplitude_affinity * 0.15;

        // Transform to attraction/repulsion range
        if raw_affinity > 0.5 {
            (raw_affinity - 0.5) * 2.0
        } else {
            (raw_affinity - 0.5) * 2.0
        }
    }

    /// Perturb the resonance slightly (drift over time)
    pub fn drift(&mut self, amount: f64) {
        self.phase = (self.phase + amount * 0.1).rem_euclid(1.0);
        // Amplitude decays slowly
        self.amplitude = (self.amplitude - self.decay * amount).max(0.01);
    }

    /// Boost amplitude (when packet is "used" or reinforced)
    pub fn boost(&mut self, amount: f64) {
        self.amplitude = (self.amplitude + amount).min(1.0);
    }

    /// Calculate a "distance" in resonance space
    pub fn distance(&self, other: &Resonance) -> f64 {
        let affinity = self.affinity(other);
        (1.0 - affinity) / 2.0
    }

    /// Get a display-friendly representation
    pub fn signature(&self) -> String {
        let speed = self.velocity.magnitude();
        format!(
            "λ{:.3}α{:.2}v{:.2}",
            self.wavelength,
            self.amplitude,
            speed
        )
    }

    /// Detailed signature for debugging
    pub fn signature_full(&self) -> String {
        format!(
            "λ{:.3} φ{:.2} α{:.2} v({:.2},{:.2},{:.2})",
            self.wavelength,
            self.phase,
            self.amplitude,
            self.velocity.x,
            self.velocity.y,
            self.velocity.z
        )
    }
}

impl Default for Resonance {
    fn default() -> Self {
        Self::new(0.5)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn identical_resonance_attracts() {
        let r1 = Resonance::new(0.5);
        let r2 = Resonance::new(0.5);
        assert!(r1.affinity(&r2) > 0.5);
    }

    #[test]
    fn opposite_resonance_repels() {
        let r1 = Resonance::new(0.0);
        let r2 = Resonance::new(0.5);
        assert!(r1.affinity(&r2) < 0.0);
    }

    #[test]
    fn velocity_alignment_matters() {
        let r1 = Resonance::full(0.5, 1.0, Velocity::new(1.0, 0.0, 0.0));
        let r2_same = Resonance::full(0.5, 1.0, Velocity::new(1.0, 0.0, 0.0));
        let r2_opp = Resonance::full(0.5, 1.0, Velocity::new(-1.0, 0.0, 0.0));
        
        assert!(r1.affinity(&r2_same) > r1.affinity(&r2_opp));
    }
}
