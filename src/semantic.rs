//! Semantic Resonance — Deriving resonance from content
//! 
//! The key insight: resonance must EMERGE from content, not be assigned randomly.
//! This creates the mapping: content → resonance → node → supernode
//! 
//! Like quantum mechanics, similar content "tunnels" to similar locations
//! because it resonates at similar frequencies.

use crate::resonance::{Resonance, Velocity};
use serde::{Deserialize, Serialize};
use std::f64::consts::PI;

/// Base resonance values for character classes
/// These are the "atomic" resonances that combine into molecular patterns
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CharacterResonance {
    pub wavelength: f64,   // Base frequency of this character class
    pub amplitude: f64,    // Energy contribution
    pub phase_shift: f64,  // How it shifts when combined
}

impl CharacterResonance {
    pub fn new(wavelength: f64, amplitude: f64, phase_shift: f64) -> Self {
        Self { wavelength, amplitude, phase_shift }
    }
}

/// The resonance table — maps characters to their base resonance
/// This is like the periodic table of semantic elements
pub struct ResonanceTable {
    // Vowels resonate in warm frequencies (longer wavelengths)
    // Consonants in cooler frequencies (shorter wavelengths)
    // Numbers have their own spectral band
    // Punctuation affects phase and velocity
}

impl ResonanceTable {
    /// Get the base resonance for a character
    /// Characters are grouped into spectral bands like light
    pub fn char_resonance(c: char) -> CharacterResonance {
        match c {
            // Vowels — the "warm" spectrum (red/infrared wavelengths)
            // Vowels carry meaning, they're the energy carriers
            'a' | 'A' => CharacterResonance::new(0.70, 1.0, 0.0),
            'e' | 'E' => CharacterResonance::new(0.75, 0.9, PI / 5.0),
            'i' | 'I' => CharacterResonance::new(0.80, 0.85, PI / 4.0),
            'o' | 'O' => CharacterResonance::new(0.72, 0.95, PI / 6.0),
            'u' | 'U' => CharacterResonance::new(0.78, 0.8, PI / 3.0),
            'y' | 'Y' => CharacterResonance::new(0.76, 0.7, PI / 7.0), // Sometimes vowel
            
            // Common consonants — the "cool" spectrum (blue/violet)
            // These shape meaning, like how blue light refracts more
            't' | 'T' => CharacterResonance::new(0.30, 0.8, 0.1),
            'n' | 'N' => CharacterResonance::new(0.32, 0.85, 0.15),
            's' | 'S' => CharacterResonance::new(0.28, 0.75, 0.2),
            'r' | 'R' => CharacterResonance::new(0.35, 0.9, 0.12),
            'l' | 'L' => CharacterResonance::new(0.38, 0.82, 0.18),
            'h' | 'H' => CharacterResonance::new(0.25, 0.5, 0.05), // Aspirated, low energy
            
            // Hard consonants — the "ultraviolet" (high frequency, sharp)
            'k' | 'K' => CharacterResonance::new(0.15, 0.7, 0.3),
            'p' | 'P' => CharacterResonance::new(0.18, 0.65, 0.25),
            'b' | 'B' => CharacterResonance::new(0.20, 0.72, 0.28),
            'd' | 'D' => CharacterResonance::new(0.22, 0.68, 0.22),
            'g' | 'G' => CharacterResonance::new(0.17, 0.66, 0.32),
            
            // Fricatives and sibilants — "X-ray" range
            'f' | 'F' => CharacterResonance::new(0.12, 0.55, 0.4),
            'v' | 'V' => CharacterResonance::new(0.14, 0.58, 0.38),
            'z' | 'Z' => CharacterResonance::new(0.10, 0.6, 0.45),
            'x' | 'X' => CharacterResonance::new(0.08, 0.5, 0.5),
            
            // Other consonants
            'c' | 'C' => CharacterResonance::new(0.24, 0.62, 0.35),
            'j' | 'J' => CharacterResonance::new(0.19, 0.48, 0.42),
            'm' | 'M' => CharacterResonance::new(0.40, 0.88, 0.08),
            'w' | 'W' => CharacterResonance::new(0.45, 0.6, 0.1),
            'q' | 'Q' => CharacterResonance::new(0.11, 0.45, 0.48),
            
            // Numbers — the "radio wave" spectrum (very long wavelengths)
            // Numbers are stable, foundational — low frequency, high amplitude
            '0' => CharacterResonance::new(0.90, 1.0, 0.0),   // Zero — maximum wavelength
            '1' => CharacterResonance::new(0.91, 0.95, 0.1),
            '2' => CharacterResonance::new(0.92, 0.90, 0.2),
            '3' => CharacterResonance::new(0.93, 0.85, 0.3),
            '4' => CharacterResonance::new(0.94, 0.80, 0.4),
            '5' => CharacterResonance::new(0.95, 0.75, 0.5),
            '6' => CharacterResonance::new(0.96, 0.70, 0.6),
            '7' => CharacterResonance::new(0.97, 0.65, 0.7),
            '8' => CharacterResonance::new(0.98, 0.60, 0.8),
            '9' => CharacterResonance::new(0.99, 0.55, 0.9),
            
            // Punctuation — phase modulators (affect rhythm, not frequency)
            '.' => CharacterResonance::new(0.50, 0.3, PI),        // Full stop — phase flip
            ',' => CharacterResonance::new(0.50, 0.2, PI / 2.0),  // Pause
            '!' => CharacterResonance::new(0.50, 0.9, PI * 1.5),  // Exclamation — energy spike
            '?' => CharacterResonance::new(0.50, 0.7, PI * 0.75), // Question — uncertainty
            ':' => CharacterResonance::new(0.50, 0.4, PI / 3.0),
            ';' => CharacterResonance::new(0.50, 0.35, PI / 2.5),
            '-' => CharacterResonance::new(0.50, 0.1, 0.0),       // Continuation
            '_' => CharacterResonance::new(0.50, 0.15, 0.0),
            
            // Brackets — containment (affect velocity/direction)
            '(' | ')' => CharacterResonance::new(0.55, 0.25, PI / 4.0),
            '[' | ']' => CharacterResonance::new(0.56, 0.28, PI / 4.0),
            '{' | '}' => CharacterResonance::new(0.57, 0.30, PI / 4.0),
            '<' | '>' => CharacterResonance::new(0.58, 0.32, PI / 4.0),
            
            // Quotes — resonance dampeners
            '"' | '\'' => CharacterResonance::new(0.52, 0.2, PI / 8.0),
            '`' => CharacterResonance::new(0.53, 0.18, PI / 8.0),
            
            // Math/logic operators — high precision frequencies
            '+' => CharacterResonance::new(0.60, 0.5, 0.0),
            '=' => CharacterResonance::new(0.61, 0.55, 0.0),
            '*' => CharacterResonance::new(0.62, 0.6, PI / 6.0),
            '/' => CharacterResonance::new(0.63, 0.45, PI / 5.0),
            '%' => CharacterResonance::new(0.64, 0.40, PI / 4.0),
            '&' => CharacterResonance::new(0.65, 0.52, PI / 7.0),
            '|' => CharacterResonance::new(0.66, 0.48, PI / 8.0),
            '^' => CharacterResonance::new(0.67, 0.42, PI / 3.0),
            
            // Whitespace — the silence between notes
            ' ' => CharacterResonance::new(0.50, 0.05, 0.0),
            '\t' => CharacterResonance::new(0.50, 0.08, 0.0),
            '\n' => CharacterResonance::new(0.50, 0.1, PI / 2.0),
            '\r' => CharacterResonance::new(0.50, 0.1, PI / 2.0),
            
            // Special symbols
            '@' => CharacterResonance::new(0.42, 0.65, PI / 5.0),
            '#' => CharacterResonance::new(0.43, 0.62, PI / 4.5),
            '$' => CharacterResonance::new(0.44, 0.70, PI / 4.0),
            
            // Unknown — mid-spectrum, low energy
            _ => CharacterResonance::new(0.50, 0.3, 0.0),
        }
    }
}

/// Compute the resonance signature for arbitrary content
/// This is the core algorithm: content → resonance
/// 
/// The metaphor: characters are like atoms, they combine into molecules (words),
/// which combine into materials (documents). The resonance emerges from the
/// interference pattern of all the component frequencies.
/// 
/// KEY INSIGHT: Order matters! "cat" ≠ "tac" because:
/// - 'c' creates wave C
/// - 'c' + 'a' = C modulated by A = CA (different from A alone)
/// - 'ca' + 't' = CA modulated by T = CAT
/// vs:
/// - 't' creates wave T  
/// - 't' + 'a' = T modulated by A = TA
/// - 'ta' + 'c' = TA modulated by C = TAC
/// 
/// Each character MODULATES the existing wave, creating interference patterns.
/// This is the quantum tunneling mechanism — sequential phase accumulation.
pub fn content_to_resonance(content: &str) -> Resonance {
    if content.is_empty() {
        return Resonance::full(0.5, 1.0, Velocity::zero());
    }
    
    let chars: Vec<char> = content.chars().collect();
    let n = chars.len() as f64;
    
    // === SEQUENTIAL WAVE MODULATION ===
    // Each character modulates the CURRENT wave state, not just adds to it
    // This is how "cat" ≠ "tac"
    
    let mut wave_state = WaveState::new();
    
    for (i, &c) in chars.iter().enumerate() {
        let cr = ResonanceTable::char_resonance(c);
        wave_state.modulate(c, &cr, i, chars.len());
    }
    
    wave_state.to_resonance(n)
}

/// Internal wave state for sequential modulation
struct WaveState {
    wavelength: f64,
    amplitude: f64,
    phase: f64,
    velocity: (f64, f64, f64),
    // Running interference pattern — this is what makes order matter
    interference_hash: f64,
}

impl WaveState {
    fn new() -> Self {
        Self {
            wavelength: 0.5,  // Start at mid-spectrum
            amplitude: 1.0,
            phase: 0.0,
            velocity: (0.0, 0.0, 0.0),
            interference_hash: 0.0,
        }
    }
    
    /// Modulate the current wave state with a new character
    /// This is where ORDER matters — each char changes the existing state
    fn modulate(&mut self, c: char, cr: &CharacterResonance, position: usize, total_len: usize) {
        let pos_weight = 1.0 / (1.0 + position as f64 * 0.1);
        let pos_ratio = position as f64 / total_len.max(1) as f64;
        let char_code = c as u32 as f64;
        
        // === WAVELENGTH: Order-dependent interference ===
        // The interference pattern depends on the CURRENT wavelength
        // So c→a→t produces different result than t→a→c
        let current_influence = self.wavelength * (1.0 + self.interference_hash * 0.0001);
        let interference = ((current_influence * cr.wavelength).sqrt() + 
                           (self.phase.sin() * 0.1)).abs();
        self.wavelength = self.wavelength * (1.0 - pos_weight * 0.4) 
                        + cr.wavelength * pos_weight * 0.3
                        + interference * pos_weight * 0.15;
        self.wavelength = self.wavelength.clamp(0.01, 0.99);
        
        // === AMPLITUDE: Multiplicative interference ===
        let harmony = 1.0 - (self.wavelength - cr.wavelength).abs();
        self.amplitude *= 0.9 + cr.amplitude * 0.1 * harmony;
        self.amplitude = self.amplitude.clamp(0.1, 2.0);
        
        // === PHASE: Strong order dependency ===
        // Phase accumulates differently based on:
        // 1. The character's intrinsic phase shift
        // 2. The CURRENT phase (nonlinear feedback)
        // 3. Position-dependent scaling
        let current_phase_influence = 1.0 + self.phase.sin() * 0.5 + self.phase.cos() * 0.3;
        let position_factor = 1.0 + pos_ratio * 0.5;
        let phase_modulation = cr.phase_shift * current_phase_influence * position_factor
                             + char_code * 0.01 * (1.0 + self.wavelength);
        self.phase += phase_modulation;
        
        // === VELOCITY: Direction emerges from sequence ===
        let angle_xy = char_code * 0.1 + self.phase * 0.5;  // Phase affects direction!
        let angle_z = cr.wavelength * PI;
        
        self.velocity.0 += angle_xy.cos() * pos_weight * cr.amplitude;
        self.velocity.1 += angle_xy.sin() * pos_weight * cr.amplitude;
        self.velocity.2 += (angle_z.cos() - 0.5) * pos_weight;
        
        // === INTERFERENCE HASH: The "memory" of the sequence ===
        // This number uniquely encodes the ORDER of characters
        // Like a rolling hash that's sensitive to position
        self.interference_hash = self.interference_hash * 31.0 
                               + char_code * (1.0 + pos_ratio)
                               + self.wavelength * 100.0;
        self.interference_hash = self.interference_hash % 10000.0;
    }
    
    fn to_resonance(&self, char_count: f64) -> Resonance {
        // Normalize phase
        let phase = self.phase % (2.0 * PI);
        
        // Normalize velocity
        let vel_mag = (self.velocity.0.powi(2) + self.velocity.1.powi(2) + self.velocity.2.powi(2)).sqrt();
        let velocity = if vel_mag > 0.001 {
            Velocity {
                x: self.velocity.0 / vel_mag * 0.5,
                y: self.velocity.1 / vel_mag * 0.5,
                z: self.velocity.2 / vel_mag * 0.5,
            }
        } else {
            Velocity::zero()
        };
        
        // Use interference hash to create unique harmonics
        let h1 = ((self.interference_hash * 0.0001) % 1.0) + 0.5;
        let h2 = ((self.interference_hash * 0.00001) % 1.0) + 0.25;
        let h3 = ((self.interference_hash * 0.000001) % 1.0) + 0.125;
        
        let mut resonance = Resonance::full(self.wavelength, self.amplitude, velocity);
        resonance.phase = phase;
        resonance.harmonics = [h1, h2, h3];
        resonance.decay = 1.0 / (1.0 + char_count * 0.01);
        
        resonance
    }
}

/// Word-level resonance — treats words as fundamental units
/// More useful for semantic similarity than character-level
pub fn words_to_resonance(content: &str) -> Resonance {
    let words: Vec<&str> = content.split_whitespace().collect();
    
    if words.is_empty() {
        return content_to_resonance(content);
    }
    
    // Each WORD has its own resonance, then words combine
    let mut wave_state = WaveState::new();
    
    for (i, word) in words.iter().enumerate() {
        // Get the word's intrinsic resonance
        let word_resonance = content_to_resonance(word);
        
        // Convert to a CharacterResonance-like structure for modulation
        let word_cr = CharacterResonance {
            wavelength: word_resonance.wavelength,
            amplitude: word_resonance.amplitude,
            phase_shift: word_resonance.phase,
        };
        
        // Use first char of word as the "character" for velocity calculation
        let first_char = word.chars().next().unwrap_or(' ');
        wave_state.modulate(first_char, &word_cr, i, words.len());
    }
    
    wave_state.to_resonance(words.len() as f64)
}

/// Hybrid resonance — combines character and word level
/// This captures both local patterns (spelling) and global patterns (meaning)
pub fn hybrid_resonance(content: &str) -> Resonance {
    let char_res = content_to_resonance(content);
    let word_res = words_to_resonance(content);
    
    // Blend: word-level dominates for semantics, char-level for uniqueness
    Resonance::full(
        word_res.wavelength * 0.7 + char_res.wavelength * 0.3,
        word_res.amplitude * 0.6 + char_res.amplitude * 0.4,
        Velocity {
            x: word_res.velocity.x * 0.7 + char_res.velocity.x * 0.3,
            y: word_res.velocity.y * 0.7 + char_res.velocity.y * 0.3,
            z: word_res.velocity.z * 0.7 + char_res.velocity.z * 0.3,
        },
    )
}

/// Compute resonance similarity — how likely two contents will cluster
pub fn resonance_similarity(content_a: &str, content_b: &str) -> f64 {
    let res_a = content_to_resonance(content_a);
    let res_b = content_to_resonance(content_b);
    res_a.affinity(&res_b)
}

/// Token-level resonance for larger content
/// Breaks content into words and computes hierarchical resonance
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TokenResonance {
    pub token: String,
    pub resonance: Resonance,
}

impl TokenResonance {
    pub fn from_token(token: &str) -> Self {
        Self {
            token: token.to_string(),
            resonance: content_to_resonance(token),
        }
    }
}

/// Document-level resonance — hierarchical composition
/// Words → Sentences → Paragraphs → Document
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DocumentResonance {
    pub tokens: Vec<TokenResonance>,
    pub aggregate: Resonance,
}

impl DocumentResonance {
    pub fn from_content(content: &str) -> Self {
        // Tokenize on whitespace and punctuation
        let tokens: Vec<TokenResonance> = content
            .split_whitespace()
            .map(|word| TokenResonance::from_token(word))
            .collect();
        
        // Aggregate resonance from tokens
        let aggregate = if tokens.is_empty() {
            content_to_resonance("")
        } else {
            // Combine token resonances with interference
            let mut wavelength = 0.0;
            let mut amplitude = 0.0;
            let mut vx = 0.0;
            let mut vy = 0.0;
            let mut vz = 0.0;
            
            for (i, tr) in tokens.iter().enumerate() {
                let weight = 1.0 / (1.0 + i as f64 * 0.05);
                wavelength += tr.resonance.wavelength * weight;
                amplitude += tr.resonance.amplitude * weight;
                vx += tr.resonance.velocity.x * weight;
                vy += tr.resonance.velocity.y * weight;
                vz += tr.resonance.velocity.z * weight;
            }
            
            let total_weight: f64 = (0..tokens.len())
                .map(|i| 1.0 / (1.0 + i as f64 * 0.05))
                .sum();
            
            Resonance::full(
                wavelength / total_weight,
                amplitude / total_weight,
                Velocity {
                    x: vx / total_weight,
                    y: vy / total_weight,
                    z: vz / total_weight,
                },
            )
        };
        
        Self { tokens, aggregate }
    }
    
    /// How well does this document resonate with a node's frequency?
    pub fn affinity_to(&self, node_resonance: &Resonance) -> f64 {
        self.aggregate.affinity(node_resonance)
    }
}

/// Content types derive their packet types from resonance patterns
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum ContentClass {
    Numeric,      // Dominated by numbers
    Textual,      // Dominated by letters
    Code,         // Mix of symbols and alphanumerics
    Data,         // Structured (high punctuation)
    Query,        // Questions (contains ?)
    Command,      // Imperatives (contains !)
}

impl ContentClass {
    pub fn classify(content: &str) -> Self {
        let mut letters = 0;
        let mut digits = 0;
        let mut symbols = 0;
        let mut questions = 0;
        let mut exclamations = 0;
        
        for c in content.chars() {
            if c.is_alphabetic() { letters += 1; }
            else if c.is_numeric() { digits += 1; }
            else if c == '?' { questions += 1; }
            else if c == '!' { exclamations += 1; }
            else if !c.is_whitespace() { symbols += 1; }
        }
        
        let total = (letters + digits + symbols).max(1);
        
        if questions > 0 {
            ContentClass::Query
        } else if exclamations > 0 {
            ContentClass::Command
        } else if digits as f64 / total as f64 > 0.5 {
            ContentClass::Numeric
        } else if symbols as f64 / total as f64 > 0.3 {
            if letters > digits { ContentClass::Code } else { ContentClass::Data }
        } else {
            ContentClass::Textual
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    
    #[test]
    fn test_order_matters_cat_vs_tac() {
        // THE KEY TEST: Order must matter!
        // "cat" and "tac" have the same characters but different resonance
        let cat = content_to_resonance("cat");
        let tac = content_to_resonance("tac");
        
        // They should have DIFFERENT wavelengths
        assert!(
            (cat.wavelength - tac.wavelength).abs() > 0.001,
            "cat (λ={:.4}) and tac (λ={:.4}) must have different wavelengths!",
            cat.wavelength, tac.wavelength
        );
        
        // They should have DIFFERENT phases
        assert!(
            (cat.phase - tac.phase).abs() > 0.01,
            "cat (φ={:.4}) and tac (φ={:.4}) must have different phases!",
            cat.phase, tac.phase
        );
        
        // Different velocity directions
        let vel_diff = ((cat.velocity.x - tac.velocity.x).powi(2) 
                      + (cat.velocity.y - tac.velocity.y).powi(2)).sqrt();
        assert!(
            vel_diff > 0.01,
            "cat and tac must have different velocity vectors!"
        );
        
        println!("cat: λ={:.4}, α={:.4}, φ={:.4}, v=({:.3},{:.3},{:.3})",
            cat.wavelength, cat.amplitude, cat.phase,
            cat.velocity.x, cat.velocity.y, cat.velocity.z);
        println!("tac: λ={:.4}, α={:.4}, φ={:.4}, v=({:.3},{:.3},{:.3})",
            tac.wavelength, tac.amplitude, tac.phase,
            tac.velocity.x, tac.velocity.y, tac.velocity.z);
    }
    
    #[test]
    fn test_similar_content_resonates() {
        let sim1 = resonance_similarity("hello", "hallo");
        let sim2 = resonance_similarity("hello", "world");
        let sim3 = resonance_similarity("hello", "12345");
        
        // Similar words should resonate more than different ones
        assert!(sim1 > sim2, "hello/hallo should resonate more than hello/world");
        assert!(sim2 > sim3, "hello/world should resonate more than hello/12345");
    }
    
    #[test]
    fn test_word_order_matters() {
        // "the dog bit the man" vs "the man bit the dog"
        let a = words_to_resonance("the dog bit the man");
        let b = words_to_resonance("the man bit the dog");
        
        // Same words, different order = different resonance
        assert!(
            (a.wavelength - b.wavelength).abs() > 0.001 ||
            (a.phase - b.phase).abs() > 0.01,
            "Word order must affect resonance!"
        );
    }
    
    #[test]
    fn test_content_classification() {
        assert_eq!(ContentClass::classify("12345"), ContentClass::Numeric);
        assert_eq!(ContentClass::classify("hello world"), ContentClass::Textual);
        assert_eq!(ContentClass::classify("fn main() {}"), ContentClass::Code);
        assert_eq!(ContentClass::classify("what is this?"), ContentClass::Query);
    }
}
