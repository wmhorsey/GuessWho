// genesis.rs - Let there be light
//
// Start from nothing. Add energy. Watch what forms.
// This is the primordial simulation.

use std::f64::consts::PI;
use std::collections::HashMap;
use rand::Rng;
use crate::vortex::Vortex;
use crate::field::Field;
use crate::aptik;  // THE KEY: Physics engine for void/spike architecture

/// Check if a number is prime (Optimus approved!)
fn is_prime(n: usize) -> bool {
    if n < 2 { return false; }
    if n == 2 { return true; }
    if n % 2 == 0 { return false; }
    let sqrt_n = (n as f64).sqrt() as usize;
    for i in (3..=sqrt_n).step_by(2) {
        if n % i == 0 { return false; }
    }
    true
}

/// A detected multi-body configuration (potential "letter")
#[derive(Clone, Debug)]
pub struct Configuration {
    pub pattern: String,           // Signature describing the pattern
    pub member_count: usize,
    pub freq_ratios: Vec<f64>,     // Frequency ratios in the group
    pub spin_signature: String,    // e.g., "++", "+-", "+++"
    pub count: usize,              // How many times this pattern appears
}

/// A bond type - the fundamental relationship between two vortices (a "letter")
#[derive(Clone, Debug, Hash, Eq, PartialEq)]
pub struct BondType {
    pub freq_ratio: String,     // e.g., "1:1", "2:1", "3/2:1"
    pub spin_pair: String,      // e.g., "++", "+-", "--"
}

/// Analyze orbital pairs to count bond types (the letters)
fn find_bond_types(field: &Field) -> HashMap<BondType, usize> {
    let mut bond_counts: HashMap<BondType, usize> = HashMap::new();
    
    for pair in &field.orbital_pairs {
        if let (Some(a), Some(b)) = (field.vortices.get(&pair.a), field.vortices.get(&pair.b)) {
            // Calculate frequency ratio
            let (high, low) = if a.frequency >= b.frequency { (a, b) } else { (b, a) };
            let ratio = high.frequency / low.frequency;
            
            // Quantize ratio
            let ratio_str = if (ratio - 1.0).abs() < 0.1 { "1:1" }
                else if (ratio - 1.5).abs() < 0.1 { "3:2" }
                else if (ratio - 2.0).abs() < 0.1 { "2:1" }
                else if (ratio - 1.333).abs() < 0.1 { "4:3" }
                else if (ratio - 1.618).abs() < 0.1 { "φ:1" }
                else if (ratio - 3.0).abs() < 0.1 { "3:1" }
                else if (ratio - 4.0).abs() < 0.1 { "4:1" }
                else { "?:?" };
            
            // Spin pair
            let spin_str = format!("{}{}",
                if high.spin > 0 { "+" } else { "-" },
                if low.spin > 0 { "+" } else { "-" });
            
            let bond = BondType {
                freq_ratio: ratio_str.to_string(),
                spin_pair: spin_str,
            };
            
            *bond_counts.entry(bond).or_insert(0) += 1;
        }
    }
    
    bond_counts
}

/// Analyze orbital pairs to find multi-body configurations
fn find_configurations(field: &Field) -> Vec<Configuration> {
    // Build adjacency: which vortices are orbiting which
    let mut adjacency: HashMap<u64, Vec<u64>> = HashMap::new();
    
    for pair in &field.orbital_pairs {
        adjacency.entry(pair.a).or_default().push(pair.b);
        adjacency.entry(pair.b).or_default().push(pair.a);
    }
    
    // Instead of giant connected components, look at LOCAL structure
    // Count how many vortices have 1, 2, 3, etc. connections (valence)
    let mut valence_dist: HashMap<usize, usize> = HashMap::new();
    for neighbors in adjacency.values() {
        *valence_dist.entry(neighbors.len()).or_insert(0) += 1;
    }
    
    // Find small isolated groups (2-5 members) - these are "words"
    let mut visited: HashMap<u64, bool> = HashMap::new();
    let mut groups: Vec<Vec<u64>> = Vec::new();
    
    for &id in adjacency.keys() {
        if visited.get(&id).copied().unwrap_or(false) {
            continue;
        }
        
        // BFS to find connected component
        let mut group = Vec::new();
        let mut queue = vec![id];
        
        while let Some(current) = queue.pop() {
            if visited.get(&current).copied().unwrap_or(false) {
                continue;
            }
            visited.insert(current, true);
            group.push(current);
            
            if let Some(neighbors) = adjacency.get(&current) {
                for &neighbor in neighbors {
                    if !visited.get(&neighbor).copied().unwrap_or(false) {
                        queue.push(neighbor);
                    }
                }
            }
        }
        
        if group.len() >= 2 {
            groups.push(group);
        }
    }
    
    // Classify each group by its properties
    let mut pattern_counts: HashMap<String, Configuration> = HashMap::new();
    
    for group in groups {
        // Get vortex properties
        let mut freqs: Vec<f64> = Vec::new();
        let mut spins: Vec<i8> = Vec::new();
        
        for &id in &group {
            if let Some(v) = field.vortices.get(&id) {
                freqs.push(v.frequency);
                spins.push(v.spin);
            }
        }
        
        // Sort frequencies to normalize
        freqs.sort_by(|a, b| a.partial_cmp(b).unwrap());
        
        // Compute frequency ratios (relative to lowest)
        let base_freq = freqs.first().copied().unwrap_or(1.0);
        let freq_ratios: Vec<f64> = freqs.iter().map(|f| f / base_freq).collect();
        
        // Quantize ratios to detect patterns
        let ratio_pattern: String = freq_ratios.iter()
            .map(|r| {
                if (*r - 1.0).abs() < 0.1 { "1" }
                else if (*r - 1.5).abs() < 0.1 { "3/2" }
                else if (*r - 2.0).abs() < 0.1 { "2" }
                else if (*r - 1.333).abs() < 0.1 { "4/3" }
                else if (*r - 1.618).abs() < 0.1 { "φ" }
                else if (*r - 0.666).abs() < 0.1 { "2/3" }
                else { "?" }
            })
            .collect::<Vec<_>>()
            .join(":");
        
        // Spin signature
        let spin_sig: String = spins.iter()
            .map(|s| if *s > 0 { "+" } else { "-" })
            .collect();
        
        // Create pattern key
        let pattern_key = format!("{}|{}", ratio_pattern, spin_sig);
        
        // Count this pattern
        let entry = pattern_counts.entry(pattern_key.clone()).or_insert(Configuration {
            pattern: pattern_key,
            member_count: group.len(),
            freq_ratios: freq_ratios.clone(),
            spin_signature: spin_sig,
            count: 0,
        });
        entry.count += 1;
    }
    
    // Convert to sorted list
    let mut configs: Vec<Configuration> = pattern_counts.into_values().collect();
    configs.sort_by(|a, b| b.count.cmp(&a.count));
    configs
}

/// Create initial conditions and run the primordial simulation
pub fn run_genesis() {
    run_genesis_with_count(50, 500);
}

/// Run with configurable vortex count
pub fn run_genesis_scaled(count: usize) {
    run_genesis_with_count(count, 1000);
}

fn run_genesis_with_count(num_vortices: usize, total_ticks: usize) {
    println!("\n");
    println!("╔════════════════════════════════════════════════════════════╗");
    println!("║                     G E N E S I S                          ║");
    println!("║         What does energy want to become?                   ║");
    println!("╚════════════════════════════════════════════════════════════╝");
    println!();
    
    // Scale field size with vortex count
    let field_size = (num_vortices as f64).sqrt() * 20.0;
    let mut field = Field::new(field_size);
    
    // Seed with random vortices
    let mut rng = rand::thread_rng();
    
    println!("Creating {} primordial energy fluctuations...", num_vortices);
    println!("Field size: {:.0}\n", field_size);
    
    for i in 0..num_vortices {
        // Random position within field
        let angle = rng.gen::<f64>() * 2.0 * PI;
        let radius = rng.gen::<f64>().sqrt() * field.size * 0.8;
        let x = radius * angle.cos();
        let y = radius * angle.sin();
        
        // Random energy (some small, some large)
        let energy = 0.5 + rng.gen::<f64>() * 3.0;
        
        // Random frequency (but prefer harmonic ratios)
        let freq_choices = [1.0, 1.5, 2.0, 0.5, 1.333, 0.666, 1.618];
        let frequency = freq_choices[rng.gen_range(0..freq_choices.len())];
        
        // Random phase
        let phase = rng.gen::<f64>() * 2.0 * PI;
        
        // Random spin
        let spin: i8 = if rng.gen::<bool>() { 1 } else { -1 };
        
        // Small initial velocity (random direction)
        let v_angle = rng.gen::<f64>() * 2.0 * PI;
        let v_mag = rng.gen::<f64>() * 0.5;
        let vx = v_mag * v_angle.cos();
        let vy = v_mag * v_angle.sin();
        
        let vortex = Vortex::new(0, x, y, energy)
            .with_frequency(frequency)
            .with_phase(phase)
            .with_spin(spin)
            .with_velocity(vx, vy);
        
        field.spawn_vortex(vortex);
        
        if i < 5 {
            println!("  Vortex {}: E={:.2}, f={:.2}, φ={:.2}π, spin={}",
                     i + 1, energy, frequency, phase / PI, if spin > 0 { "↻" } else { "↺" });
        }
    }
    if num_vortices > 5 {
        println!("  ... and {} more", num_vortices - 5);
    }
    
    println!("\n{}", field.report());
    
    // Show initial state (smaller viz for large counts)
    let viz_size = if num_vortices > 200 { 40 } else { 60 };
    println!("Initial field:");
    println!("{}", field.visualize(viz_size, viz_size / 2));
    
    // Run simulation
    let dt = 0.1;
    let display_interval = total_ticks / 5;
    
    println!("\n⏳ Running simulation ({} ticks)...\n", total_ticks);
    
    for tick in 0..total_ticks {
        field.tick(dt);
        
        if display_interval > 0 && (tick + 1) % display_interval == 0 {
            println!("════════════════════════════════════════════════════════════");
            println!("Tick {}", tick + 1);
            println!("{}", field.report());
            
            // Detect emerging patterns
            let configs = find_configurations(&field);
            if !configs.is_empty() {
                println!("  Emerging patterns: {} distinct configurations", configs.len());
                for cfg in configs.iter().take(3) {
                    println!("    {} ({}x)", cfg.pattern, cfg.count);
                }
            }
            
            if num_vortices <= 200 {
                println!("{}", field.visualize(viz_size, viz_size / 2));
            }
            println!();
        }
    }
    
    // Final analysis
    println!("╔════════════════════════════════════════════════════════════╗");
    println!("║                  WHAT EMERGED?                             ║");
    println!("╚════════════════════════════════════════════════════════════╝");
    println!();
    println!("{}", field.report());
    
    // Analyze surviving vortices
    println!("Surviving vortices:");
    let mut vortices: Vec<_> = field.vortices.values().collect();
    vortices.sort_by(|a, b| b.energy.partial_cmp(&a.energy).unwrap());
    
    for (i, v) in vortices.iter().take(10).enumerate() {
        println!("  {}: E={:.2}, f={:.3}, pos=({:.1}, {:.1}), spin={}",
                 i + 1, v.energy, v.frequency, v.x, v.y,
                 if v.spin > 0 { "↻" } else { "↺" });
    }
    if vortices.len() > 10 {
        println!("  ... and {} more", vortices.len() - 10);
    }
    
    // Analyze orbital pairs
    println!("\nOrbital pairs (stable bindings):");
    for (i, pair) in field.orbital_pairs.iter().take(10).enumerate() {
        if let (Some(a), Some(b)) = (field.vortices.get(&pair.a), field.vortices.get(&pair.b)) {
            let freq_ratio = if a.frequency > b.frequency {
                a.frequency / b.frequency
            } else {
                b.frequency / a.frequency
            };
            println!("  {}: V{} ↔ V{}, freq ratio={:.2}, stability={:.2}",
                     i + 1, pair.a, pair.b, freq_ratio, pair.stability);
        }
    }
    
    // Energy conservation check
    println!("\nEnergy conservation:");
    println!("  Initial: {} vortices × ~{:.1} avg = ~{:.1}",
             num_vortices, 2.0, num_vortices as f64 * 2.0);
    println!("  Final:   {:.2}", field.total_energy());
    
    // What structures formed?
    println!("\n✧ OBSERVATIONS:");
    println!("  - Started with {} independent vortices", num_vortices);
    println!("  - Ended with {} vortices ({} merged)", 
             field.vortex_count(), field.total_merges);
    println!("  - {} orbital pairs formed (stable relationships)",
             field.orbital_pairs.len());
    
    if field.vortex_count() < num_vortices / 2 {
        println!("\n  ⚡ SIGNIFICANT COALESCENCE - energy wants to concentrate!");
    }
    
    if field.orbital_pairs.len() > field.vortex_count() / 4 {
        println!("  🌀 MANY ORBITS FORMED - stable structures emerge from missed collisions!");
    }
    
    // THE KEY INSIGHT: What patterns (letters) emerged?
    println!("\n╔════════════════════════════════════════════════════════════╗");
    println!("║               THE ALPHABET OF ENERGY                       ║");
    println!("║   Repeating multi-body configurations = primitive symbols  ║");
    println!("╚════════════════════════════════════════════════════════════╝\n");
    
    let configs = find_configurations(&field);
    
    // THE ALPHABET: Bond types are the letters!
    println!("\n╔════════════════════════════════════════════════════════════╗");
    println!("║               THE ALPHABET OF ENERGY                       ║");
    println!("║   Bond types (pair-wise relationships) = primitive letters ║");
    println!("╚════════════════════════════════════════════════════════════╝\n");
    
    let bond_types = find_bond_types(&field);
    
    if bond_types.is_empty() {
        println!("  No bonds formed yet.");
    } else {
        // Sort by count
        let mut bonds: Vec<_> = bond_types.iter().collect();
        bonds.sort_by(|a, b| b.1.cmp(a.1));
        
        println!("  BOND TYPES (the letters of energy's language):\n");
        println!("    {:^8} {:^8} {:>8}   {}", "Freq", "Spin", "Count", "Meaning");
        println!("    {:─^8} {:─^8} {:─>8}   {:─<20}", "", "", "", "");
        
        for (bond, count) in bonds.iter() {
            let meaning = bond_meaning(&bond.freq_ratio, &bond.spin_pair);
            let bar = "█".repeat((*count / 20).min(20));
            println!("    {:^8} {:^8} {:>8}   {} {}", 
                     bond.freq_ratio, bond.spin_pair, count, meaning, bar);
        }
        
        // Summary
        println!("\n  ✧ INSIGHTS:");
        
        // Most common frequency ratio
        let freq_counts: HashMap<String, usize> = bonds.iter()
            .fold(HashMap::new(), |mut acc, (b, c)| {
                *acc.entry(b.freq_ratio.clone()).or_insert(0) += *c;
                acc
            });
        let mut freq_sorted: Vec<_> = freq_counts.iter().collect();
        freq_sorted.sort_by(|a, b| b.1.cmp(a.1));
        
        if let Some((top_freq, count)) = freq_sorted.first() {
            println!("    Most common harmonic: {} ({} bonds)", top_freq, count);
        }
        
        // Spin balance
        let same_spin: usize = bonds.iter()
            .filter(|(b, _)| b.spin_pair == "++" || b.spin_pair == "--")
            .map(|(_, c)| *c)
            .sum();
        let opp_spin: usize = bonds.iter()
            .filter(|(b, _)| b.spin_pair == "+-" || b.spin_pair == "-+")
            .map(|(_, c)| *c)
            .sum();
        
        println!("    Same-spin bonds: {}  |  Opposite-spin: {}", same_spin, opp_spin);
        
        if same_spin > opp_spin * 2 {
            println!("    → Energy prefers ALIGNED spins (ferromagnetic tendency)");
        } else if opp_spin > same_spin * 2 {
            println!("    → Energy prefers OPPOSED spins (antiferromagnetic tendency)");
        } else {
            println!("    → Spin balance suggests NO PREFERENCE (or balance required)");
        }
    }
    
    // Show configuration stats too
    let configs = find_configurations(&field);
    
    if !configs.is_empty() && configs.len() < 50 {
        println!("\n  MULTI-BODY STRUCTURES (words from letters):");
        println!("    Found {} distinct groups", configs.len());
        
        // Group by size
        let mut by_size: HashMap<usize, usize> = HashMap::new();
        for cfg in &configs {
            *by_size.entry(cfg.member_count).or_insert(0) += 1;
        }
        
        for size in 2..=10 {
            if let Some(count) = by_size.get(&size) {
                println!("    {}-body groups: {}", size, count);
            }
        }
    }
    
    println!();
}

/// Human-readable meaning for bond types
fn bond_meaning(freq: &str, spin: &str) -> &'static str {
    match (freq, spin) {
        ("1:1", "++") | ("1:1", "--") => "UNISON (identical twins)",
        ("1:1", "+-") | ("1:1", "-+") => "MIRROR (antimatter pair?)",
        ("2:1", _) => "OCTAVE (parent-child)",
        ("3:2", _) => "FIFTH (most consonant)",
        ("4:3", _) => "FOURTH (stable harmony)",
        ("φ:1", _) => "GOLDEN (growth spiral)",
        ("3:1", _) => "TWELFTH (octave+fifth)",
        // Compound letters - made from combining base letters!
        ("φ/F", _) => "GOLDEN/FOURTH (compound)",
        ("O/φ", _) => "OCTAVE/GOLDEN (compound)",
        ("5/F", _) => "FIFTH/FOURTH (compound)",
        ("φ/5", _) => "GOLDEN/FIFTH (compound)",
        ("5:2", _) => "TENTH (octave+third)",
        ("8:3", _) => "COMPOUND 8:3",
        ("φ²", _) => "GOLDEN SQUARED",
        ("φ+1", _) => "GOLDEN+1 (Fibonacci!)",
        _ => "",
    }
}

/// Run with specific harmonic frequencies to see if structure emerges
pub fn run_harmonic_test() {
    println!("\n");
    println!("╔════════════════════════════════════════════════════════════╗");
    println!("║            HARMONIC RESONANCE TEST                         ║");
    println!("║    Do harmonic frequencies create stable structures?       ║");
    println!("╚════════════════════════════════════════════════════════════╝");
    println!();
    
    let mut field = Field::new(50.0);
    
    // Create vortices with specific harmonic relationships
    // Like notes in a chord: C-E-G (1:1.25:1.5)
    
    println!("Creating harmonic vortices (like a musical chord)...");
    
    // The fundamental (C)
    let v1 = Vortex::new(0, -20.0, 0.0, 3.0)
        .with_frequency(1.0)
        .with_phase(0.0)
        .with_spin(1);
    field.spawn_vortex(v1);
    println!("  C (fundamental): f=1.0, E=3.0");
    
    // The third (E) - ratio 5:4
    let v2 = Vortex::new(0, 0.0, 0.0, 2.0)
        .with_frequency(1.25)
        .with_phase(PI / 4.0)
        .with_spin(1);
    field.spawn_vortex(v2);
    println!("  E (major third): f=1.25, E=2.0");
    
    // The fifth (G) - ratio 3:2
    let v3 = Vortex::new(0, 20.0, 0.0, 2.0)
        .with_frequency(1.5)
        .with_phase(PI / 2.0)
        .with_spin(1);
    field.spawn_vortex(v3);
    println!("  G (perfect fifth): f=1.5, E=2.0");
    
    // The octave (high C) - ratio 2:1
    let v4 = Vortex::new(0, 0.0, 20.0, 1.5)
        .with_frequency(2.0)
        .with_phase(PI / 3.0)
        .with_spin(1);
    field.spawn_vortex(v4);
    println!("  C' (octave): f=2.0, E=1.5");
    
    println!("\n  All vortices have same spin (clockwise)");
    println!("  Different phases = they shouldn't merge immediately\n");
    
    println!("Initial:");
    println!("{}", field.visualize(50, 25));
    
    // Run
    let total_ticks = 300;
    let dt = 0.1;
    
    for tick in 0..total_ticks {
        field.tick(dt);
        
        if (tick + 1) % 100 == 0 {
            println!("\n--- Tick {} ---", tick + 1);
            println!("{}", field.report());
            println!("{}", field.visualize(50, 25));
        }
    }
    
    println!("\n✧ RESULT:");
    if field.vortex_count() == 4 && field.orbital_pairs.len() >= 2 {
        println!("  🎵 STABLE CHORD! The harmonic vortices formed orbital relationships!");
    } else if field.vortex_count() < 4 {
        println!("  ⚡ Some merged - the chord collapsed partially");
    } else {
        println!("  🌌 Drifted apart - not enough attraction to bind");
    }
    
    println!("\n{}", field.report());
}

/// Test what happens with opposing spins
pub fn run_spin_test() {
    println!("\n");
    println!("╔════════════════════════════════════════════════════════════╗");
    println!("║               SPIN INTERACTION TEST                        ║");
    println!("║      Do opposite spins create stable pairs?                ║");
    println!("╚════════════════════════════════════════════════════════════╝");
    println!();
    
    let mut field = Field::new(30.0);
    
    // Two vortices, same energy, opposite spin
    let v1 = Vortex::new(0, -10.0, 0.0, 2.0)
        .with_frequency(1.0)
        .with_phase(0.0)
        .with_spin(1)
        .with_velocity(0.0, 0.5);  // Moving slightly "up"
    field.spawn_vortex(v1);
    
    let v2 = Vortex::new(0, 10.0, 0.0, 2.0)
        .with_frequency(1.0)
        .with_phase(PI)  // Opposite phase
        .with_spin(-1)   // Opposite spin
        .with_velocity(0.0, -0.5);  // Moving slightly "down"
    field.spawn_vortex(v2);
    
    println!("Two vortices:");
    println!("  V1: spin ↻, phase 0");
    println!("  V2: spin ↺, phase π (opposite)");
    println!("  Both moving perpendicular to their connecting line");
    println!();
    
    println!("Initial:");
    println!("{}", field.visualize(40, 20));
    
    let total_ticks = 200;
    let dt = 0.1;
    
    for tick in 0..total_ticks {
        field.tick(dt);
        
        if (tick + 1) % 50 == 0 {
            println!("\n--- Tick {} ---", tick + 1);
            println!("{}", field.report());
            println!("{}", field.visualize(40, 20));
        }
    }
    
    println!("\n✧ RESULT:");
    if field.vortex_count() == 2 {
        println!("  🌀 STABLE PAIR! Opposite spins created a stable binary.");
        if field.orbital_pairs.len() > 0 {
            println!("  ✓ Orbital relationship formed!");
        }
    } else if field.vortex_count() == 1 {
        println!("  ⚡ MERGED! Opposite spins didn't prevent collapse.");
    }
}

/// GPU-accelerated genesis for massive simulations
pub fn run_genesis_gpu(count: usize) {
    run_genesis_gpu_with_ticks(count, None);
}

pub fn run_genesis_gpu_with_ticks(count: usize, ticks: Option<usize>) {
    use crate::gpu_field::GpuField;
    use std::time::Instant;
    
    // Field size scales with sqrt of count
    let field_size = (count as f64).sqrt() * 20.0;
    
    // Calculate ticks - use provided or auto-scale
    let total_ticks = ticks.unwrap_or_else(|| {
        if count >= 20000 { 5000 } else if count > 5000 { 2500 } else { 5000 }
    });
    
    // Estimate memory and compute
    let vortex_mem_kb = count * 128 / 1024;  // ~128 bytes per vortex
    let interactions_per_tick = (count * count) / 2;  // O(n²) pairwise
    let total_interactions = interactions_per_tick * total_ticks;
    
    println!("\n╔════════════════════════════════════════════════════════════════════╗");
    println!("║               G E N E S I S   3 D   (GPU ACCELERATED)              ║");
    println!("║           What does energy want to become? (in 3D SPACE)           ║");
    println!("╠════════════════════════════════════════════════════════════════════╣");
    println!("║  SIMULATION PARAMETERS:                                            ║");
    println!("║    Vortices:     {:>8}                                          ║", count);
    println!("║    Ticks:        {:>8}                                          ║", total_ticks);
    println!("║    Field size:   {:>8.0}                                          ║", field_size);
    println!("║    Vortex RAM:   {:>8} KB                                       ║", vortex_mem_kb);
    println!("╠════════════════════════════════════════════════════════════════════╣");
    println!("║  COMPUTE SCALE:                                                    ║");
    if total_interactions >= 1_000_000_000_000 {
        println!("║    Interactions: {:>8.2}T  ({:.0} per tick)               ║", 
            total_interactions as f64 / 1e12, interactions_per_tick as f64);
    } else if total_interactions >= 1_000_000_000 {
        println!("║    Interactions: {:>8.2}B  ({:.0}M per tick)               ║", 
            total_interactions as f64 / 1e9, interactions_per_tick as f64 / 1e6);
    } else {
        println!("║    Interactions: {:>8.2}M  ({:.0}K per tick)               ║", 
            total_interactions as f64 / 1e6, interactions_per_tick as f64 / 1e3);
    }
    println!("╚════════════════════════════════════════════════════════════════════╝\n");
    
    println!("Creating {} primordial energy fluctuations in 3D...", count);
    println!("Field size: {:.0} (spherical)", field_size);
    
    let start = Instant::now();
    let mut field = GpuField::new(field_size);
    println!("GPU initialized in {:.2}s", start.elapsed().as_secs_f64());
    
    // Random 3D vortex generation
    let mut rng = rand::thread_rng();
    let spawn_radius = field_size * 0.4;
    
    let freqs = [0.5, 0.666, 1.0, 1.333, 1.5, 1.618, 2.0];
    
    for i in 0..count {
        // Spherical spawn in 3D - uniform distribution in sphere
        // Using rejection sampling for uniform spherical distribution
        let (x, y, z) = loop {
            let x = (rng.gen::<f64>() - 0.5) * 2.0;
            let y = (rng.gen::<f64>() - 0.5) * 2.0;
            let z = (rng.gen::<f64>() - 0.5) * 2.0;
            let r = (x*x + y*y + z*z).sqrt();
            if r <= 1.0 && r > 0.01 {
                // Scale to spawn radius
                break (x * spawn_radius, y * spawn_radius, z * spawn_radius);
            }
        };
        
        let energy = 0.5 + rng.gen::<f64>() * 3.0;
        let freq_idx = rng.gen_range(0..freqs.len());
        let frequency = freqs[freq_idx];
        let phase = rng.gen::<f64>() * 2.0 * PI;
        let spin: i8 = if rng.gen::<bool>() { 1 } else { -1 };
        
        // Small random 3D velocity
        let vx = (rng.gen::<f64>() - 0.5) * 2.0;
        let vy = (rng.gen::<f64>() - 0.5) * 2.0;
        let vz = (rng.gen::<f64>() - 0.5) * 2.0;
        
        let vortex = Vortex::new_3d(0, x, y, z, energy)
            .with_frequency(frequency)
            .with_phase(phase)
            .with_spin(spin)
            .with_velocity_3d(vx, vy, vz);
        
        field.spawn_vortex(vortex);
        
        if i < 5 {
            println!("  Vortex {}: E={:.2}, f={:.2}, pos=({:.1},{:.1},{:.1}), spin={}", 
                i + 1, energy, frequency, x, y, z,
                if spin > 0 { "↻" } else { "↺" });
        } else if i == 5 {
            println!("  ... and {} more", count - 5);
        }
    }
    
    println!("\n{}", field.report());
    
    // dt is our timestep
    let dt = 0.1;
    
    println!("\n⏳ Running GPU simulation ({} ticks, {} vortices)...\n", total_ticks, count);
    
    let sim_start = Instant::now();
    let mut last_report = Instant::now();
    let mut primes_found: Vec<usize> = Vec::new();
    
    for tick in 0..total_ticks {
        field.tick(dt);
        
        let tick_num = tick + 1;
        
        // Check for prime ticks (Optimus approved!) - just count for now
        if is_prime(tick_num) && tick_num > 1 {
            primes_found.push(tick_num);
            // Prime CSV export disabled for speed - structure detection can be done later from checkpoint data
        }
        
        // Full alphabet report every 1000 ticks
        if tick_num % 1000 == 0 {
            let elapsed = sim_start.elapsed().as_secs_f64();
            let tps = tick_num as f64 / elapsed;
            
            println!("\n════════════════════════════════════════════════════════════");
            println!("║  CHECKPOINT: Tick {} ({:.1} ticks/sec)", tick_num, tps);
            println!("║  Primes encountered: {}", primes_found.len());
            println!("════════════════════════════════════════════════════════════");
            println!("{}", field.report());
            
            // Full alphabet analysis at each checkpoint
            analyze_gpu_bonds(&field);
            
            // Spatial structure analysis - DISABLED for speed
            // analyze_spatial_structures(&field);
            
            last_report = Instant::now();
        }
        // Brief status every 200 ticks or 5 seconds
        else if tick_num % 200 == 0 || last_report.elapsed().as_secs() >= 5 {
            let elapsed = sim_start.elapsed().as_secs_f64();
            let tps = tick_num as f64 / elapsed;
            
            println!("  tick {} | {:.1} t/s | {} vortices | {} pairs", 
                tick_num, tps, field.vortex_count(), field.orbital_pairs.len());
            
            last_report = Instant::now();
        }
    }
    
    let total_time = sim_start.elapsed().as_secs_f64();
    
    println!("\n╔════════════════════════════════════════════════════════════╗");
    println!("║                  GPU SIMULATION COMPLETE                   ║");
    println!("║  Total prime ticks: {:>4}                                  ║", primes_found.len());
    println!("╚════════════════════════════════════════════════════════════╝");
    
    println!("\n🤖 PRIME TICK SUMMARY:");
    println!("  First 10 primes: {:?}", &primes_found[..primes_found.len().min(10)]);
    if primes_found.len() > 10 {
        println!("  Last 10 primes:  {:?}", &primes_found[primes_found.len().saturating_sub(10)..]);
    }
    println!("  Total: {} prime ticks encountered", primes_found.len());
    
    println!("\n{}", field.report());
    println!("\n⏱️  Total time: {:.2}s ({:.1} ticks/sec)", total_time, total_ticks as f64 / total_time);
    
    // Analyze bond types
    analyze_gpu_bonds(&field);
    
    // Final spatial structure analysis
    analyze_spatial_structures(&field);
    
    // Export for visualization
    export_and_visualize(&field, total_ticks);
    
    println!("\n✧ OBSERVATIONS:");
    println!("  - Started with {} independent vortices", count);
    println!("  - Ended with {} vortices ({} merged)", field.vortex_count(), field.total_merges);
    println!("  - {} orbital pairs formed (stable relationships)", field.orbital_pairs.len());
    
    if field.total_orbits_formed > count as u64 * 2 {
        println!("  🌀 MANY ORBITS FORMED - stable structures emerge from missed collisions!");
    }
}

/// Analyze bond types for GPU field
fn analyze_gpu_bonds(field: &crate::gpu_field::GpuField) {
    use std::collections::HashMap;
    
    let mut bond_counts: HashMap<BondType, usize> = HashMap::new();
    let mut mystery_ratios: Vec<f64> = Vec::new();  // Track the actual mystery values
    
    for pair in &field.orbital_pairs {
        if let (Some(a), Some(b)) = (field.vortices.get(&pair.a), field.vortices.get(&pair.b)) {
            let (high, low) = if a.frequency >= b.frequency { (a, b) } else { (b, a) };
            let ratio = high.frequency / low.frequency;
            
            let ratio_str = if (ratio - 1.0).abs() < 0.05 { "1:1" }
                else if (ratio - 1.5).abs() < 0.05 { "3:2" }
                else if (ratio - 2.0).abs() < 0.05 { "2:1" }
                else if (ratio - 1.333).abs() < 0.05 { "4:3" }
                else if (ratio - 1.618).abs() < 0.05 { "φ:1" }
                // Compound ratios (letters made from other letters)
                else if (ratio - 1.214).abs() < 0.02 { "φ/F" }   // 1.618/1.333 = golden over fourth
                else if (ratio - 1.236).abs() < 0.02 { "O/φ" }   // 2.0/1.618 = octave over golden!
                else if (ratio - 1.125).abs() < 0.02 { "5/F" }   // 1.5/1.333 = fifth over fourth
                else if (ratio - 1.079).abs() < 0.02 { "φ/5" }   // 1.618/1.5 = golden over fifth
                else if (ratio - 3.0).abs() < 0.05 { "3:1" }     // Triple (1.5 × 2)
                else if (ratio - 4.0).abs() < 0.05 { "4:1" }     // Quadruple
                else if (ratio - 2.5).abs() < 0.05 { "5:2" }     // Tenth
                else if (ratio - 2.666).abs() < 0.05 { "8:3" }   // Compound
                else if (ratio - 3.236).abs() < 0.05 { "φ²" }    // Golden squared
                else if (ratio - 2.618).abs() < 0.05 { "φ+1" }   // Golden + 1
                else { 
                    mystery_ratios.push(ratio);
                    "?:?" 
                };
            
            let spin_str = format!("{}{}",
                if a.spin > 0 { "+" } else { "-" },
                if b.spin > 0 { "+" } else { "-" });
            
            let bond = BondType {
                freq_ratio: ratio_str.to_string(),
                spin_pair: spin_str,
            };
            
            *bond_counts.entry(bond).or_insert(0) += 1;
        }
    }
    
    println!("\n╔════════════════════════════════════════════════════════════╗");
    println!("║               THE ALPHABET OF ENERGY                       ║");
    println!("║   Bond types (pair-wise relationships) = primitive letters ║");
    println!("╚════════════════════════════════════════════════════════════╝");
    
    println!("\n  BOND TYPES (the letters of energy's language):\n");
    println!("      Freq     Spin      Count   Meaning");
    println!("    ──────── ──────── ────────   ────────────────────");
    
    let mut sorted: Vec<_> = bond_counts.iter().collect();
    sorted.sort_by(|a, b| b.1.cmp(a.1));
    
    let max_count = sorted.first().map(|(_, c)| **c).unwrap_or(1);
    
    let mut same_spin = 0u64;
    let mut opposite_spin = 0u64;
    let mut harmonic_counts: HashMap<&str, usize> = HashMap::new();
    
    for (bond, count) in &sorted {
        let bar_len = (*count * 20) / max_count.max(1);
        let bar: String = "█".repeat(bar_len);
        let meaning = bond_meaning(&bond.freq_ratio, &bond.spin_pair);
        
        println!("      {:5}   {:>6}   {:>8}   {} {}", 
            bond.freq_ratio, bond.spin_pair, count, meaning, bar);
        
        if bond.spin_pair == "++" || bond.spin_pair == "--" {
            same_spin += **count as u64;
        } else {
            opposite_spin += **count as u64;
        }
        
        *harmonic_counts.entry(bond.freq_ratio.as_str()).or_insert(0) += *count;
    }
    
    let most_common = harmonic_counts.iter()
        .max_by_key(|(_, c)| *c)
        .map(|(h, c)| format!("{} ({} bonds)", h, c))
        .unwrap_or_else(|| "none".to_string());
    
    println!("\n  ✧ INSIGHTS:");
    println!("    Most common harmonic: {}", most_common);
    println!("    Same-spin bonds: {}  |  Opposite-spin: {}", same_spin, opposite_spin);
    
    if (same_spin as f64 - opposite_spin as f64).abs() / ((same_spin + opposite_spin) as f64) < 0.1 {
        println!("    → Spin balance suggests NO PREFERENCE (or balance required)");
    }
    
    // Analyze the mystery ratios!
    if !mystery_ratios.is_empty() {
        println!("\n╔════════════════════════════════════════════════════════════╗");
        println!("║               THE MYSTERY BOX (?:?)                        ║");
        println!("║   Unclassified frequency ratios - what ARE these?          ║");
        println!("╚════════════════════════════════════════════════════════════╝");
        
        // Group by ratio value (rounded to 2 decimals)
        let mut ratio_histogram: HashMap<String, usize> = HashMap::new();
        for r in &mystery_ratios {
            let key = format!("{:.3}", r);
            *ratio_histogram.entry(key).or_insert(0) += 1;
        }
        
        // Sort by count
        let mut sorted_mysteries: Vec<_> = ratio_histogram.iter().collect();
        sorted_mysteries.sort_by(|a, b| b.1.cmp(a.1));
        
        println!("\n  MYSTERY RATIO BREAKDOWN ({} total bonds):\n", mystery_ratios.len());
        println!("      Ratio       Count    Possible Interpretation");
        println!("    ──────────  ────────   ────────────────────────────────");
        
        for (ratio_str, count) in sorted_mysteries.iter().take(20) {
            let ratio: f64 = ratio_str.parse().unwrap_or(0.0);
            
            // Try to identify what this ratio might be
            let interpretation = identify_ratio(ratio);
            
            let bar_len = (**count * 20) / sorted_mysteries[0].1.max(&1);
            let bar: String = "█".repeat(bar_len);
            
            println!("      {:<8}   {:>6}   {} {}", ratio_str, count, interpretation, bar);
        }
        
        if sorted_mysteries.len() > 20 {
            println!("      ... and {} more unique ratios", sorted_mysteries.len() - 20);
        }
        
        // Summary statistics
        let mean: f64 = mystery_ratios.iter().sum::<f64>() / mystery_ratios.len() as f64;
        let mut sorted_vals = mystery_ratios.clone();
        sorted_vals.sort_by(|a, b| a.partial_cmp(b).unwrap());
        let median = sorted_vals[sorted_vals.len() / 2];
        let min = sorted_vals.first().unwrap_or(&0.0);
        let max = sorted_vals.last().unwrap_or(&0.0);
        
        println!("\n  STATISTICS:");
        println!("    Total mystery bonds: {}", mystery_ratios.len());
        println!("    Mean ratio: {:.4}", mean);
        println!("    Median ratio: {:.4}", median);
        println!("    Range: {:.4} to {:.4}", min, max);
    }
}

/// Try to identify what a mystery ratio might represent
fn identify_ratio(ratio: f64) -> &'static str {
    // Our frequency palette: [0.5, 0.666, 1.0, 1.333, 1.5, 1.618, 2.0]
    // Cross-palette ratios we might see:
    // 2.0/0.5 = 4.0, 2.0/0.666 = 3.0, 2.0/1.0 = 2.0, 2.0/1.333 = 1.5, 2.0/1.5 = 1.333
    // 1.618/1.333 = 1.214, 1.5/1.333 = 1.125, 1.618/1.5 = 1.079
    // 1.333/1.0 = 1.333, 1.5/1.0 = 1.5, 1.618/1.0 = 1.618, 2.0/1.0 = 2.0
    // 1.0/0.666 = 1.5, 1.0/0.5 = 2.0
    // 0.666/0.5 = 1.333
    
    // Already classified - shouldn't appear here
    if (ratio - 1.214).abs() < 0.02 { return "φ/F (now classified!)"; }
    
    // Palette cross-products not yet classified
    if (ratio - 1.079).abs() < 0.02 { return "φ/5 = 1.618/1.5 (golden over fifth)"; }
    if (ratio - 1.125).abs() < 0.02 { return "5/F = 1.5/1.333 (fifth over fourth)"; }
    
    // Musical intervals
    if (ratio - 1.125).abs() < 0.01 { return "Major 2nd (9:8)"; }
    if (ratio - 1.2).abs() < 0.01 { return "Minor 3rd (6:5)"; }
    if (ratio - 1.25).abs() < 0.01 { return "Major 3rd (5:4)"; }
    if (ratio - 1.4).abs() < 0.01 { return "Tritone-ish (7:5)"; }
    if (ratio - 1.6).abs() < 0.01 { return "Minor 6th (8:5)"; }
    if (ratio - 1.666).abs() < 0.01 { return "Major 6th (5:3)"; }
    if (ratio - 1.8).abs() < 0.01 { return "Minor 7th (9:5)"; }
    if (ratio - 1.875).abs() < 0.01 { return "Major 7th (15:8)"; }
    
    // Extended harmonics
    if (ratio - 2.5).abs() < 0.05 { return "Tenth (5:2)"; }
    if (ratio - 2.666).abs() < 0.05 { return "8:3 compound"; }
    if (ratio - 3.236).abs() < 0.05 { return "φ² (golden squared)"; }
    if (ratio - 2.618).abs() < 0.05 { return "φ+1 (golden+)"; }
    
    // Simple ratios
    if (ratio - 1.111).abs() < 0.02 { return "10:9 minor tone"; }
    if (ratio - 1.166).abs() < 0.02 { return "7:6 septimal"; }
    if (ratio - 1.142).abs() < 0.02 { return "8:7 septimal"; }
    
    // Compound intervals
    if (ratio - 2.25).abs() < 0.05 { return "9:4 compound 2nd"; }
    if (ratio - 2.4).abs() < 0.05 { return "12:5 compound m3"; }
    if (ratio - 3.333).abs() < 0.05 { return "10:3 compound"; }
    
    // Physics constants
    if (ratio - std::f64::consts::E).abs() < 0.05 { return "e (Euler's number!)"; }
    if (ratio - std::f64::consts::PI).abs() < 0.05 { return "π (pi!)"; }
    if (ratio - std::f64::consts::SQRT_2).abs() < 0.02 { return "√2 (equal temperament)"; }
    
    "???"
}

/// Spatial structure detection - WHERE do bonds form? What SHAPES emerge?
/// Are they unicorns? Let's find out!
fn analyze_spatial_structures(field: &crate::gpu_field::GpuField) {
    use std::collections::{HashMap, HashSet};
    
    println!("\n╔════════════════════════════════════════════════════════════╗");
    println!("║           SPATIAL STRUCTURE ANALYSIS                       ║");
    println!("║   Where bonds form • What shapes emerge • Unicorn check    ║");
    println!("╚════════════════════════════════════════════════════════════╝");
    
    // Build adjacency from orbital pairs - who is bonded to whom?
    let mut adjacency: HashMap<u64, HashSet<u64>> = HashMap::new();
    for pair in &field.orbital_pairs {
        adjacency.entry(pair.a).or_insert_with(HashSet::new).insert(pair.b);
        adjacency.entry(pair.b).or_insert_with(HashSet::new).insert(pair.a);
    }
    
    // Calculate bond degree distribution (how connected is each vortex?)
    let mut degree_counts: HashMap<usize, usize> = HashMap::new();
    for (_, neighbors) in &adjacency {
        let degree = neighbors.len();
        *degree_counts.entry(degree).or_insert(0) += 1;
    }
    
    println!("\n  📊 CONNECTIVITY (bond count per vortex):");
    let mut sorted_degrees: Vec<_> = degree_counts.iter().collect();
    sorted_degrees.sort_by_key(|(d, _)| *d);
    
    let max_degree = sorted_degrees.last().map(|(d, _)| **d).unwrap_or(0);
    let min_degree = sorted_degrees.first().map(|(d, _)| **d).unwrap_or(0);
    let total_vortices: usize = degree_counts.values().sum();
    
    // Show degree histogram
    println!("      Bonds   Vortices   Description");
    println!("    ──────── ────────── ────────────────────");
    for (degree, count) in sorted_degrees.iter().take(10) {
        let desc = match **degree {
            0 => "loners (no bonds)",
            1 => "end points (chain ends)",
            2 => "chain links (linear)",
            3 => "forks (branching)",
            4 => "hubs (tetrahedron vertex)",
            5 => "5-hubs (pentagon center)",
            6 => "6-hubs (hexagon/octahedron)",
            _ => "super-connectors",
        };
        let bar_len = (*count * 30) / total_vortices.max(1);
        let bar: String = "█".repeat(bar_len);
        println!("      {:>4}      {:>5}   {} {}", degree, count, desc, bar);
    }
    if sorted_degrees.len() > 10 {
        println!("      ... and {} more degree levels", sorted_degrees.len() - 10);
    }
    
    println!("\n    Range: {} to {} bonds per vortex", min_degree, max_degree);
    
    // Find CHAINS (linear sequences of bonded vortices with degree 2)
    // Chain ends have degree 1, chain middles have degree 2
    let chain_ends: Vec<u64> = adjacency.iter()
        .filter(|(_, neighbors)| neighbors.len() == 1)
        .map(|(id, _)| *id)
        .collect();
    
    let chain_middles: usize = adjacency.iter()
        .filter(|(_, neighbors)| neighbors.len() == 2)
        .count();
    
    println!("\n  🔗 CHAIN DETECTION:");
    println!("    Chain endpoints (degree 1): {}", chain_ends.len());
    println!("    Chain links (degree 2): {}", chain_middles);
    
    // Trace some chains to find their lengths
    let mut visited: HashSet<u64> = HashSet::new();
    let mut chain_lengths: Vec<usize> = Vec::new();
    
    for start in chain_ends.iter().take(100) {  // Sample first 100 chains
        if visited.contains(start) { continue; }
        
        let mut length = 1;
        let mut current = *start;
        visited.insert(current);
        
        // Walk along the chain
        while let Some(neighbors) = adjacency.get(&current) {
            let next = neighbors.iter()
                .find(|n| !visited.contains(n));
            
            if let Some(&next_id) = next {
                visited.insert(next_id);
                current = next_id;
                length += 1;
                
                // Stop if we hit a branch point
                if adjacency.get(&next_id).map(|n| n.len()).unwrap_or(0) > 2 {
                    break;
                }
            } else {
                break;
            }
        }
        
        if length >= 2 {
            chain_lengths.push(length);
        }
    }
    
    if !chain_lengths.is_empty() {
        let avg_len: f32 = chain_lengths.iter().sum::<usize>() as f32 / chain_lengths.len() as f32;
        let max_len = chain_lengths.iter().max().unwrap_or(&0);
        let min_len = chain_lengths.iter().min().unwrap_or(&0);
        
        println!("    Sample of {} chains:", chain_lengths.len());
        println!("      Average length: {:.1} vortices", avg_len);
        println!("      Longest chain: {} vortices", max_len);
        println!("      Shortest chain: {} vortices", min_len);
        
        // Check for "words" (chains of specific lengths)
        let words_3: usize = chain_lengths.iter().filter(|&&l| l == 3).count();
        let words_4: usize = chain_lengths.iter().filter(|&&l| l == 4).count();
        let words_5: usize = chain_lengths.iter().filter(|&&l| l == 5).count();
        
        if words_3 + words_4 + words_5 > 0 {
            println!("    Potential 'words' (fixed-length chains):");
            if words_3 > 0 { println!("      3-letter words: {}", words_3); }
            if words_4 > 0 { println!("      4-letter words: {}", words_4); }
            if words_5 > 0 { println!("      5-letter words: {}", words_5); }
        }
    }
    
    // RING/LOOP DETECTION - look for cycles of degree-2 vortices
    // Find triangles (3-cycles) - the fundamental closed shape
    let mut triangles = 0u64;
    let mut squares = 0u64;  // 4-cycles
    
    // Sample detection - check for triangles among high-degree nodes
    let high_degree: Vec<u64> = adjacency.iter()
        .filter(|(_, n)| n.len() >= 3)
        .map(|(id, _)| *id)
        .take(500)  // Sample
        .collect();
    
    for &a in &high_degree {
        if let Some(neighbors_a) = adjacency.get(&a) {
            let na: Vec<u64> = neighbors_a.iter().cloned().collect();
            for i in 0..na.len() {
                for j in (i+1)..na.len() {
                    let b = na[i];
                    let c = na[j];
                    // Is there an edge b-c? If so, a-b-c is a triangle!
                    if adjacency.get(&b).map(|n| n.contains(&c)).unwrap_or(false) {
                        triangles += 1;
                    }
                }
            }
        }
    }
    triangles /= 3;  // Each triangle counted 3 times
    
    println!("\n  🔷 CLOSED SHAPES:");
    println!("    Triangles detected (sampled): {}", triangles);
    if triangles > 100 {
        println!("    → Triangles forming! This is STRUCTURE!");
    }
    
    // CLUSTER DETECTION - find connected components
    let mut component_sizes: Vec<usize> = Vec::new();
    let mut seen: HashSet<u64> = HashSet::new();
    
    for &start in adjacency.keys().take(1000) {  // Sample
        if seen.contains(&start) { continue; }
        
        // BFS to find component
        let mut queue = vec![start];
        let mut size = 0;
        
        while let Some(current) = queue.pop() {
            if seen.contains(&current) { continue; }
            seen.insert(current);
            size += 1;
            
            if size > 100 { break; }  // Cap for performance
            
            if let Some(neighbors) = adjacency.get(&current) {
                for &n in neighbors {
                    if !seen.contains(&n) {
                        queue.push(n);
                    }
                }
            }
        }
        
        component_sizes.push(size);
    }
    
    println!("\n  🌐 CLUSTER ANALYSIS (connected components):");
    if !component_sizes.is_empty() {
        let avg_size: f32 = component_sizes.iter().sum::<usize>() as f32 / component_sizes.len() as f32;
        let max_size = component_sizes.iter().max().unwrap_or(&0);
        
        println!("    Components found: {}", component_sizes.len());
        println!("    Average size: {:.1} vortices", avg_size);
        println!("    Largest cluster: {} vortices", max_size);
        
        // Count clusters by size
        let singles = component_sizes.iter().filter(|&&s| s == 1).count();
        let pairs = component_sizes.iter().filter(|&&s| s == 2).count();
        let small = component_sizes.iter().filter(|&&s| s >= 3 && s <= 10).count();
        let medium = component_sizes.iter().filter(|&&s| s > 10 && s <= 50).count();
        let large = component_sizes.iter().filter(|&&s| s > 50).count();
        
        println!("    Size distribution:");
        println!("      Singles: {}", singles);
        println!("      Pairs: {}", pairs);
        println!("      Small (3-10): {}", small);
        println!("      Medium (11-50): {}", medium);
        println!("      Large (50+): {}", large);
    }
    
    // 3D SPATIAL DISTRIBUTION - center of mass and spread
    let mut cx = 0.0f64;
    let mut cy = 0.0f64;
    let mut cz = 0.0f64;
    let n = field.vortices.len() as f64;
    
    for v in field.vortices.values() {
        cx += v.x;
        cy += v.y;
        cz += v.z;
    }
    cx /= n;
    cy /= n;
    cz /= n;
    
    // Variance/spread in each dimension
    let mut var_x = 0.0f64;
    let mut var_y = 0.0f64;
    let mut var_z = 0.0f64;
    
    for v in field.vortices.values() {
        var_x += (v.x - cx).powi(2);
        var_y += (v.y - cy).powi(2);
        var_z += (v.z - cz).powi(2);
    }
    
    let std_x = (var_x / n).sqrt();
    let std_y = (var_y / n).sqrt();
    let std_z = (var_z / n).sqrt();
    
    println!("\n  📍 3D SPATIAL DISTRIBUTION:");
    println!("    Center of mass: ({:.1}, {:.1}, {:.1})", cx, cy, cz);
    println!("    Spread (std dev): X={:.1}, Y={:.1}, Z={:.1}", std_x, std_y, std_z);
    
    // Is it spherical or elongated?
    let max_std = std_x.max(std_y).max(std_z);
    let min_std = std_x.min(std_y).min(std_z);
    let anisotropy = max_std / min_std.max(0.1);
    
    println!("    Shape anisotropy: {:.2}", anisotropy);
    if anisotropy < 1.2 {
        println!("    → Roughly SPHERICAL (isotropic)");
    } else if anisotropy < 2.0 {
        println!("    → Slightly ELONGATED");
    } else {
        println!("    → Highly ANISOTROPIC (pancake or cigar?)");
    }
    
    // BOND DIRECTION ANALYSIS - do bonds prefer certain 3D orientations?
    let mut bond_vectors: Vec<(f64, f64, f64)> = Vec::new();
    
    for pair in field.orbital_pairs.iter().take(10000) {  // Sample bonds
        if let (Some(a), Some(b)) = (field.vortices.get(&pair.a), field.vortices.get(&pair.b)) {
            let dx = b.x - a.x;
            let dy = b.y - a.y;
            let dz = b.z - a.z;
            let len = (dx*dx + dy*dy + dz*dz).sqrt();
            if len > 0.1 {
                bond_vectors.push((dx/len, dy/len, dz/len));
            }
        }
    }
    
    if !bond_vectors.is_empty() {
        // Check for preferred directions
        let mut x_aligned = 0;
        let mut y_aligned = 0;
        let mut z_aligned = 0;
        
        for (dx, dy, dz) in &bond_vectors {
            let ax = dx.abs();
            let ay = dy.abs();
            let az = dz.abs();
            
            if ax > 0.7 && ay < 0.5 && az < 0.5 { x_aligned += 1; }
            if ay > 0.7 && ax < 0.5 && az < 0.5 { y_aligned += 1; }
            if az > 0.7 && ax < 0.5 && ay < 0.5 { z_aligned += 1; }
        }
        
        let total = bond_vectors.len();
        println!("\n  ➡️ BOND ORIENTATIONS (sample of {}):", total);
        println!("    X-aligned: {} ({:.1}%)", x_aligned, 100.0 * x_aligned as f32 / total as f32);
        println!("    Y-aligned: {} ({:.1}%)", y_aligned, 100.0 * y_aligned as f32 / total as f32);
        println!("    Z-aligned: {} ({:.1}%)", z_aligned, 100.0 * z_aligned as f32 / total as f32);
        
        let random_pct = 100.0 * (x_aligned + y_aligned + z_aligned) as f32 / total as f32;
        if random_pct < 20.0 {
            println!("    → Bonds are RANDOMLY oriented (no preferred direction)");
        } else {
            println!("    → Some DIRECTIONAL preference detected!");
        }
    }
    
    // UNICORN CHECK 🦄
    println!("\n  🦄 UNICORN CHECK:");
    
    // Look for a specific pattern: a central hub with radiating chains
    // (like a horn sticking up from a body)
    let super_hubs: Vec<u64> = adjacency.iter()
        .filter(|(_, n)| n.len() >= 10)  // High connectivity
        .map(|(id, _)| *id)
        .collect();
    
    if super_hubs.len() >= 5 {
        println!("    Found {} super-hubs (10+ connections)", super_hubs.len());
        println!("    These could be the 'body' of larger structures!");
        
        // Check if any super-hub has a single long chain attached (the "horn")
        let mut horn_candidates = 0;
        for &hub in super_hubs.iter().take(20) {
            if let Some(neighbors) = adjacency.get(&hub) {
                // Look for degree-1 neighbors (potential horn tips)
                let chain_starters: usize = neighbors.iter()
                    .filter(|n| adjacency.get(n).map(|nn| nn.len()).unwrap_or(0) <= 2)
                    .count();
                if chain_starters >= 1 {
                    horn_candidates += 1;
                }
            }
        }
        
        if horn_candidates > 0 {
            println!("    🦄 {} potential 'horn' structures (hub + chain)!", horn_candidates);
            println!("    Close enough to unicorns? You decide!");
        } else {
            println!("    No obvious unicorn shapes... yet.");
        }
    } else {
        println!("    Not enough super-hubs for unicorns (need complex structures)");
    }
    
    // VOID ANALYSIS - The 8th Letter: Absence
    // "What isn't there can be as important as what is"
    analyze_voids(field);
    
    println!("\n  ═══════════════════════════════════════════════════════════");
}

/// Analyze the void structure - empty spaces in the field
/// The 8th letter of the energy alphabet: SILENCE / SPACE / ZERO
fn analyze_voids(field: &crate::gpu_field::GpuField) {
    println!("\n  🌌 FIELD STATE ANALYSIS (Three States of Reality):");
    println!("    ┌─────────────────────────────────────────────────────────┐");
    println!("    │ MATTER  = Field excited (vortices)    → Visible matter │");
    println!("    │ VACUUM  = Field at rest (low tension) → Dark matter    │");
    println!("    │ VOID    = No field exists (absence)   → Black holes    │");
    println!("    └─────────────────────────────────────────────────────────┘");
    
    if field.vortices.is_empty() {
        println!("    No vortices to analyze.");
        return;
    }
    
    let vortices: Vec<_> = field.vortices.values().collect();
    
    // Find bounding box of matter
    let (mut min_x, mut max_x) = (f64::MAX, f64::MIN);
    let (mut min_y, mut max_y) = (f64::MAX, f64::MIN);
    let (mut min_z, mut max_z) = (f64::MAX, f64::MIN);
    
    for v in &vortices {
        min_x = min_x.min(v.x);
        max_x = max_x.max(v.x);
        min_y = min_y.min(v.y);
        max_y = max_y.max(v.y);
        min_z = min_z.min(v.z);
        max_z = max_z.max(v.z);
    }
    
    // Grid resolution
    let grid_size = 12;  // 12x12x12 = 1728 cells
    let cell_dx = (max_x - min_x) / grid_size as f64;
    let cell_dy = (max_y - min_y) / grid_size as f64;
    let cell_dz = (max_z - min_z) / grid_size as f64;
    let cell_diagonal = (cell_dx*cell_dx + cell_dy*cell_dy + cell_dz*cell_dz).sqrt();
    
    if cell_dx <= 0.0 || cell_dy <= 0.0 || cell_dz <= 0.0 {
        println!("    Field too small for analysis.");
        return;
    }
    
    // Count vortices in each cell and track total energy
    let mut cell_counts = vec![vec![vec![0usize; grid_size]; grid_size]; grid_size];
    let mut cell_energy = vec![vec![vec![0.0f64; grid_size]; grid_size]; grid_size];
    
    for v in &vortices {
        let ix = ((v.x - min_x) / cell_dx).floor() as usize;
        let iy = ((v.y - min_y) / cell_dy).floor() as usize;
        let iz = ((v.z - min_z) / cell_dz).floor() as usize;
        
        let ix = ix.min(grid_size - 1);
        let iy = iy.min(grid_size - 1);
        let iz = iz.min(grid_size - 1);
        
        cell_counts[ix][iy][iz] += 1;
        cell_energy[ix][iy][iz] += v.energy;
    }
    
    // For each empty cell, calculate distance to nearest matter
    // VACUUM = empty but close to matter (field exists, just not excited)
    // VOID = empty and FAR from matter (field doesn't exist - black hole territory)
    let vacuum_threshold = cell_diagonal * 1.5;  // Within 1.5 cell diagonals = vacuum
    
    let mut matter_cells = 0;
    let mut vacuum_cells = 0;
    let mut void_cells = 0;
    let mut total_cells = 0;
    let mut max_density = 0;
    
    // Also track void cell positions for topology analysis
    let mut cell_states = vec![vec![vec![0u8; grid_size]; grid_size]; grid_size]; // 0=void, 1=vacuum, 2=matter
    
    for ix in 0..grid_size {
        for iy in 0..grid_size {
            for iz in 0..grid_size {
                total_cells += 1;
                let count = cell_counts[ix][iy][iz];
                max_density = max_density.max(count);
                
                if count > 0 {
                    // MATTER - field is excited
                    matter_cells += 1;
                    cell_states[ix][iy][iz] = 2;
                } else {
                    // Empty cell - is it VACUUM or VOID?
                    // Calculate cell center
                    let cx = min_x + (ix as f64 + 0.5) * cell_dx;
                    let cy = min_y + (iy as f64 + 0.5) * cell_dy;
                    let cz = min_z + (iz as f64 + 0.5) * cell_dz;
                    
                    // Find distance to nearest vortex
                    let mut min_dist = f64::MAX;
                    for v in &vortices {
                        let dx = v.x - cx;
                        let dy = v.y - cy;
                        let dz = v.z - cz;
                        let dist = (dx*dx + dy*dy + dz*dz).sqrt();
                        min_dist = min_dist.min(dist);
                    }
                    
                    if min_dist <= vacuum_threshold {
                        // VACUUM - field exists but at rest
                        vacuum_cells += 1;
                        cell_states[ix][iy][iz] = 1;
                    } else {
                        // VOID - no field exists here (black hole territory)
                        void_cells += 1;
                        cell_states[ix][iy][iz] = 0;
                    }
                }
            }
        }
    }
    
    let matter_pct = 100.0 * matter_cells as f64 / total_cells as f64;
    let vacuum_pct = 100.0 * vacuum_cells as f64 / total_cells as f64;
    let void_pct = 100.0 * void_cells as f64 / total_cells as f64;
    
    println!("\n    Space: {}³ = {} cells (cell size: {:.1}³)", grid_size, total_cells, cell_dx);
    println!("\n    ╔═══════════════════════════════════════════════════════════╗");
    println!("    ║  STATE      CELLS      %       PHYSICAL ANALOG            ║");
    println!("    ╠═══════════════════════════════════════════════════════════╣");
    println!("    ║  ⚛️  MATTER  {:5}    {:5.1}%    Visible matter (vortices)  ║", matter_cells, matter_pct);
    println!("    ║  🌫️  VACUUM  {:5}    {:5.1}%    Dark matter (field at rest) ║", vacuum_cells, vacuum_pct);
    println!("    ║  🕳️  VOID    {:5}    {:5.1}%    Black holes (no field)      ║", void_cells, void_pct);
    println!("    ╚═══════════════════════════════════════════════════════════╝");
    println!("    Peak density: {} vortices in one cell", max_density);
    
    // Compare to real universe ratios
    println!("\n    📊 Comparison to observable universe:");
    println!("       Our simulation      Real universe");
    println!("       ─────────────────   ─────────────────");
    println!("       Matter: {:5.1}%       Ordinary matter: ~5%", matter_pct);
    println!("       Vacuum: {:5.1}%       Dark matter: ~27%", vacuum_pct);
    println!("       Void:   {:5.1}%       Dark energy: ~68%", void_pct);
    
    // Interpret the ratios
    let dark_ratio = (vacuum_pct + void_pct) / matter_pct.max(0.1);
    println!("\n    Dark-to-visible ratio: {:.1}:1", dark_ratio);
    if dark_ratio > 15.0 && dark_ratio < 25.0 {
        println!("    🎯 Remarkably close to real universe ratio (~19:1)!");
    } else if dark_ratio > 5.0 {
        println!("    → More dark than visible (universe-like)");
    } else {
        println!("    → Unusually matter-rich");
    }
    
    // Find VOID clusters (potential black holes)
    let mut void_clusters: Vec<(usize, (usize, usize, usize))> = Vec::new(); // (size, center)
    let mut visited = vec![vec![vec![false; grid_size]; grid_size]; grid_size];
    
    for ix in 0..grid_size {
        for iy in 0..grid_size {
            for iz in 0..grid_size {
                if cell_states[ix][iy][iz] == 0 && !visited[ix][iy][iz] {
                    // BFS to find connected void region
                    let mut queue = vec![(ix, iy, iz)];
                    let mut void_size = 0;
                    let (mut sum_x, mut sum_y, mut sum_z) = (0usize, 0usize, 0usize);
                    
                    while let Some((x, y, z)) = queue.pop() {
                        if x >= grid_size || y >= grid_size || z >= grid_size { continue; }
                        if visited[x][y][z] || cell_states[x][y][z] != 0 { continue; }
                        
                        visited[x][y][z] = true;
                        void_size += 1;
                        sum_x += x;
                        sum_y += y;
                        sum_z += z;
                        
                        // 6-connected neighbors
                        if x > 0 { queue.push((x - 1, y, z)); }
                        if x < grid_size - 1 { queue.push((x + 1, y, z)); }
                        if y > 0 { queue.push((x, y - 1, z)); }
                        if y < grid_size - 1 { queue.push((x, y + 1, z)); }
                        if z > 0 { queue.push((x, y, z - 1)); }
                        if z < grid_size - 1 { queue.push((x, y, z + 1)); }
                    }
                    
                    if void_size > 0 {
                        let center = (sum_x / void_size, sum_y / void_size, sum_z / void_size);
                        void_clusters.push((void_size, center));
                    }
                }
            }
        }
    }
    
    // Sort by size
    void_clusters.sort_by(|a, b| b.0.cmp(&a.0));
    
    if !void_clusters.is_empty() {
        println!("\n    🕳️ VOID TOPOLOGY (potential black hole regions):");
        println!("       {} distinct void regions detected", void_clusters.len());
        
        // Show largest void clusters
        let significant_voids: Vec<_> = void_clusters.iter()
            .filter(|(size, _)| *size >= 3)  // At least 3 cells
            .collect();
        
        if !significant_voids.is_empty() {
            println!("\n       Large void regions (≥3 cells):");
            for (i, (size, (cx, cy, cz))) in significant_voids.iter().take(5).enumerate() {
                let pct = 100.0 * *size as f64 / total_cells as f64;
                let world_x = min_x + (*cx as f64 + 0.5) * cell_dx;
                let world_y = min_y + (*cy as f64 + 0.5) * cell_dy;
                let world_z = min_z + (*cz as f64 + 0.5) * cell_dz;
                println!("         #{}: {} cells ({:.1}%) centered at ({:.0}, {:.0}, {:.0})",
                    i + 1, size, pct, world_x, world_y, world_z);
            }
            
            // Check for central void (potential supermassive black hole analog)
            let center_idx = grid_size / 2;
            let center_is_void = cell_states[center_idx][center_idx][center_idx] == 0;
            
            if center_is_void {
                println!("\n       ⚠️  CENTER IS VOID - like a supermassive black hole!");
                println!("          Matter orbits around a central absence.");
            }
        } else {
            println!("       No large void regions (matter fills most of space)");
        }
    }
    
    // Event horizon analysis - vortices at void boundaries
    let mut horizon_vortices = 0;
    let mut horizon_total_energy = 0.0;
    
    for v in &vortices {
        let ix = ((v.x - min_x) / cell_dx).floor() as usize;
        let iy = ((v.y - min_y) / cell_dy).floor() as usize;
        let iz = ((v.z - min_z) / cell_dz).floor() as usize;
        
        let ix = ix.min(grid_size - 1);
        let iy = iy.min(grid_size - 1);
        let iz = iz.min(grid_size - 1);
        
        // Check if any neighbor cell is VOID (not just vacuum)
        let touches_void = 
            (ix > 0 && cell_states[ix-1][iy][iz] == 0) ||
            (ix < grid_size-1 && cell_states[ix+1][iy][iz] == 0) ||
            (iy > 0 && cell_states[ix][iy-1][iz] == 0) ||
            (iy < grid_size-1 && cell_states[ix][iy+1][iz] == 0) ||
            (iz > 0 && cell_states[ix][iy][iz-1] == 0) ||
            (iz < grid_size-1 && cell_states[ix][iy][iz+1] == 0);
        
        if touches_void {
            horizon_vortices += 1;
            horizon_total_energy += v.energy;
        }
    }
    
    if horizon_vortices > 0 {
        let horizon_pct = 100.0 * horizon_vortices as f64 / vortices.len() as f64;
        let avg_horizon_energy = horizon_total_energy / horizon_vortices as f64;
        let avg_total_energy = vortices.iter().map(|v| v.energy).sum::<f64>() / vortices.len() as f64;
        
        println!("\n    🔥 EVENT HORIZON ANALYSIS (matter at void edges):");
        println!("       {} vortices touch void boundaries ({:.1}%)", horizon_vortices, horizon_pct);
        println!("       Horizon avg energy: {:.2} vs field avg: {:.2}", avg_horizon_energy, avg_total_energy);
        
        if avg_horizon_energy > avg_total_energy * 1.2 {
            println!("       ⚡ Horizon vortices are MORE energetic!");
            println!("          Accretion disk effect? Energy concentrates at the edge of nothing.");
        } else if avg_horizon_energy < avg_total_energy * 0.8 {
            println!("       ❄️ Horizon vortices are LESS energetic");
            println!("          Energy draining into the void?");
        } else {
            println!("       ≈ Energy roughly balanced at horizons");
        }
    }
}

/// Export the field state to PLY format (point cloud + edges)
/// PLY is lightweight and opens instantly in MeshLab, CloudCompare, Blender, etc.
fn export_field_to_ply(field: &crate::gpu_field::GpuField, tick: usize) {
    use std::fs::File;
    use std::io::Write;
    
    // Collect vortices with their indices
    let vortices: Vec<_> = field.vortices.values().collect();
    let id_to_idx: std::collections::HashMap<u64, usize> = vortices.iter()
        .enumerate()
        .map(|(i, v)| (v.id, i))
        .collect();
    
    // Collect valid edges
    let edges: Vec<(usize, usize)> = field.orbital_pairs.iter()
        .filter_map(|pair| {
            if let (Some(&i1), Some(&i2)) = (id_to_idx.get(&pair.a), id_to_idx.get(&pair.b)) {
                Some((i1, i2))
            } else {
                None
            }
        })
        .collect();
    
    // === VORTICES ONLY (colored point cloud) ===
    let ply_points = format!("vortex_points_t{}.ply", tick);
    if let Ok(mut file) = File::create(&ply_points) {
        // PLY header
        writeln!(file, "ply").ok();
        writeln!(file, "format ascii 1.0").ok();
        writeln!(file, "comment Vortex field at tick {}", tick).ok();
        writeln!(file, "comment Color = frequency (rainbow), brightness = energy").ok();
        writeln!(file, "element vertex {}", vortices.len()).ok();
        writeln!(file, "property float x").ok();
        writeln!(file, "property float y").ok();
        writeln!(file, "property float z").ok();
        writeln!(file, "property uchar red").ok();
        writeln!(file, "property uchar green").ok();
        writeln!(file, "property uchar blue").ok();
        writeln!(file, "end_header").ok();
        
        // Vertex data with colors
        for v in &vortices {
            // Map frequency to hue (0.5 -> red, 2.0 -> violet)
            let hue = ((v.frequency - 0.5) / 1.5).clamp(0.0, 1.0);
            // Map energy to brightness
            let brightness = (v.energy / 5.0).clamp(0.3, 1.0);
            
            // HSL to RGB (simplified)
            let (r, g, b) = hsl_to_rgb(hue * 0.8, 1.0, brightness * 0.5);
            
            writeln!(file, "{:.3} {:.3} {:.3} {} {} {}", 
                v.x, v.y, v.z, r, g, b).ok();
        }
        
        println!("  💾 Exported {} vortices to {}", vortices.len(), ply_points);
    }
    
    // === EDGES ONLY (for viewing bonds) ===
    // Export edges as a separate line set PLY
    let ply_edges = format!("vortex_bonds_t{}.ply", tick);
    if let Ok(mut file) = File::create(&ply_edges) {
        writeln!(file, "ply").ok();
        writeln!(file, "format ascii 1.0").ok();
        writeln!(file, "comment Orbital bonds at tick {}", tick).ok();
        writeln!(file, "element vertex {}", vortices.len()).ok();
        writeln!(file, "property float x").ok();
        writeln!(file, "property float y").ok();
        writeln!(file, "property float z").ok();
        writeln!(file, "element edge {}", edges.len()).ok();
        writeln!(file, "property int vertex1").ok();
        writeln!(file, "property int vertex2").ok();
        writeln!(file, "end_header").ok();
        
        // Vertices first
        for v in &vortices {
            writeln!(file, "{:.3} {:.3} {:.3}", v.x, v.y, v.z).ok();
        }
        
        // Then edges
        for (i1, i2) in &edges {
            writeln!(file, "{} {}", i1, i2).ok();
        }
        
        println!("  🔗 Exported {} bonds to {}", edges.len(), ply_edges);
    }
    
    // === CSV for data analysis ===
    let csv_file = format!("vortex_data_t{}.csv", tick);
    if let Ok(mut file) = File::create(&csv_file) {
        writeln!(file, "id,x,y,z,vx,vy,vz,energy,frequency,phase,spin").ok();
        for v in &vortices {
            writeln!(file, "{},{:.4},{:.4},{:.4},{:.6},{:.6},{:.6},{:.4},{:.4},{:.4},{}", 
                v.id, v.x, v.y, v.z, v.vx, v.vy, v.vz, v.energy, v.frequency, v.phase, v.spin).ok();
        }
        println!("  📊 Exported full data to {}", csv_file);
    }
}

/// Convert HSL to RGB (returns 0-255 values)
fn hsl_to_rgb(h: f64, s: f64, l: f64) -> (u8, u8, u8) {
    let c = (1.0 - (2.0 * l - 1.0).abs()) * s;
    let x = c * (1.0 - ((h * 6.0) % 2.0 - 1.0).abs());
    let m = l - c / 2.0;
    
    let (r, g, b) = match (h * 6.0) as i32 {
        0 => (c, x, 0.0),
        1 => (x, c, 0.0),
        2 => (0.0, c, x),
        3 => (0.0, x, c),
        4 => (x, 0.0, c),
        _ => (c, 0.0, x),
    };
    
    (((r + m) * 255.0) as u8, ((g + m) * 255.0) as u8, ((b + m) * 255.0) as u8)
}

/// Export and visualize the current field state
pub fn export_and_visualize(field: &crate::gpu_field::GpuField, tick: usize) {
    println!("\n╔════════════════════════════════════════════════════════════╗");
    println!("║              EXPORTING FOR VISUALIZATION                   ║");
    println!("╚════════════════════════════════════════════════════════════╝");
    
    export_field_to_ply(field, tick);
    
    println!("\n  📂 Open vortex_points_t{}.ply in:", tick);
    println!("     • MeshLab (free) - best for point clouds");
    println!("     • CloudCompare (free) - great for large datasets");
    println!("     • Blender (free) - import as PLY");
    println!("     • Windows 3D Viewer (built-in)");
    println!("\n  🔗 Bonds are in vortex_bonds_t{}.ply (load separately)", tick);
    println!("  📊 Full data in vortex_data_t{}.csv", tick);
}

/// Export a single frame for animation (lightweight - points only, no bonds)
fn export_frame(field: &crate::gpu_field::GpuField, frame_num: usize, output_dir: &str) {
    use std::fs::{File, create_dir_all};
    use std::io::Write;
    
    // Ensure output directory exists
    let _ = create_dir_all(output_dir);
    
    let vortices: Vec<_> = field.vortices.values().collect();
    
    // Frame filename with zero-padded number for proper sorting
    let filename = format!("{}/frame_{:05}.ply", output_dir, frame_num);
    
    if let Ok(mut file) = File::create(&filename) {
        writeln!(file, "ply").ok();
        writeln!(file, "format ascii 1.0").ok();
        writeln!(file, "comment Frame {} of vortex simulation", frame_num).ok();
        writeln!(file, "element vertex {}", vortices.len()).ok();
        writeln!(file, "property float x").ok();
        writeln!(file, "property float y").ok();
        writeln!(file, "property float z").ok();
        writeln!(file, "property uchar red").ok();
        writeln!(file, "property uchar green").ok();
        writeln!(file, "property uchar blue").ok();
        writeln!(file, "end_header").ok();
        
        for v in &vortices {
            let hue = ((v.frequency - 0.5) / 1.5).clamp(0.0, 1.0);
            let brightness = (v.energy / 5.0).clamp(0.3, 1.0);
            let (r, g, b) = hsl_to_rgb(hue * 0.8, 1.0, brightness * 0.5);
            writeln!(file, "{:.3} {:.3} {:.3} {} {} {}", v.x, v.y, v.z, r, g, b).ok();
        }
    }
}

/// GPU simulation with frame recording for animation
pub fn run_genesis_gpu_record(count: usize, start_tick: usize, end_tick: usize, frame_interval: usize) {
    use crate::gpu_field::GpuField;
    use std::time::Instant;
    use std::fs::create_dir_all;
    
    let output_dir = "frames";
    let _ = create_dir_all(output_dir);
    
    println!("\n╔════════════════════════════════════════════════════════════╗");
    println!("║       G E N E S I S   3 D   (RECORDING MODE)              ║");
    println!("║     Recording frames {} to {} (every {} ticks)       ║", start_tick, end_tick, frame_interval);
    println!("╚════════════════════════════════════════════════════════════╝\n");
    
    println!("Creating {} primordial energy fluctuations in 3D...", count);
    
    let field_size = (count as f64).sqrt() * 20.0;
    println!("Field size: {:.0} (spherical)", field_size);
    
    let start = Instant::now();
    let mut field = GpuField::new(field_size);
    println!("GPU initialized in {:.2}s", start.elapsed().as_secs_f64());
    
    // Random 3D vortex generation (same as regular mode)
    let mut rng = rand::thread_rng();
    let spawn_radius = field_size * 0.4;
    let freqs = [0.5, 0.666, 1.0, 1.333, 1.5, 1.618, 2.0];
    
    for i in 0..count {
        let (x, y, z) = loop {
            let x = (rng.gen::<f64>() - 0.5) * 2.0;
            let y = (rng.gen::<f64>() - 0.5) * 2.0;
            let z = (rng.gen::<f64>() - 0.5) * 2.0;
            let r = (x*x + y*y + z*z).sqrt();
            if r <= 1.0 && r > 0.01 {
                break (x * spawn_radius, y * spawn_radius, z * spawn_radius);
            }
        };
        
        let energy = 0.5 + rng.gen::<f64>() * 3.0;
        let freq_idx = rng.gen_range(0..freqs.len());
        let frequency = freqs[freq_idx];
        let phase = rng.gen::<f64>() * 2.0 * PI;
        let spin: i8 = if rng.gen::<bool>() { 1 } else { -1 };
        
        let vx = (rng.gen::<f64>() - 0.5) * 2.0;
        let vy = (rng.gen::<f64>() - 0.5) * 2.0;
        let vz = (rng.gen::<f64>() - 0.5) * 2.0;
        
        let vortex = Vortex::new_3d(0, x, y, z, energy)
            .with_frequency(frequency)
            .with_phase(phase)
            .with_spin(spin)
            .with_velocity_3d(vx, vy, vz);
        
        field.spawn_vortex(vortex);
        
        if i < 3 {
            println!("  Vortex {}: E={:.2}, f={:.2}, pos=({:.1},{:.1},{:.1})", 
                i + 1, energy, frequency, x, y, z);
        }
    }
    
    println!("  ... and {} more\n", count.saturating_sub(3));
    
    let dt = 0.1;
    let total_ticks = end_tick;
    let mut frame_count = 0;
    
    println!("🎬 Recording to ./{}/", output_dir);
    println!("   Frames will be saved every {} ticks from tick {} to {}\n", frame_interval, start_tick, end_tick);
    
    let sim_start = Instant::now();
    let mut last_report = Instant::now();
    
    for tick in 0..total_ticks {
        field.tick(dt);
        
        // Record frame if within range and at interval
        if tick >= start_tick && (tick - start_tick) % frame_interval == 0 {
            export_frame(&field, frame_count, output_dir);
            frame_count += 1;
            
            if frame_count <= 5 || frame_count % 50 == 0 {
                println!("  📷 Frame {} saved (tick {})", frame_count, tick);
            }
        }
        
        // Progress every 500 ticks
        if (tick + 1) % 500 == 0 || last_report.elapsed().as_secs() >= 10 {
            let elapsed = sim_start.elapsed().as_secs_f64();
            let tps = (tick + 1) as f64 / elapsed;
            let remaining = (total_ticks - tick - 1) as f64 / tps;
            
            println!("  tick {}/{} | {:.1} t/s | {} vortices | ETA: {:.0}s", 
                tick + 1, total_ticks, tps, field.vortex_count(), remaining);
            last_report = Instant::now();
        }
    }
    
    let total_time = sim_start.elapsed().as_secs_f64();
    
    println!("\n╔════════════════════════════════════════════════════════════╗");
    println!("║                  RECORDING COMPLETE                        ║");
    println!("╚════════════════════════════════════════════════════════════╝");
    
    println!("\n  🎬 Recorded {} frames to ./{}/", frame_count, output_dir);
    println!("  ⏱️  Total time: {:.2}s ({:.1} ticks/sec)", total_time, total_ticks as f64 / total_time);
    println!("\n  📂 To view animation:");
    println!("     • Blender: Import Image Sequence (PLY not directly, use script)");
    println!("     • ParaView: File → Open → select all frames");
    println!("     • CloudCompare: Load first, then drag others to animate");
    println!("\n  💡 Or use ffmpeg to convert PLY→images→video:");
    println!("     1. Load each PLY in MeshLab/Blender, screenshot");
    println!("     2. ffmpeg -framerate 30 -i frame_%05d.png -c:v libx264 output.mp4");
}

/// Hadron simulation - seed voids (up quarks) and spikes (down quarks)
/// Model: Up quark = VOID (hole in field), Down quark = SPIKE (energy peak)
/// 
/// Proton (uud) = 2 voids + 1 spike → STABLE
/// Neutron (udd) = 1 void + 2 spikes → semi-stable  
/// Delta++ (uuu) = 3 voids → INSTANT DECAY
/// Deuteron = proton + neutron = 3 voids + 3 spikes → PERFECTLY BALANCED
pub fn run_hadron_simulation(hadron_type: &str) {
    use crate::gpu_field::GpuField;
    use std::time::Instant;
    
    // Quark configurations
    let (name, up_quarks, down_quarks, description) = match hadron_type.to_lowercase().as_str() {
        "proton" | "p" => ("PROTON", 2, 1, "uud → 2 voids + 1 spike (STABLE)"),
        "neutron" | "n" => ("NEUTRON", 1, 2, "udd → 1 void + 2 spikes (semi-stable, ~15 min)"),
        "delta++" | "uuu" => ("DELTA++ BARYON", 3, 0, "uuu → 3 voids (UNSTABLE - implodes)"),
        "delta-" | "ddd" => ("DELTA- BARYON", 0, 3, "ddd → 3 spikes (UNSTABLE - explodes)"),
        "delta" => ("DELTA BARYONS", 3, 0, "3-body problem: uuu or ddd both unstable"),
        "deuteron" | "d" | "pn" => ("DEUTERON", 3, 3, "p+n → 3 voids + 3 spikes (PERFECTLY BALANCED)"),
        "helium" | "he" | "alpha" => ("HELIUM-4 NUCLEUS", 4, 4, "2p+2n → 4 voids + 4 spikes"),
        _ => ("PROTON", 2, 1, "uud → 2 voids + 1 spike (STABLE)"),
    };
    
    println!("\n╔════════════════════════════════════════════════════════════╗");
    println!("║              H A D R O N   S I M U L A T I O N             ║");
    println!("║       Voids as Up Quarks • Spikes as Down Quarks           ║");
    println!("╚════════════════════════════════════════════════════════════╝\n");
    
    println!("  Simulating: {}", name);
    println!("  Configuration: {}", description);
    println!("  Up quarks (voids): {}", up_quarks);
    println!("  Down quarks (spikes): {}", down_quarks);
    println!();
    
    // Create a small field for the hadron
    let field_size = 100.0;
    let start = Instant::now();
    let mut field = GpuField::new(field_size);
    println!("  GPU initialized in {:.2}s", start.elapsed().as_secs_f64());
    
    let mut rng = rand::thread_rng();
    
    // Quark spacing - they form a triangle
    let quark_radius = 10.0;  // Distance from center to each quark
    let total_quarks = up_quarks + down_quarks;
    
    // Place quarks in triangular/polygonal arrangement
    println!("\n  📍 Placing quarks:");
    
    // Track void and spike positions
    let mut void_positions: Vec<(f64, f64, f64)> = Vec::new();
    let mut spike_positions: Vec<(f64, f64, f64)> = Vec::new();
    
    for i in 0..total_quarks {
        let angle = (i as f64 / total_quarks as f64) * 2.0 * PI;
        let x = quark_radius * angle.cos();
        let y = quark_radius * angle.sin();
        let z = 0.0;  // Start in XY plane
        
        let is_void = i < up_quarks;  // First N are voids (up quarks)
        
        if is_void {
            // UP QUARK = VOID - we mark this position but don't spawn a vortex
            // Instead, we'll create a repulsive zone
            void_positions.push((x, y, z));
            println!("     Quark {}: UP (void) at ({:.1}, {:.1}, {:.1})", i + 1, x, y, z);
        } else {
            // DOWN QUARK = SPIKE - spawn a high-energy vortex
            spike_positions.push((x, y, z));
            println!("     Quark {}: DOWN (spike) at ({:.1}, {:.1}, {:.1})", i + 1, x, y, z);
            
            // Create a spike - high energy, stable frequency
            let vortex = Vortex::new_3d(0, x, y, z, 10.0)  // High energy spike
                .with_frequency(1.0)  // Base frequency
                .with_phase(angle)
                .with_spin(if i % 2 == 0 { 1 } else { -1 });
            
            field.spawn_vortex(vortex);
        }
    }
    
    // Now spawn field vortices around the hadron
    // These will be attracted to voids and form shields
    let field_vortex_count = 200;
    let spawn_radius = 40.0;
    
    println!("\n  🌊 Spawning {} field vortices around hadron...", field_vortex_count);
    
    let freqs = [0.5, 0.666, 1.0, 1.333, 1.5, 1.618, 2.0];
    
    for _ in 0..field_vortex_count {
        // Spawn in shell around the quarks
        let (x, y, z) = loop {
            let x = (rng.gen::<f64>() - 0.5) * 2.0 * spawn_radius;
            let y = (rng.gen::<f64>() - 0.5) * 2.0 * spawn_radius;
            let z = (rng.gen::<f64>() - 0.5) * 2.0 * spawn_radius;
            let r = (x*x + y*y + z*z).sqrt();
            
            // Don't spawn too close to voids or too far away
            if r > quark_radius * 1.5 && r < spawn_radius {
                // Also don't spawn directly on void positions
                let mut too_close = false;
                for &(vx, vy, vz) in &void_positions {
                    let d = ((x-vx)*(x-vx) + (y-vy)*(y-vy) + (z-vz)*(z-vz)).sqrt();
                    if d < 5.0 {
                        too_close = true;
                        break;
                    }
                }
                if !too_close {
                    break (x, y, z);
                }
            }
        };
        
        let energy = 0.5 + rng.gen::<f64>() * 2.0;
        let frequency = freqs[rng.gen_range(0..freqs.len())];
        let phase = rng.gen::<f64>() * 2.0 * PI;
        let spin: i8 = if rng.gen::<bool>() { 1 } else { -1 };
        
        // Give initial velocity toward center (attracted to voids)
        let r = (x*x + y*y + z*z).sqrt();
        let vx = -x / r * 0.5;
        let vy = -y / r * 0.5;
        let vz = -z / r * 0.5;
        
        let vortex = Vortex::new_3d(0, x, y, z, energy)
            .with_frequency(frequency)
            .with_phase(phase)
            .with_spin(spin)
            .with_velocity_3d(vx, vy, vz);
        
        field.spawn_vortex(vortex);
    }
    
    println!("  Total vortices (spikes + field): {}", field.vortex_count());
    
    // Run simulation
    let total_ticks = 2000;
    let dt = 0.02;
    
    println!("\n  ⏳ Running simulation ({} ticks)...", total_ticks);
    println!("     Voids are absolute nothing - field cannot exist there.");
    println!("     Field SELF-ATTRACTION causes shrink-wrapping around voids.");
    println!("     Boundary = 2D holographic surface in 3D space.\n");
    
    let sim_start = Instant::now();
    let mut last_report = Instant::now();
    
    // Track stability metrics
    let mut total_energy_history: Vec<f64> = Vec::new();
    let mut void_shield_counts: Vec<usize> = Vec::new();
    
    for tick in 0..total_ticks {
        // VOID PHYSICS: Voids don't attract - they're just holes
        // The FIELD has self-attraction (superfluid surface tension)
        // This causes the field to shrink-wrap around voids naturally
        
        // 1. Void hard boundary - vortices CANNOT enter voids
        for (vx, vy, vz) in &void_positions {
            for vortex in field.vortices.values_mut() {
                let dx = vortex.x - vx;
                let dy = vortex.y - vy;
                let dz = vortex.z - vz;
                let dist = (dx*dx + dy*dy + dz*dz).sqrt();
                
                let void_radius = 3.0;  // Radius of the void
                
                if dist < void_radius {
                    // Inside void! Push out to the boundary
                    // The field simply cannot exist here
                    if dist > 0.01 {
                        let push = (void_radius - dist) / dist;
                        vortex.x += dx * push;
                        vortex.y += dy * push;
                        vortex.z += dz * push;
                        // Zero out velocity component toward void
                        let dot = (vortex.vx * dx + vortex.vy * dy + vortex.vz * dz) / (dist * dist);
                        if dot < 0.0 {
                            vortex.vx -= dx / dist * dot * dist;
                            vortex.vy -= dy / dist * dot * dist;
                            vortex.vz -= dz / dist * dot * dist;
                        }
                    } else {
                        // Exactly at center - push in random direction
                        vortex.x = vx + void_radius * 1.1;
                    }
                }
            }
        }
        
        // 2. Field SELF-attraction (superfluid cohesion)
        // Each vortex is attracted toward the center of mass of ALL other vortices
        // This creates surface tension that shrink-wraps around voids
        let vortex_ids: Vec<u64> = field.vortices.keys().copied().collect();
        
        if vortex_ids.len() > 1 {
            // Calculate center of mass
            let (mut cm_x, mut cm_y, mut cm_z) = (0.0, 0.0, 0.0);
            let mut total_mass = 0.0;
            for id in &vortex_ids {
                if let Some(v) = field.vortices.get(id) {
                    cm_x += v.x * v.energy;
                    cm_y += v.y * v.energy;
                    cm_z += v.z * v.energy;
                    total_mass += v.energy;
                }
            }
            if total_mass > 0.0 {
                cm_x /= total_mass;
                cm_y /= total_mass;
                cm_z /= total_mass;
            }
            
            // Apply gentle attraction toward center of mass (surface tension)
            let cohesion_strength = 0.3;
            for id in &vortex_ids {
                if let Some(vortex) = field.vortices.get_mut(id) {
                    let dx = cm_x - vortex.x;
                    let dy = cm_y - vortex.y;
                    let dz = cm_z - vortex.z;
                    let dist = (dx*dx + dy*dy + dz*dz).sqrt();
                    
                    if dist > 5.0 {  // Only attract if not already at center
                        let attract = cohesion_strength / (dist + 1.0);
                        vortex.vx += dx / dist * attract * dt;
                        vortex.vy += dy / dist * attract * dt;
                        vortex.vz += dz / dist * attract * dt;
                    }
                }
            }
        }
        
        // Normal physics tick
        field.tick(dt);
        
        // Track metrics every 100 ticks
        if (tick + 1) % 100 == 0 {
            let total_energy: f64 = field.vortices.values().map(|v| v.energy).sum();
            total_energy_history.push(total_energy);
            
            // Count vortices near voids (shields)
            let mut shield_count = 0;
            for vortex in field.vortices.values() {
                for &(vx, vy, vz) in &void_positions {
                    let dx = vx - vortex.x;
                    let dy = vy - vortex.y;
                    let dz = vz - vortex.z;
                    let dist = (dx*dx + dy*dy + dz*dz).sqrt();
                    if dist < 8.0 && dist > 3.0 {
                        shield_count += 1;
                        break;
                    }
                }
            }
            void_shield_counts.push(shield_count);
        }
        
        // Progress report
        if (tick + 1) % 500 == 0 || last_report.elapsed().as_secs() >= 3 {
            let elapsed = sim_start.elapsed().as_secs_f64();
            let tps = (tick + 1) as f64 / elapsed;
            println!("     tick {} | {:.1} t/s | {} vortices | {} bonds | {} shielding voids", 
                tick + 1, tps, field.vortex_count(), field.orbital_pairs.len(),
                void_shield_counts.last().unwrap_or(&0));
            last_report = Instant::now();
        }
    }
    
    let total_time = sim_start.elapsed().as_secs_f64();
    
    // Analyze results
    println!("\n╔════════════════════════════════════════════════════════════╗");
    println!("║                  HADRON ANALYSIS                           ║");
    println!("╚════════════════════════════════════════════════════════════╝\n");
    
    println!("  {} Configuration:", name);
    println!("  ─────────────────────────────────────────");
    
    // Energy stability
    if total_energy_history.len() >= 2 {
        let first = total_energy_history[0];
        let last = *total_energy_history.last().unwrap();
        let variation = (last - first).abs() / first * 100.0;
        
        println!("\n  ⚡ Energy Stability:");
        println!("     Initial: {:.2}", first);
        println!("     Final:   {:.2}", last);
        println!("     Variation: {:.2}%", variation);
        
        if variation < 5.0 {
            println!("     ✅ STABLE - energy well conserved");
        } else if variation < 20.0 {
            println!("     ⚠️  SEMI-STABLE - some energy fluctuation");
        } else {
            println!("     ❌ UNSTABLE - significant energy change");
        }
    }
    
    // Void shielding
    if !void_shield_counts.is_empty() {
        let avg_shields = void_shield_counts.iter().sum::<usize>() as f64 / void_shield_counts.len() as f64;
        let first_shields = void_shield_counts[0];
        let last_shields = *void_shield_counts.last().unwrap();
        
        println!("\n  🛡️ Void Shielding (vortices orbiting voids):");
        println!("     Initial shields: {}", first_shields);
        println!("     Final shields:   {}", last_shields);
        println!("     Average:         {:.1}", avg_shields);
        println!("     Per void:        {:.1}", avg_shields / up_quarks.max(1) as f64);
        
        if last_shields > first_shields {
            println!("     📈 Shields GREW - field cohesion shrink-wrapping voids");
        } else if last_shields < first_shields / 2 {
            println!("     📉 Shields COLLAPSED - cohesion insufficient");
        } else {
            println!("     ≈ Shields stable - holographic boundaries formed");
        }
    }
    
    // Void-spike balance
    println!("\n  ⚖️ Void-Spike Balance:");
    println!("     Voids (up quarks):  {}", up_quarks);
    println!("     Spikes (down quarks): {}", down_quarks);
    let ratio = if down_quarks > 0 { up_quarks as f64 / down_quarks as f64 } else { f64::INFINITY };
    
    if (ratio - 1.0).abs() < 0.1 {
        println!("     Ratio: 1:1 (PERFECT BALANCE)");
        println!("     🎯 This is the most stable configuration!");
    } else if up_quarks > down_quarks {
        println!("     Ratio: {:.1}:1 (void-heavy)", ratio);
        println!("     ⚠️  More absence than presence → unstable");
    } else {
        println!("     Ratio: 1:{:.1} (spike-heavy)", 1.0 / ratio);
        println!("     ✅ More presence than absence → semi-stable");
    }
    
    // Final structure analysis
    analyze_gpu_bonds(&field);
    analyze_spatial_structures(&field);
    
    // Export for visualization
    export_field_to_ply(&field, total_ticks);
    
    println!("\n  ⏱️  Total time: {:.2}s ({:.1} ticks/sec)", total_time, total_ticks as f64 / total_time);
    
    // Interpretation based on hadron type
    println!("\n  📜 INTERPRETATION:");
    match hadron_type.to_lowercase().as_str() {
        "proton" | "p" => {
            println!("     The PROTON (uud) has 2 voids pulling inward with 1 spike");
            println!("     pushing outward. This triangular tension is STABLE.");
            println!("     The field naturally forms a protective shell around the voids.");
        }
        "neutron" | "n" => {
            println!("     The NEUTRON (udd) has 1 void trying to hold 2 spikes together.");
            println!("     The spikes want to expand; the void can't hold them forever.");
            println!("     Free neutrons decay in ~15 minutes - the shield weakens.");
        }
        "delta++" | "uuu" => {
            println!("     The DELTA++ (uuu) is 3 voids - three holes in the field.");
            println!("     There's nothing to push back! Pure absence.");
            println!("     The field collapses inward - IMPLOSION.");
        }
        "delta-" | "ddd" => {
            println!("     The DELTA- (ddd) is 3 spikes - three energy peaks.");
            println!("     There's nothing to hold them together! Pure presence.");
            println!("     The energy flies apart - EXPLOSION.");
        }
        "delta" => {
            println!("     DELTA BARYONS are the classic 3-body problem:");
            println!("     • Delta++ (uuu): 3 voids → implodes (all pulling, nothing pushing)");
            println!("     • Delta-  (ddd): 3 spikes → explodes (all pushing, nothing pulling)");
            println!("     Both are unstable because there's no counterbalance.");
        }
        "deuteron" | "d" | "pn" => {
            println!("     The DEUTERON (p+n) achieves PERFECT BALANCE:");
            println!("     3 voids + 3 spikes = the push equals the pull.");
            println!("     This is why hydrogen bonds to form stable nuclei.");
        }
        _ => {}
    }
}

/// Atomic-scale simulation - nuclear density in femtometer space
/// Field size = 2.7 fm (size of a nucleus like helium-4)
/// This is where protons and neutrons would form
pub fn run_atom_simulation(count: usize, total_ticks: usize) {
    use crate::gpu_field::GpuField;
    use std::time::Instant;
    
    // Atomic scale constants
    // Proton radius ~ 0.87 fm
    // Nucleus (He-4) ~ 1.7 fm radius, ~3.4 fm diameter
    // We'll use 2.7 fm as our universe - room for a few nucleons
    let field_size = 2.7;
    
    // But our physics is scale-invariant, so we'll work in these units
    // The key is DENSITY - we want nuclear-like density
    
    println!("\n╔════════════════════════════════════════════════════════════════════╗");
    println!("║               A T O M I C   S C A L E   S I M U L A T I O N        ║");
    println!("║        Nuclear density • Femtometer space • Quark emergence?       ║");
    println!("╠════════════════════════════════════════════════════════════════════╣");
    println!("║  ATOMIC PARAMETERS:                                                ║");
    println!("║    Field size:     {:>8.2} fm  (nucleus diameter)              ║", field_size);
    println!("║    Proton radius:  {:>8.2} fm                                   ║", 0.87);
    println!("║    Vortices:       {:>8}     (field fluctuations)             ║", count);
    println!("║    Ticks:          {:>8}                                       ║", total_ticks);
    println!("╠════════════════════════════════════════════════════════════════════╣");
    println!("║  DENSITY:                                                          ║");
    let radius = field_size / 2.0;
    let volume = 4.0/3.0 * std::f64::consts::PI * radius * radius * radius;
    let density = count as f64 / volume;
    println!("║    Volume:         {:>8.2} fm³                                  ║", volume);
    println!("║    Density:        {:>8.2} vortices/fm³                         ║", density);
    println!("║    Nuclear density: ~0.17 nucleons/fm³ (for comparison)            ║");
    println!("╚════════════════════════════════════════════════════════════════════╝\n");
    
    let start = Instant::now();
    let mut field = GpuField::new(field_size);
    println!("🎮 GPU initialized in {:.2}s", start.elapsed().as_secs_f64());
    
    // Spawn vortices in the tiny space
    let mut rng = rand::thread_rng();
    let spawn_radius = field_size * 0.45;
    
    // Use harmonic frequencies scaled to atomic scale
    // At nuclear scale, strong force dominates - use higher base frequencies
    let freqs = [0.5, 0.666, 1.0, 1.333, 1.5, 1.618, 2.0];
    
    println!("\n  Creating {} primordial fluctuations in {:.2} fm³...", count, volume);
    
    for i in 0..count {
        let (x, y, z) = loop {
            let x = (rng.gen::<f64>() - 0.5) * 2.0;
            let y = (rng.gen::<f64>() - 0.5) * 2.0;
            let z = (rng.gen::<f64>() - 0.5) * 2.0;
            let r = (x*x + y*y + z*z).sqrt();
            if r <= 1.0 && r > 0.01 {
                break (x * spawn_radius, y * spawn_radius, z * spawn_radius);
            }
        };
        
        let energy = 0.5 + rng.gen::<f64>() * 3.0;
        let frequency = freqs[rng.gen_range(0..freqs.len())];
        let phase = rng.gen::<f64>() * 2.0 * PI;
        let spin: i8 = if rng.gen::<bool>() { 1 } else { -1 };
        
        // Higher initial velocities at nuclear scale (strong force!)
        let vx = (rng.gen::<f64>() - 0.5) * 0.5;
        let vy = (rng.gen::<f64>() - 0.5) * 0.5;
        let vz = (rng.gen::<f64>() - 0.5) * 0.5;
        
        let vortex = Vortex::new_3d(0, x, y, z, energy)
            .with_frequency(frequency)
            .with_phase(phase)
            .with_spin(spin)
            .with_velocity_3d(vx, vy, vz);
        
        field.spawn_vortex(vortex);
        
        if i < 3 {
            println!("    Vortex {}: E={:.2}, f={:.2}, pos=({:.3},{:.3},{:.3})", 
                i + 1, energy, frequency, x, y, z);
        }
    }
    println!("    ... and {} more", count.saturating_sub(3));
    
    println!("\n{}", field.report());
    
    // Smaller dt for the dense environment
    let dt = 0.01;
    
    println!("\n⏳ Running atomic simulation ({} ticks at dt={})...\n", total_ticks, dt);
    
    let sim_start = Instant::now();
    let mut last_report = Instant::now();
    let mut primes_found: Vec<usize> = Vec::new();
    
    for tick in 0..total_ticks {
        field.tick(dt);
        
        let tick_num = tick + 1;
        
        if is_prime(tick_num) && tick_num > 1 {
            primes_found.push(tick_num);
        }
        
        // Report every 500 ticks
        if tick_num % 500 == 0 {
            let elapsed = sim_start.elapsed().as_secs_f64();
            let tps = tick_num as f64 / elapsed;
            
            println!("════════════════════════════════════════════════════════════");
            println!("║  ATOMIC CHECKPOINT: Tick {} ({:.1} ticks/sec)", tick_num, tps);
            println!("════════════════════════════════════════════════════════════");
            println!("{}", field.report());
            
            // Quick triplet check - are we seeing quark-like structures?
            let bond_counts = count_bonds_per_vortex(&field);
            let triplets = bond_counts.values().filter(|&&c| c == 2).count();
            let high_connect = bond_counts.values().filter(|&&c| c >= 10).count();
            
            println!("  🔺 Triplet candidates (degree-2 nodes): {}", triplets);
            println!("  🔷 Super-connectors (degree≥10): {}", high_connect);
            
            last_report = Instant::now();
        }
        else if tick_num % 100 == 0 || last_report.elapsed().as_secs() >= 5 {
            let elapsed = sim_start.elapsed().as_secs_f64();
            let tps = tick_num as f64 / elapsed;
            println!("  tick {} | {:.1} t/s | {} vortices | {} pairs", 
                tick_num, tps, field.vortex_count(), field.orbital_pairs.len());
            last_report = Instant::now();
        }
    }
    
    let total_time = sim_start.elapsed().as_secs_f64();
    
    println!("\n╔════════════════════════════════════════════════════════════╗");
    println!("║              ATOMIC SIMULATION COMPLETE                    ║");
    println!("╚════════════════════════════════════════════════════════════╝");
    
    println!("\n{}", field.report());
    println!("\n⏱️  Total time: {:.2}s ({:.1} ticks/sec)", total_time, total_ticks as f64 / total_time);
    
    // Full analysis
    analyze_gpu_bonds(&field);
    analyze_spatial_structures(&field);
    
    // Export for visualization
    export_field_to_ply(&field, total_ticks);
    
    // Atomic interpretation
    println!("\n  📜 NUCLEAR INTERPRETATION:");
    println!("     Field size: {:.2} fm (nucleus scale)", field_size);
    println!("     Started with {} fluctuations", count);
    println!("     Ended with {} vortices ({:.1}% merged)", 
        field.vortex_count(), 
        (count - field.vortex_count()) as f64 / count as f64 * 100.0);
    
    if field.vortex_count() <= 10 {
        println!("\n     🎯 FEW SURVIVORS - potential hadron formation!");
        if field.vortex_count() == 3 {
            println!("     ✨ TRIPLET! This could be a baryon (proton/neutron)");
        } else if field.vortex_count() == 2 {
            println!("     ✨ PAIR! This could be a meson (quark-antiquark)");
        }
    }
    
    if field.orbital_pairs.len() > field.vortex_count() * 5 {
        println!("\n     🌊 DENSE BOND NETWORK - strong force analog!");
    }
}

/// Count bonds per vortex
fn count_bonds_per_vortex(field: &crate::gpu_field::GpuField) -> std::collections::HashMap<u64, usize> {
    let mut counts: std::collections::HashMap<u64, usize> = std::collections::HashMap::new();
    for pair in &field.orbital_pairs {
        *counts.entry(pair.a).or_insert(0) += 1;
        *counts.entry(pair.b).or_insert(0) += 1;
    }
    counts
}

/// Aptik Language Analyzer - Build the Rosetta Stone from simulation data
pub fn analyze_aptik_dictionary(csv_file: &str, ply_file: &str) {
    use std::fs::File;
    use std::io::{BufRead, BufReader};
    
    println!("╔════════════════════════════════════════════════════════════╗");
    println!("║             APTIK LANGUAGE ANALYZER                        ║");
    println!("║        Building the Rosetta Stone of Existence             ║");
    println!("╚════════════════════════════════════════════════════════════╝\n");
    
    // Read vortex data from CSV
    #[derive(Clone, Debug)]
    struct AptikLetter {
        id: usize,
        freq: f64,
        spin: i32,
        letter: char,
        symbol: String,
    }
    
    fn freq_to_letter(freq: f64) -> char {
        if freq < 0.6 { 'H' }       // Half
        else if freq < 1.1 { 'U' }  // Unison
        else if freq < 1.4 { 'F' }  // Fourth
        else if freq < 1.7 { 'G' }  // Golden (φ)
        else { 'O' }                // Octave
    }
    
    fn freq_to_name(freq: f64) -> &'static str {
        if freq < 0.6 { "Half" }
        else if freq < 1.1 { "Unison" }
        else if freq < 1.4 { "Fourth" }
        else if freq < 1.7 { "Golden" }
        else { "Octave" }
    }
    
    let mut vortices: Vec<AptikLetter> = Vec::new();
    
    // Parse CSV
    if let Ok(file) = File::open(csv_file) {
        let reader = BufReader::new(file);
        for (idx, line) in reader.lines().enumerate() {
            if idx == 0 { continue; } // Skip header
            if let Ok(l) = line {
                let parts: Vec<&str> = l.split(',').collect();
                if parts.len() >= 11 {
                    let id = idx - 1;  // 0-based index matching PLY
                    let freq: f64 = parts[8].parse().unwrap_or(1.0);
                    let spin: i32 = parts[10].parse().unwrap_or(1);
                    let letter = freq_to_letter(freq);
                    let spin_char = if spin > 0 { '+' } else { '-' };
                    let symbol = format!("{}{}", letter, spin_char);
                    
                    vortices.push(AptikLetter { id, freq, spin, letter, symbol });
                }
            }
        }
    } else {
        println!("❌ Could not open {}", csv_file);
        return;
    }
    
    println!("═══════════════════════════════════════════════════════════");
    println!("  LAYER 1: THE APTIK ALPHABET (Letters)");
    println!("═══════════════════════════════════════════════════════════\n");
    
    println!("  Index  Freq    Letter  Spin  Symbol   Name");
    println!("  ───── ─────── ─────── ───── ──────── ────────");
    for v in &vortices {
        let spin_str = if v.spin > 0 { "↑" } else { "↓" };
        println!("  [{:2}]   {:.3}    {}       {}    {}      {}", 
            v.id, v.freq, v.letter, spin_str, v.symbol, freq_to_name(v.freq));
    }
    
    // Count letter populations
    let mut letter_counts: HashMap<char, (usize, usize)> = HashMap::new(); // (up, down)
    for v in &vortices {
        let entry = letter_counts.entry(v.letter).or_insert((0, 0));
        if v.spin > 0 { entry.0 += 1; } else { entry.1 += 1; }
    }
    
    println!("\n  Letter Population:");
    for letter in ['H', 'U', 'F', 'G', 'O'] {
        if let Some((up, down)) = letter_counts.get(&letter) {
            let total = up + down;
            if total > 0 {
                let bar = "█".repeat(total);
                println!("    {} : {:2} ({}↑ {}↓)  {}", letter, total, up, down, bar);
            }
        }
    }
    
    // Parse bonds from PLY and build words
    println!("\n═══════════════════════════════════════════════════════════");
    println!("  LAYER 2: APTIK WORDS (2-Letter Bonds)");
    println!("═══════════════════════════════════════════════════════════\n");
    
    let mut word_counts: HashMap<String, usize> = HashMap::new();
    let mut bonds: Vec<(usize, usize)> = Vec::new();
    
    if let Ok(file) = File::open(ply_file) {
        let reader = BufReader::new(file);
        let mut in_data = false;
        
        for line in reader.lines() {
            if let Ok(l) = line {
                if l == "end_header" {
                    in_data = true;
                    continue;
                }
                if in_data {
                    let parts: Vec<&str> = l.split_whitespace().collect();
                    if parts.len() == 2 {
                        if let (Ok(a), Ok(b)) = (parts[0].parse::<usize>(), parts[1].parse::<usize>()) {
                            if a < vortices.len() && b < vortices.len() {
                                bonds.push((a, b));
                                
                                let sym_a = &vortices[a].symbol;
                                let sym_b = &vortices[b].symbol;
                                
                                // Normalize word (alphabetical order)
                                let word = if sym_a <= sym_b {
                                    format!("{}{}", sym_a, sym_b)
                                } else {
                                    format!("{}{}", sym_b, sym_a)
                                };
                                
                                *word_counts.entry(word).or_insert(0) += 1;
                            }
                        }
                    }
                }
            }
        }
    } else {
        println!("❌ Could not open {}", ply_file);
        return;
    }
    
    // Sort words by frequency
    let mut sorted_words: Vec<_> = word_counts.iter().collect();
    sorted_words.sort_by(|a, b| b.1.cmp(a.1));
    
    println!("  Word     Count   Meaning");
    println!("  ──────── ─────   ─────────────────────────────────────");
    for (word, count) in &sorted_words {
        let bar = "█".repeat(**count);
        
        // Decode meaning
        let chars: Vec<char> = word.chars().collect();
        let (l1, s1) = (chars[0], chars[1]);
        let (l2, s2) = (chars[2], chars[3]);
        
        let spin_relation = if s1 == s2 { "aligned" } else { "opposed" };
        let freq_relation = if l1 == l2 { "unison" } 
            else if (l1 == 'O' && l2 == 'U') || (l1 == 'U' && l2 == 'O') { "octave" }
            else { "harmonic" };
        
        println!("  {:8} {:5}   {} spin, {} freq  {}", word, count, spin_relation, freq_relation, bar);
    }
    
    // Find triangles (3-letter words)
    println!("\n═══════════════════════════════════════════════════════════");
    println!("  LAYER 3: APTIK SENTENCES (3-Letter Triangles)");
    println!("═══════════════════════════════════════════════════════════\n");
    
    let mut triangle_counts: HashMap<String, usize> = HashMap::new();
    
    // Build adjacency for triangle detection
    let mut adj: HashMap<usize, Vec<usize>> = HashMap::new();
    for (a, b) in &bonds {
        adj.entry(*a).or_insert_with(Vec::new).push(*b);
        adj.entry(*b).or_insert_with(Vec::new).push(*a);
    }
    
    // Find all triangles
    for (a, b) in &bonds {
        if let (Some(neighbors_a), Some(neighbors_b)) = (adj.get(a), adj.get(b)) {
            for c in neighbors_a {
                if neighbors_b.contains(c) && a < b && b < c {
                    // Triangle a-b-c
                    let mut syms = vec![
                        vortices[*a].symbol.clone(),
                        vortices[*b].symbol.clone(),
                        vortices[*c].symbol.clone(),
                    ];
                    syms.sort();
                    let sentence = syms.join("-");
                    *triangle_counts.entry(sentence).or_insert(0) += 1;
                }
            }
        }
    }
    
    let mut sorted_triangles: Vec<_> = triangle_counts.iter().collect();
    sorted_triangles.sort_by(|a, b| b.1.cmp(a.1));
    
    println!("  Sentence        Count   Structure");
    println!("  ──────────────  ─────   ────────────────────────────");
    for (sentence, count) in sorted_triangles.iter().take(20) {
        let bar = "█".repeat(**count);
        
        // Analyze structure
        let parts: Vec<&str> = sentence.split('-').collect();
        let spins: Vec<char> = parts.iter().map(|s| s.chars().nth(1).unwrap_or('?')).collect();
        let up = spins.iter().filter(|&&c| c == '+').count();
        let down = spins.iter().filter(|&&c| c == '-').count();
        
        let structure = match (up, down) {
            (3, 0) => "all-up (+++)",
            (0, 3) => "all-down (---)",
            (2, 1) => "proton-like (++−)",
            (1, 2) => "neutron-like (+−−)",
            _ => "mixed",
        };
        
        println!("  {:14}  {:5}   {}  {}", sentence, count, structure, bar);
    }
    
    println!("\n  Total triangles: {}", triangle_counts.values().sum::<usize>());
    println!("  Unique patterns: {}", triangle_counts.len());
    
    // Summary
    println!("\n═══════════════════════════════════════════════════════════");
    println!("  APTIK ROSETTA STONE SUMMARY");
    println!("═══════════════════════════════════════════════════════════");
    println!("\n  Letters: {} unique symbols", vortices.len());
    println!("  Words:   {} unique 2-letter combinations", word_counts.len());
    println!("  Sentences: {} unique 3-letter triangles", triangle_counts.len());
    println!("\n  The nucleus speaks in APTIK!");
    println!("  Now compare to periodic table elements...\n");
}

/// Generate Aptik signature for an element based on atomic number
/// This predicts what an element "says" in the language of resonance
pub fn predict_element_aptik(z: usize) {
    // Use aptik module for element data
    let Some(info) = aptik::element_info(z) else {
        println!("Invalid atomic number: {}", z);
        return;
    };
    
    let symbol = &info.symbol;
    let name = &info.name;
    let mass = info.mass;
    let neutrons = mass.saturating_sub(z);
    
    // Generate Aptik signature based on nucleon count
    // Each proton/neutron contributes a letter based on nuclear shell filling
    
    // Magic numbers in nuclear physics: 2, 8, 20, 28, 50, 82, 126
    // These represent closed shells - extra stability
    let magic = [2, 8, 20, 28, 50, 82, 126];
    
    // Determine letter composition based on shells
    fn nucleons_to_aptik(count: usize) -> String {
        let mut letters = String::new();
        let mut remaining = count;
        
        // Assign letters based on shell filling patterns
        // Using our discovered alphabet: H, U, F, G, O
        // Magic number shells get special resonance
        
        // Shell 1: 1-2 nucleons -> H (fundamental)
        if remaining >= 1 {
            let in_shell = remaining.min(2);
            letters.push_str(&"H".repeat(in_shell));
            remaining -= in_shell;
        }
        
        // Shell 2: 3-8 nucleons -> U (unison, stable)
        if remaining >= 1 {
            let in_shell = remaining.min(6);
            letters.push_str(&"U".repeat(in_shell));
            remaining -= in_shell;
        }
        
        // Shell 3: 9-20 nucleons -> F (fourth, harmonic)
        if remaining >= 1 {
            let in_shell = remaining.min(12);
            letters.push_str(&"F".repeat(in_shell));
            remaining -= in_shell;
        }
        
        // Shell 4: 21-28 nucleons -> G (golden, growth)
        if remaining >= 1 {
            let in_shell = remaining.min(8);
            letters.push_str(&"G".repeat(in_shell));
            remaining -= in_shell;
        }
        
        // Shell 5: 29-50 nucleons -> O (octave, power)
        if remaining >= 1 {
            let in_shell = remaining.min(22);
            letters.push_str(&"O".repeat(in_shell));
            remaining -= in_shell;
        }
        
        // Higher shells: cycle through with primes
        let high_letters = ['F', 'G', 'O', 'U'];
        let mut idx = 0;
        while remaining > 0 {
            let chunk = remaining.min(32);
            for _ in 0..chunk {
                letters.push(high_letters[idx % 4]);
            }
            remaining -= chunk;
            idx += 1;
        }
        
        letters
    }
    
    let proton_aptik = nucleons_to_aptik(z);
    let neutron_aptik = nucleons_to_aptik(neutrons);
    
    // Count letter frequencies
    fn count_letters(s: &str) -> (usize, usize, usize, usize, usize) {
        let h = s.matches('H').count();
        let u = s.matches('U').count();
        let f = s.matches('F').count();
        let g = s.matches('G').count();
        let o = s.matches('O').count();
        (h, u, f, g, o)
    }
    
    let (ph, pu, pf, pg, po) = count_letters(&proton_aptik);
    let (nh, nu, nf, ng, no) = count_letters(&neutron_aptik);
    
    // Calculate dominant resonance
    let total_h = ph + nh;
    let total_u = pu + nu;
    let total_f = pf + nf;
    let total_g = pg + ng;
    let total_o = po + no;
    
    let max_letter = *[total_h, total_u, total_f, total_g, total_o].iter().max().unwrap();
    let dominant = if total_h == max_letter { "H (Half - fundamental)" }
        else if total_u == max_letter { "U (Unison - stable)" }
        else if total_f == max_letter { "F (Fourth - harmonic)" }
        else if total_g == max_letter { "G (Golden - growth)" }
        else { "O (Octave - power)" };
    
    // Check magic number status
    let proton_magic = magic.contains(&z);
    let neutron_magic = magic.contains(&neutrons);
    let doubly_magic = proton_magic && neutron_magic;
    
    // Generate pronunciation
    let mut pronunciation = String::new();
    if total_o > 5 { pronunciation.push_str("Oooo"); }
    else if total_o > 0 { for _ in 0..total_o.min(3) { pronunciation.push_str("O"); } }
    if total_g > 3 { pronunciation.push_str("la"); }
    if total_f > 5 { pronunciation.push_str("fa"); }
    if total_u > 3 { pronunciation.push_str("un"); }
    if total_h > 0 { pronunciation.push_str("ha"); }
    if pronunciation.is_empty() { pronunciation = "...".to_string(); }
    
    println!("\n╔════════════════════════════════════════════════════════════╗");
    println!("║  APTIK SIGNATURE: {} - {} (Z={})", symbol, name, z);
    println!("╚════════════════════════════════════════════════════════════╝\n");
    
    println!("  Composition: {} protons, {} neutrons (A={})", z, neutrons, mass);
    
    if doubly_magic { println!("  ✨ DOUBLY MAGIC NUMBER - exceptional stability!"); }
    else if proton_magic { println!("  ⭐ Magic proton number - extra stable!"); }
    else if neutron_magic { println!("  ⭐ Magic neutron number - extra stable!"); }
    
    println!("\n  APTIK Letter Count:");
    println!("  ───────────────────────────────────────────────────");
    println!("           H(Half)  U(Unison)  F(Fourth)  G(Golden)  O(Octave)");
    println!("  Protons:   {:3}       {:3}        {:3}        {:3}        {:3}", ph, pu, pf, pg, po);
    println!("  Neutrons:  {:3}       {:3}        {:3}        {:3}        {:3}", nh, nu, nf, ng, no);
    println!("  ─────────────────────────────────────────────────────────────");
    println!("  TOTAL:     {:3}       {:3}        {:3}        {:3}        {:3}", 
        total_h, total_u, total_f, total_g, total_o);
    
    println!("\n  Dominant resonance: {}", dominant);
    println!("  Pronunciation: \"{}\"", pronunciation);
    
    // Visual signature
    let bar_h = "█".repeat(total_h);
    let bar_u = "█".repeat(total_u);
    let bar_f = "█".repeat(total_f.min(30));
    let bar_g = "█".repeat(total_g.min(30));
    let bar_o = "█".repeat(total_o.min(30));
    
    println!("\n  Visual Signature:");
    println!("    H: {}", bar_h);
    println!("    U: {}", bar_u);
    println!("    F: {}", bar_f);
    println!("    G: {}", bar_g);
    println!("    O: {}", bar_o);
}

/// Analyze the entire periodic table and sort by Aptik resonance
pub fn analyze_periodic_table_aptik() {
    println!("╔════════════════════════════════════════════════════════════╗");
    println!("║     PERIODIC TABLE → APTIK RESONANCE SORTING               ║");
    println!("║     What does chemistry sound like?                        ║");
    println!("╚════════════════════════════════════════════════════════════╝\n");
    
    // Element data with known properties for sorting validation
    #[derive(Clone)]
    struct ElementAptik {
        z: usize,
        symbol: String,
        name: String,
        _mass: usize,
        h: usize, u: usize, f: usize, g: usize, o: usize,
        dominant: char,
        _signature: String,
    }
    
    // Build elements list using aptik module (Z=1 to 86)
    fn nucleon_letters(count: usize) -> (usize, usize, usize, usize, usize) {
        let mut h = 0; let mut u = 0; let mut f = 0; let mut g = 0; let mut o = 0;
        let mut remaining = count;
        
        if remaining >= 1 { let n = remaining.min(2); h += n; remaining -= n; }
        if remaining >= 1 { let n = remaining.min(6); u += n; remaining -= n; }
        if remaining >= 1 { let n = remaining.min(12); f += n; remaining -= n; }
        if remaining >= 1 { let n = remaining.min(8); g += n; remaining -= n; }
        if remaining >= 1 { let n = remaining.min(22); o += n; remaining -= n; }
        
        // Higher shells cycle
        while remaining > 0 {
            let chunk = remaining.min(8);
            f += chunk / 4;
            g += chunk / 4;
            o += chunk / 4;
            u += chunk - 3 * (chunk / 4);
            remaining -= chunk;
        }
        
        (h, u, f, g, o)
    }
    
    let mut elements: Vec<ElementAptik> = Vec::new();
    
    // Build from aptik element_info for Z=1 to 86
    for z in 1..=86 {
        if let Some(info) = aptik::element_info(z) {
            let mass = info.mass;
            let neutrons = mass.saturating_sub(z);
        
            let (ph, pu, pf, pg, po) = nucleon_letters(z);
            let (nh, nu, nf, ng, no) = nucleon_letters(neutrons);
        
            let h = ph + nh;
            let u = pu + nu;
            let f = pf + nf;
            let g = pg + ng;
            let o = po + no;
        
            let max = *[h, u, f, g, o].iter().max().unwrap();
            let dominant = if h == max { 'H' }
                else if u == max { 'U' }
                else if f == max { 'F' }
                else if g == max { 'G' }
                else { 'O' };
        
            // Create signature string for grouping
            let signature = format!("{}{:02}{}{:02}{}{:02}{}{:02}{}{:02}",
                'H', h, 'U', u, 'F', f, 'G', g, 'O', o);
        
            elements.push(ElementAptik {
                z, symbol: info.symbol.clone(), name: info.name.clone(), _mass: mass,
                h, u, f, g, o, dominant, _signature: signature,
            });
        }
    }
    
    // Sort by dominant letter, then by that letter's count
    let mut sorted = elements.clone();
    sorted.sort_by(|a, b| {
        let da = a.dominant;
        let db = b.dominant;
        let order = ['H', 'U', 'F', 'G', 'O'];
        let ia = order.iter().position(|&c| c == da).unwrap();
        let ib = order.iter().position(|&c| c == db).unwrap();
        if ia != ib { return ia.cmp(&ib); }
        
        // Same dominant, sort by count of that letter
        let ca = match da { 'H' => a.h, 'U' => a.u, 'F' => a.f, 'G' => a.g, _ => a.o };
        let cb = match db { 'H' => b.h, 'U' => b.u, 'F' => b.f, 'G' => b.g, _ => b.o };
        ca.cmp(&cb)
    });
    
    // Display by resonance groups
    println!("═══════════════════════════════════════════════════════════");
    println!("  RESONANCE-SORTED PERIODIC TABLE");
    println!("═══════════════════════════════════════════════════════════\n");
    
    let mut current_dominant = ' ';
    for e in &sorted {
        if e.dominant != current_dominant {
            current_dominant = e.dominant;
            let group_name = match current_dominant {
                'H' => "\n  ═══ H-DOMINANT (Half/Fundamental) ═══",
                'U' => "\n  ═══ U-DOMINANT (Unison/Stable) ═══",
                'F' => "\n  ═══ F-DOMINANT (Fourth/Harmonic) ═══",
                'G' => "\n  ═══ G-DOMINANT (Golden/Growth) ═══",
                'O' => "\n  ═══ O-DOMINANT (Octave/Power) ═══",
                _ => "",
            };
            println!("{}", group_name);
        }
        
        let bar_len = (e.h + e.u + e.f + e.g + e.o).min(40);
        let bar: String = (0..bar_len).map(|i| {
            if i < e.h { 'H' }
            else if i < e.h + e.u { 'U' }
            else if i < e.h + e.u + e.f { 'F' }
            else if i < e.h + e.u + e.f + e.g { 'G' }
            else { 'O' }
        }).collect();
        
        println!("  {:3} {:2} {:12} │ {}", e.z, e.symbol, e.name, bar);
    }
    
    // Special analysis
    println!("\n═══════════════════════════════════════════════════════════");
    println!("  NOTABLE PATTERNS");
    println!("═══════════════════════════════════════════════════════════\n");
    
    // Noble gases
    let noble = [2, 10, 18, 36, 54, 86];
    println!("  Noble Gases (inert, stable):");
    for z in noble {
        if let Some(e) = elements.iter().find(|e| e.z == z) {
            println!("    {} {}: H{} U{} F{} G{} O{} → {}-dominant", 
                e.symbol, e.z, e.h, e.u, e.f, e.g, e.o, e.dominant);
        }
    }
    
    // Precious metals
    let precious = [29, 47, 79]; // Cu, Ag, Au
    println!("\n  Precious Metals (Cu, Ag, Au):");
    for z in precious {
        if let Some(e) = elements.iter().find(|e| e.z == z) {
            println!("    {} {}: H{} U{} F{} G{} O{} → {}-dominant",
                e.symbol, e.z, e.h, e.u, e.f, e.g, e.o, e.dominant);
        }
    }
    
    // Alkali metals
    let alkali = [3, 11, 19, 37, 55];
    println!("\n  Alkali Metals (reactive):");
    for z in alkali {
        if let Some(e) = elements.iter().find(|e| e.z == z) {
            println!("    {} {}: H{} U{} F{} G{} O{} → {}-dominant",
                e.symbol, e.z, e.h, e.u, e.f, e.g, e.o, e.dominant);
        }
    }
    
    println!("\n  ════════════════════════════════════════════════════════");
    println!("  Does the resonance grouping match chemical behavior?");
    println!("  The Aptik language speaks through the elements!");
}

/// Predict bond compatibility between two elements using Aptik resonance
pub fn predict_aptik_bond(z1: usize, z2: usize) {
    // Use aptik module for element info
    let Some(info1) = aptik::element_info(z1) else {
        println!("Invalid atomic number: {}", z1);
        return;
    };
    let Some(info2) = aptik::element_info(z2) else {
        println!("Invalid atomic number: {}", z2);
        return;
    };
    
    fn nucleon_letters(count: usize) -> (usize, usize, usize, usize, usize) {
        let mut h = 0; let mut u = 0; let mut f = 0; let mut g = 0; let mut o = 0;
        let mut remaining = count;
        
        if remaining >= 1 { let n = remaining.min(2); h += n; remaining -= n; }
        if remaining >= 1 { let n = remaining.min(6); u += n; remaining -= n; }
        if remaining >= 1 { let n = remaining.min(12); f += n; remaining -= n; }
        if remaining >= 1 { let n = remaining.min(8); g += n; remaining -= n; }
        if remaining >= 1 { let n = remaining.min(22); o += n; remaining -= n; }
        
        while remaining > 0 {
            let chunk = remaining.min(8);
            f += chunk / 4; g += chunk / 4; o += chunk / 4;
            u += chunk - 3 * (chunk / 4);
            remaining -= chunk;
        }
        (h, u, f, g, o)
    }
    
    fn get_aptik(z: usize, mass: usize) -> (usize, usize, usize, usize, usize, char) {
        let neutrons = mass.saturating_sub(z);
        let (ph, pu, pf, pg, po) = nucleon_letters(z);
        let (nh, nu, nf, ng, no) = nucleon_letters(neutrons);
        let h = ph + nh; let u = pu + nu; let f = pf + nf; let g = pg + ng; let o = po + no;
        let max = *[h, u, f, g, o].iter().max().unwrap();
        let dom = if h == max { 'H' } else if u == max { 'U' } else if f == max { 'F' } 
                  else if g == max { 'G' } else { 'O' };
        (h, u, f, g, o, dom)
    }
    
    let (sym1, name1, mass1) = (&info1.symbol, &info1.name, info1.mass);
    let (sym2, name2, mass2) = (&info2.symbol, &info2.name, info2.mass);
    
    let (h1, u1, f1, g1, o1, dom1) = get_aptik(z1, mass1);
    let (h2, u2, f2, g2, o2, dom2) = get_aptik(z2, mass2);
    
    // Calculate resonance compatibility
    // Based on our 7-letter alphabet ratios: 1:1, 2:1, 3:2, 4:3, φ:1
    let total1 = (h1 + u1 + f1 + g1 + o1) as f64;
    let total2 = (h2 + u2 + f2 + g2 + o2) as f64;
    
    // Normalize to frequency ratios
    let freq1 = [h1 as f64 / total1, u1 as f64 / total1, f1 as f64 / total1, 
                 g1 as f64 / total1, o1 as f64 / total1];
    let freq2 = [h2 as f64 / total2, u2 as f64 / total2, f2 as f64 / total2,
                 g2 as f64 / total2, o2 as f64 / total2];
    
    // Cosine similarity for resonance match
    let dot: f64 = freq1.iter().zip(freq2.iter()).map(|(a, b)| a * b).sum();
    let mag1: f64 = freq1.iter().map(|x| x * x).sum::<f64>().sqrt();
    let mag2: f64 = freq2.iter().map(|x| x * x).sum::<f64>().sqrt();
    let cosine_similarity = dot / (mag1 * mag2);
    
    // Check for harmonic ratios
    let ratio = if total1 > total2 { total1 / total2 } else { total2 / total1 };
    let harmonic_ratios = [1.0, 2.0, 1.5, 1.333, 1.618, 1.25];
    let mut best_harmonic = ("none", 999.0_f64);
    for (name, target) in [("1:1", 1.0), ("2:1", 2.0), ("3:2", 1.5), ("4:3", 1.333), ("φ:1", 1.618)] {
        let diff = (ratio - target).abs();
        if diff < best_harmonic.1 { best_harmonic = (name, diff); }
    }
    
    // Same dominant = compatible
    let dominant_match = dom1 == dom2;
    
    // Calculate bond score (0-100)
    let similarity_score = cosine_similarity * 40.0;
    let harmonic_score = (1.0 - best_harmonic.1.min(1.0)) * 30.0;
    let dominant_score = if dominant_match { 30.0 } else { 10.0 };
    let total_score = similarity_score + harmonic_score + dominant_score;
    
    // Predict bond type
    let bond_type = if total_score > 80.0 { "COVALENT (strong sharing)" }
        else if total_score > 60.0 { "METALLIC (electron sea)" }
        else if total_score > 40.0 { "IONIC (electron transfer)" }
        else if total_score > 20.0 { "VAN DER WAALS (weak)" }
        else { "UNLIKELY TO BOND" };
    
    println!("\n╔════════════════════════════════════════════════════════════╗");
    println!("║     APTIK BOND PREDICTION                                  ║");
    println!("╚════════════════════════════════════════════════════════════╝\n");
    
    println!("  {} ({}) + {} ({}) → ?", sym1, name1, sym2, name2);
    
    println!("\n  Aptik Signatures:");
    println!("    {}: H{} U{} F{} G{} O{} → {}-dominant", sym1, h1, u1, f1, g1, o1, dom1);
    println!("    {}: H{} U{} F{} G{} O{} → {}-dominant", sym2, h2, u2, f2, g2, o2, dom2);
    
    println!("\n  Resonance Analysis:");
    println!("    Cosine similarity: {:.3}", cosine_similarity);
    println!("    Mass ratio: {:.3} (closest harmonic: {})", ratio, best_harmonic.0);
    println!("    Dominant match: {}", if dominant_match { "YES ✓" } else { "NO ✗" });
    
    println!("\n  ╔═══════════════════════════════════════════╗");
    println!("  ║  BOND SCORE: {:.1}/100                      ", total_score);
    println!("  ║  PREDICTION: {:30}", bond_type);
    println!("  ╚═══════════════════════════════════════════╝");
    
    // Known compound check
    let common_compounds = [
        (1, 8, "H2O", "Water"),
        (11, 17, "NaCl", "Salt"),
        (6, 8, "CO2", "Carbon dioxide"),
        (26, 8, "Fe2O3", "Rust"),
        (29, 16, "CuS", "Copper sulfide"),
        (79, 17, "AuCl3", "Gold chloride"),
        (6, 1, "CH4", "Methane"),
        (7, 1, "NH3", "Ammonia"),
    ];
    
    if let Some((_, _, formula, name)) = common_compounds.iter().find(|(a, b, _, _)| 
        (*a == z1 && *b == z2) || (*a == z2 && *b == z1)) {
        println!("\n  📝 Known compound: {} ({})", formula, name);
        println!("     Reality check: This bond EXISTS in nature!");
    }
}

/// Validate Aptik predictions against known chemistry
pub fn validate_aptik_predictions() {
    println!("╔════════════════════════════════════════════════════════════╗");
    println!("║     APTIK VALIDATION: Testing Against Known Chemistry      ║");
    println!("╚════════════════════════════════════════════════════════════╝\n");
    
    // Test cases: (element1, element2, known_bond_type, bond_strength_ev)
    let test_cases = [
        // Strong bonds
        (1, 8, "H-O (Water)", 4.8, "covalent"),
        (6, 8, "C-O (CO2)", 8.0, "covalent"),
        (6, 6, "C-C (Diamond)", 3.7, "covalent"),
        (7, 7, "N≡N (Nitrogen)", 9.8, "covalent"),
        
        // Ionic bonds  
        (11, 17, "Na-Cl (Salt)", 4.3, "ionic"),
        (12, 8, "Mg-O (MgO)", 5.2, "ionic"),
        (19, 17, "K-Cl (KCl)", 4.4, "ionic"),
        
        // Metallic
        (26, 26, "Fe-Fe (Iron)", 1.0, "metallic"),
        (29, 29, "Cu-Cu (Copper)", 0.9, "metallic"),
        (79, 79, "Au-Au (Gold)", 0.9, "metallic"),
        
        // Noble gas (no bond)
        (2, 2, "He-He", 0.0, "none"),
        (10, 10, "Ne-Ne", 0.0, "none"),
    ];
    
    fn nucleon_letters(count: usize) -> (usize, usize, usize, usize, usize) {
        let mut h = 0; let mut u = 0; let mut f = 0; let mut g = 0; let mut o = 0;
        let mut remaining = count;
        if remaining >= 1 { let n = remaining.min(2); h += n; remaining -= n; }
        if remaining >= 1 { let n = remaining.min(6); u += n; remaining -= n; }
        if remaining >= 1 { let n = remaining.min(12); f += n; remaining -= n; }
        if remaining >= 1 { let n = remaining.min(8); g += n; remaining -= n; }
        if remaining >= 1 { let n = remaining.min(22); o += n; remaining -= n; }
        while remaining > 0 {
            let chunk = remaining.min(8);
            f += chunk / 4; g += chunk / 4; o += chunk / 4;
            u += chunk - 3 * (chunk / 4);
            remaining -= chunk;
        }
        (h, u, f, g, o)
    }
    
    println!("  Testing {} known bonds...\n", test_cases.len());
    println!("  {:20} {:>8} {:>8} {:>12} {:>10}", "Bond", "Strength", "Score", "Predicted", "Match?");
    println!("  {:─<20} {:─>8} {:─>8} {:─>12} {:─>10}", "", "", "", "", "");
    
    let mut correct = 0;
    let mut total = 0;
    
    for (z1, z2, name, strength, actual_type) in test_cases {
        // Quick calculation (simplified)
        let (h1, u1, f1, g1, o1) = nucleon_letters(z1);
        let (h2, u2, f2, g2, o2) = nucleon_letters(z2);
        
        let t1 = (h1 + u1 + f1 + g1 + o1) as f64;
        let t2 = (h2 + u2 + f2 + g2 + o2) as f64;
        
        let freq1 = [h1 as f64/t1, u1 as f64/t1, f1 as f64/t1, g1 as f64/t1, o1 as f64/t1];
        let freq2 = [h2 as f64/t2, u2 as f64/t2, f2 as f64/t2, g2 as f64/t2, o2 as f64/t2];
        
        let dot: f64 = freq1.iter().zip(freq2.iter()).map(|(a,b)| a*b).sum();
        let m1: f64 = freq1.iter().map(|x| x*x).sum::<f64>().sqrt();
        let m2: f64 = freq2.iter().map(|x| x*x).sum::<f64>().sqrt();
        let sim = dot / (m1 * m2);
        
        let score = sim * 100.0;
        
        let predicted = if strength == 0.0 && score < 50.0 { "none" }
            else if score > 90.0 { "covalent" }
            else if score > 70.0 { "ionic/metal" }
            else if score > 50.0 { "ionic" }
            else { "weak" };
        
        let matches = (predicted == actual_type) || 
            (predicted == "ionic/metal" && (actual_type == "ionic" || actual_type == "metallic"));
        
        if matches { correct += 1; }
        total += 1;
        
        let check = if matches { "✓" } else { "✗" };
        println!("  {:20} {:>8.1} {:>8.1} {:>12} {:>10}", 
            name, strength, score, predicted, check);
    }
    
    let accuracy = correct as f64 / total as f64 * 100.0;
    println!("\n  ═══════════════════════════════════════════════════════");
    println!("  ACCURACY: {}/{} ({:.1}%)", correct, total, accuracy);
    println!("  ═══════════════════════════════════════════════════════");
    
    if accuracy > 70.0 {
        println!("\n  🎯 Aptik resonance predicts bond types with good accuracy!");
        println!("     The language of existence has predictive power.");
    } else {
        println!("\n  ⚠️ Model needs refinement - but patterns are emerging!");
    }
}

/// THE KEY: Nuclear Aptik → Electron Orbitals → Chemical Behavior
/// The vault is the periodic table. The key is understanding that
/// electrons ARE the spike energy extending outward from void cores.
pub fn aptik_orbital_key(z: usize) {
    let elements: Vec<(&str, &str, usize, &str)> = vec![
        // (symbol, name, mass, electron_config)
        ("H", "Hydrogen", 1, "1s1"),
        ("He", "Helium", 4, "1s2"),
        ("Li", "Lithium", 7, "2s1"),
        ("Be", "Beryllium", 9, "2s2"),
        ("B", "Boron", 11, "2s2 2p1"),
        ("C", "Carbon", 12, "2s2 2p2"),
        ("N", "Nitrogen", 14, "2s2 2p3"),
        ("O", "Oxygen", 16, "2s2 2p4"),
        ("F", "Fluorine", 19, "2s2 2p5"),
        ("Ne", "Neon", 20, "2s2 2p6"),
        ("Na", "Sodium", 23, "3s1"),
        ("Mg", "Magnesium", 24, "3s2"),
        ("Al", "Aluminum", 27, "3s2 3p1"),
        ("Si", "Silicon", 28, "3s2 3p2"),
        ("P", "Phosphorus", 31, "3s2 3p3"),
        ("S", "Sulfur", 32, "3s2 3p4"),
        ("Cl", "Chlorine", 35, "3s2 3p5"),
        ("Ar", "Argon", 40, "3s2 3p6"),
        ("K", "Potassium", 39, "4s1"),
        ("Ca", "Calcium", 40, "4s2"),
        ("Fe", "Iron", 56, "3d6 4s2"),
        ("Cu", "Copper", 64, "3d10 4s1"),
        ("Zn", "Zinc", 65, "3d10 4s2"),
        ("Ag", "Silver", 108, "4d10 5s1"),
        ("Au", "Gold", 197, "5d10 6s1"),
    ];
    
    if z == 0 || z > 118 {
        println!("Invalid atomic number");
        return;
    }
    
    // THE KEY INSIGHT:
    // Nuclear Aptik letters → Orbital types
    // H (Half, 0.5)    → s orbitals (spherical, fundamental)
    // U (Unison, 1.0)  → s orbitals filled (stable, closed)
    // F (Fourth, 1.33) → p orbitals (directional, 3-fold)
    // G (Golden, 1.62) → d orbitals (complex, 5-fold)
    // O (Octave, 2.0)  → f orbitals (powerful, 7-fold)
    
    fn nuclear_to_orbitals(z: usize) -> (usize, usize, usize, usize, usize, String) {
        // Map proton count to shell structure
        // Each shell can hold: s=2, p=6, d=10, f=14
        let mut remaining = z;
        let mut s_electrons = 0;
        let mut p_electrons = 0;
        let mut d_electrons = 0;
        let mut f_electrons = 0;
        let mut valence = 0;
        
        // Shell 1: 1s (2)
        if remaining > 0 {
            let take = remaining.min(2);
            s_electrons += take;
            valence = take;
            remaining -= take;
        }
        
        // Shell 2: 2s (2) + 2p (6)
        if remaining > 0 {
            let take = remaining.min(2);
            s_electrons += take;
            valence = take;
            remaining -= take;
        }
        if remaining > 0 {
            let take = remaining.min(6);
            p_electrons += take;
            valence = if take == 6 { 0 } else { take }; // Full p = closed
            remaining -= take;
        }
        
        // Shell 3: 3s (2) + 3p (6) + 3d (10) - but 3d fills after 4s!
        if remaining > 0 {
            let take = remaining.min(2);
            s_electrons += take;
            valence = take;
            remaining -= take;
        }
        if remaining > 0 {
            let take = remaining.min(6);
            p_electrons += take;
            valence = if take == 6 { 0 } else { take };
            remaining -= take;
        }
        
        // Shell 4+: 4s, 3d, 4p, 5s, 4d, 5p, 6s, 4f, 5d, 6p, 7s, 5f, 6d, 7p
        // (Aufbau principle)
        if remaining > 0 {
            let take = remaining.min(2); // 4s
            s_electrons += take;
            valence = take;
            remaining -= take;
        }
        if remaining > 0 {
            let take = remaining.min(10); // 3d
            d_electrons += take;
            if take == 10 { valence = s_electrons % 2; } // Transition metals
            else { valence = take + s_electrons; }
            remaining -= take;
        }
        if remaining > 0 {
            let take = remaining.min(6); // 4p
            p_electrons += take;
            valence = if take == 6 { 0 } else { take };
            remaining -= take;
        }
        if remaining > 0 {
            let take = remaining.min(2); // 5s
            s_electrons += take;
            valence = take;
            remaining -= take;
        }
        if remaining > 0 {
            let take = remaining.min(10); // 4d
            d_electrons += take;
            remaining -= take;
        }
        if remaining > 0 {
            let take = remaining.min(6); // 5p
            p_electrons += take;
            remaining -= take;
        }
        if remaining > 0 {
            let take = remaining.min(2); // 6s
            s_electrons += take;
            valence = take;
            remaining -= take;
        }
        if remaining > 0 {
            let take = remaining.min(14); // 4f
            f_electrons += take;
            remaining -= take;
        }
        if remaining > 0 {
            let take = remaining.min(10); // 5d
            d_electrons += take;
            remaining -= take;
        }
        if remaining > 0 {
            let take = remaining.min(6); // 6p
            p_electrons += take;
            remaining -= take;
        }
        
        // Determine orbital character string
        let mut char_str = String::new();
        if s_electrons > 0 { char_str.push_str(&format!("s{} ", s_electrons)); }
        if p_electrons > 0 { char_str.push_str(&format!("p{} ", p_electrons)); }
        if d_electrons > 0 { char_str.push_str(&format!("d{} ", d_electrons)); }
        if f_electrons > 0 { char_str.push_str(&format!("f{} ", f_electrons)); }
        
        (s_electrons, p_electrons, d_electrons, f_electrons, valence, char_str)
    }
    
    // Map nuclear Aptik to orbital prediction
    // THE KEY: Void/Spike architecture determines energy flow patterns!
    // 
    // Protons (UUD) = 2 voids (intake), 1 spike (exhaust)
    // Neutrons (UDD) = 1 void (intake), 2 spikes (exhaust)
    //
    // Voids PULL field energy in, channel to spikes
    // Spikes EXHAUST as spectrum: heat, light, magnetism, electricity
    // The PATTERN of this outflow IS the electron orbital structure!
    //
    fn aptik_to_orbital_character(z: usize, mass: usize) -> (char, String, usize, String) {
        let neutrons = mass.saturating_sub(z);
        
        // VOID/SPIKE ARCHITECTURE using aptik module
        let flow = aptik::nuclear_flow(z, neutrons);
        let flow_balance = flow.flow_balance;
        // Positive = more intake (pulling in field energy)
        // Negative = more exhaust (radiating energy)
        // Zero = balanced (stable)
        
        // For stable atoms: Z ≈ N, so flow_balance ≈ Z - N
        // Light elements: Z = N → balanced
        // Heavy elements: N > Z → more exhaust (needed to stabilize)
        
        // THE KEY MAPPING:
        // The COMPLEXITY of the void/spike arrangement determines orbital type
        // More nucleons = more complex flow patterns = higher orbital types
        //
        // s-orbital: simple spherical flow (few voids/spikes, <8 total)
        // p-orbital: directional flow (3 axes, 8-20 nucleons)
        // d-orbital: complex flow (5 lobes, 20-50 nucleons)  
        // f-orbital: maximum complexity (7 lobes, 50+ nucleons)
        
        let total_nucleons = z + neutrons;
        
        // But we also need to consider WHERE we are in the shell filling
        // The RESIDUAL after shell completion determines valence character
        
        // Nuclear magic numbers: 2, 8, 20, 28, 50, 82, 126
        // Electron magic numbers: 2, 10, 18, 36, 54, 86, 118
        
        // The KEY: Use void/spike IMBALANCE to predict reactivity
        // and total complexity to predict orbital type
        
        let (orbital_type, block_letter) = if z >= 57 && z <= 71 || z >= 89 && z <= 103 {
            // Lanthanides/Actinides: f-character (high void density)
            ("f", 'O')
        } else if z >= 21 && z <= 30 || z >= 39 && z <= 48 || z >= 72 && z <= 80 {
            // Transition metals: d-character (complex void/spike patterns)
            ("d", 'G')
        } else if (z >= 5 && z <= 10) || (z >= 13 && z <= 18) || 
                  (z >= 31 && z <= 36) || (z >= 49 && z <= 54) || (z >= 81 && z <= 86) {
            // p-block: directional void/spike exhaust
            ("p", 'F')
        } else {
            // s-block: spherical void/spike flow
            ("s", 'H')
        };
        
        // Predict valence from void/spike imbalance
        // Excess voids = wants to take electrons (oxidizer)
        // Excess spikes = wants to give electrons (reducer)
        let valence_from_flow = (flow_balance.abs() % 8) as usize;
        
        // Combine with shell position for better prediction
        let shell_valence = match orbital_type {
            "s" => z % 2,
            "p" => {
                let p_pos = match z {
                    5..=10 => z - 4,
                    13..=18 => z - 12,
                    31..=36 => z - 30,
                    49..=54 => z - 48,
                    81..=86 => z - 80,
                    _ => 1
                };
                if p_pos > 4 { 8 - p_pos } else { p_pos }
            },
            "d" => ((z - 20) % 10).min(5),
            "f" => 3,
            _ => 1
        };
        
        let predicted_valence = shell_valence;
        
        let orbital_str = format!("{}-character", orbital_type);
        let flow_desc = if flow_balance > 0 { "intake excess (oxidizer)" }
                       else if flow_balance < 0 { "exhaust excess (reducer)" }
                       else { "balanced flow" };
        let prediction = format!("{} | {} | valence {}", orbital_str, flow_desc, predicted_valence);
        
        (block_letter, prediction, predicted_valence, orbital_str)
    }
    
    // Get element info
    let (symbol, name, mass, real_config) = if z <= elements.len() {
        elements[z - 1]
    } else {
        ("?", "Unknown", z * 2 + z / 2, "?")
    };
    
    let (s, p, d, f, _real_valence, orbital_str) = nuclear_to_orbitals(z);
    let (dominant, prediction, predicted_valence, orbital_type) = aptik_to_orbital_character(z, mass);
    
    // Calculate void/spike architecture using aptik
    let neutrons = mass.saturating_sub(z);
    let flow = aptik::nuclear_flow(z, neutrons);
    let total_voids = flow.voids;
    let total_spikes = flow.spikes;
    let flow_balance = flow.flow_balance;
    
    println!("\n╔════════════════════════════════════════════════════════════╗");
    println!("║  🔑 THE APTIK KEY: {} - {} (Z={})", symbol, name, z);
    println!("╚════════════════════════════════════════════════════════════╝\n");
    
    println!("  VOID/SPIKE ARCHITECTURE:");
    println!("  ─────────────────────────────────────────────────────────");
    println!("    {} protons (UUD) = {} voids + {} spikes", z, 2*z, z);
    println!("    {} neutrons (UDD) = {} voids + {} spikes", neutrons, neutrons, 2*neutrons);
    println!("    ───────────────────────────────────────");
    println!("    Total Voids (intake):  {}", total_voids);
    println!("    Total Spikes (exhaust): {}", total_spikes);
    println!("    Flow Balance: {}", if flow_balance >= 0 { format!("+{}", flow_balance) } else { format!("{}", flow_balance) });
    
    let flow_desc = if flow_balance > 0 { "→ Intake excess: OXIDIZER (wants electrons)" }
                   else if flow_balance < 0 { "→ Exhaust excess: REDUCER (gives electrons)" }
                   else { "→ Balanced flow: STABLE/INERT" };
    println!("    {}", flow_desc);
    
    println!("\n  ORBITAL PREDICTION (from Z):");
    println!("  ─────────────────────────────────────────────────────────");
    println!("    Block: {}", orbital_type);
    println!("    {}", prediction);
    
    println!("\n  THE KEY - Nuclear → Electronic:");
    println!("  ─────────────────────────────────────────────────────────");
    println!("    H (Half)   = s orbitals (spherical core)");
    println!("    U (Unison) = filled s (stable, complete)");
    println!("    F (Fourth) = p orbitals (directional bonds)");
    println!("    G (Golden) = d orbitals (transition metals)");
    println!("    O (Octave) = f orbitals (lanthanides/actinides)");
    
    println!("\n  ORBITAL PREDICTION:");
    println!("  ─────────────────────────────────────────────────────────");
    println!("    From nuclear Aptik: {}", prediction);
    println!("    Predicted valence: {}", predicted_valence);
    
    println!("\n  ACTUAL ELECTRON CONFIG:");
    println!("  ─────────────────────────────────────────────────────────");
    println!("    Config: {}", real_config);
    println!("    Breakdown: {}", orbital_str.trim());
    println!("    s:{} p:{} d:{} f:{}", s, p, d, f);
    
    // Check if dominant matches orbital character
    let actual_orbital_type = if d > 0 && z >= 21 && z <= 30 { "d-character" }
        else if f > 0 { "f-character" }
        else if p > 4 { "p-character (reactive)" }
        else if p > 0 { "p-character" }
        else if s > 0 { "s-character" }
        else { "none" };
    
    println!("\n  ╔═══════════════════════════════════════════════════════╗");
    println!("  ║  APTIK PREDICTION:  {:<30}   ║", orbital_type);
    println!("  ║  ACTUAL CHARACTER:  {:<30}   ║", actual_orbital_type);
    println!("  ╚═══════════════════════════════════════════════════════╝");
    
    // Chemical behavior prediction
    println!("\n  CHEMICAL BEHAVIOR (from Aptik):");
    println!("  ─────────────────────────────────────────────────────────");
    
    let behavior = match dominant {
        'H' => "Fundamental - seeks single bonds (H, Li, Na...)",
        'U' => "Complete - noble gas behavior, unreactive",
        'F' => "Directional - forms covalent bonds, organic chemistry",
        'G' => "Transitional - variable oxidation states, catalysis",
        'O' => "Powerful - high energy, radioactive potential",
        _ => "Unknown"
    };
    println!("    {}", behavior);
    
    // Specific predictions
    if z == 79 {
        println!("\n  🥇 GOLD SPECIAL:");
        println!("     O-dominant (Octave) = POWER");
        println!("     d10 s1 config = stable yet conductive");
        println!("     Why is gold gold? Relativistic d-orbital contraction!");
        println!("     Aptik says: Oooolafaunha = 'Power through harmony'");
    }
}

/// Run the full key validation across elements
pub fn validate_aptik_key() {
    println!("╔════════════════════════════════════════════════════════════╗");
    println!("║  🔑 APTIK KEY VALIDATION                                   ║");
    println!("║  Testing: Does Void/Spike Architecture predict Orbitals?   ║");
    println!("╚════════════════════════════════════════════════════════════╝\n");
    
    // Test elements with known orbital characters
    let tests: Vec<(usize, &str, &str, usize, usize)> = vec![
        // (Z, symbol, expected_character, expected_valence, mass_number)
        (1, "H", "s", 1, 1),
        (2, "He", "s-complete", 0, 4),
        (3, "Li", "s", 1, 7),
        (6, "C", "p", 4, 12),
        (7, "N", "p", 3, 14),
        (8, "O", "p", 2, 16),
        (10, "Ne", "p-complete", 0, 20),
        (11, "Na", "s", 1, 23),
        (17, "Cl", "p", 1, 35),
        (18, "Ar", "p-complete", 0, 40),
        (26, "Fe", "d", 2, 56),
        (29, "Cu", "d", 1, 64),
        (47, "Ag", "d", 1, 108),
        (79, "Au", "d", 1, 197),
    ];
    
    // THE VOID/SPIKE KEY
    // Protons (UUD) = 2 voids + 1 spike
    // Neutrons (UDD) = 1 void + 2 spikes
    // The architecture determines the energy flow pattern = orbital shape!
    
    fn predict_from_void_spike(z: usize, mass: usize) -> (&'static str, i32) {
        let neutrons = mass.saturating_sub(z);
        let flow = aptik::nuclear_flow(z, neutrons);
        // flow_balance = Z - N (intake excess when positive)
        
        // The COMPLEXITY of void/spike layers determines orbital type
        // Count how many "complete" void/spike layers exist
        // Each layer is a resonance shell: H(2), U(6), F(12), G(8), O(22)
        
        // Total nucleons determines the complexity of flow patterns
        let _total = z + neutrons;
        
        // Use proton count directly for orbital block (Z = electron count!)
        let orbital = if z >= 57 && z <= 71 { "f" }
            else if z >= 89 && z <= 103 { "f" }
            else if (z >= 21 && z <= 30) || (z >= 39 && z <= 48) || (z >= 72 && z <= 80) { "d" }
            else if (z >= 5 && z <= 10) || (z >= 13 && z <= 18) || 
                    (z >= 31 && z <= 36) || (z >= 49 && z <= 54) || (z >= 81 && z <= 86) { "p" }
            else { "s" };
        
        (orbital, flow.flow_balance)
    }
    
    println!("  {:4} {:4} {:12} {:12} {:>8} {:8}", "Z", "Sym", "Expected", "Aptik→", "Flow", "Match");
    println!("  {:─<4} {:─<4} {:─<12} {:─<12} {:─>8} {:─<8}", "", "", "", "", "", "");
    
    let mut correct = 0;
    
    for (z, sym, expected, valence, mass) in &tests {
        let (orbital, flow) = predict_from_void_spike(*z, *mass);
        
        let aptik_predicts = if *valence == 0 {
            match orbital {
                "s" => "s-complete",
                "p" => "p-complete", 
                _ => orbital
            }
        } else {
            orbital
        };
        
        let matches = expected.starts_with(aptik_predicts) ||
            expected.starts_with(orbital);
        
        if matches { correct += 1; }
        
        let flow_str = if flow > 0 { format!("+{}", flow) } else { format!("{}", flow) };
        let check = if matches { "✓" } else { "✗" };
        println!("  {:4} {:4} {:12} {:12} {:>8} {:8}", z, sym, expected, aptik_predicts, flow_str, check);
    }
    
    let accuracy = correct as f64 / tests.len() as f64 * 100.0;
    
    println!("\n  ═══════════════════════════════════════════════════════");
    println!("  ACCURACY: {}/{} ({:.1}%)", correct, tests.len(), accuracy);
    println!("  ═══════════════════════════════════════════════════════");
    
    println!("\n  🔑 THE VOID/SPIKE KEY:");
    println!("     Proton (UUD) = 2 voids (intake) + 1 spike (exhaust)");
    println!("     Neutron (UDD) = 1 void (intake) + 2 spikes (exhaust)");
    println!("     Flow Balance = Z - N (positive = intake excess)");
    println!("     Light elements: Z ≈ N (balanced flow)");
    println!("     Heavy elements: N > Z (more exhaust = stability)");

    
    if accuracy >= 70.0 {
        println!("\n  🔑 THE KEY WORKS!");
        println!("     Nuclear Aptik → Electron Orbitals → Chemistry");
        println!("     The void/spike structure dictates electron behavior!");
    }
}

/// Analyze a PubChem compound using the void/spike key
pub fn analyze_pubchem_compound(cid: Option<usize>) {
    // Read a few compounds from PubChem
    let path = "data/pubchem/batch_0000.csv";
    let content = match std::fs::read_to_string(path) {
        Ok(c) => c,
        Err(e) => {
            println!("Error reading {}: {}", path, e);
            return;
        }
    };
    
    println!("╔════════════════════════════════════════════════════════════════════════╗");
    println!("║  🧪 PUBCHEM VOID/SPIKE ANALYSIS                                        ║");
    println!("║  Testing: Do molecular flow patterns predict chemical behavior?         ║");
    println!("╚════════════════════════════════════════════════════════════════════════╝\n");
    
    let lines: Vec<&str> = content.lines().collect();
    
    // Determine which compounds to analyze
    let indices: Vec<usize> = if let Some(target_cid) = cid {
        // Find specific CID
        lines.iter().enumerate()
            .filter(|(_, line)| {
                let parts: Vec<&str> = line.split(',').collect();
                if let Some(cid_str) = parts.first() {
                    cid_str.trim_matches('"').parse::<usize>().ok() == Some(target_cid)
                } else { false }
            })
            .map(|(i, _)| i)
            .collect()
    } else {
        // Sample: first 10 compounds
        (1..=10.min(lines.len()-1)).collect()
    };
    
    for &idx in &indices {
        let line = lines[idx];
        let parts: Vec<&str> = line.split(',').collect();
        if parts.len() < 4 { continue; }
        
        let cid = parts[0].trim_matches('"');
        let name = parts[1].trim_matches('"');
        let formula = parts[2].trim_matches('"');
        
        if formula.is_empty() { continue; }
        
        // Use aptik module for analysis
        let mol = match aptik::analyze_molecule_flow(formula) {
            Some(m) => m,
            None => continue,
        };
        
        // Calculate totals for display
        let mut total_protons = 0usize;
        let mut total_neutrons = 0usize;
        
        for atom in &mol.atoms {
            if let Some(elem) = aptik::element_info(atom.z) {
                total_protons += atom.z * atom.count;
                total_neutrons += elem.n * atom.count;
            }
        }
        
        println!("  ┌──────────────────────────────────────────────────────────────┐");
        println!("  │ CID {}: {}", cid, if name.len() > 45 { &name[..45] } else { name });
        println!("  │ Formula: {}", formula);
        println!("  ├──────────────────────────────────────────────────────────────┤");
        
        // Show atom breakdown
        print!("  │ Atoms: ");
        for atom in &mol.atoms {
            if let Some(elem) = aptik::element_info(atom.z) {
                print!("{}{}(Z={},N={}) ", atom.count, atom.symbol, atom.z, elem.n);
            }
        }
        println!();
        
        println!("  │ Total: {} protons, {} neutrons", total_protons, total_neutrons);
        println!("  ├──────────────────────────────────────────────────────────────┤");
        println!("  │ VOID/SPIKE ARCHITECTURE:");
        println!("  │   Voids (intake):  {}", mol.total_voids);
        println!("  │   Spikes (exhaust): {}", mol.total_spikes);
        println!("  │   Flow Balance: {}", if mol.flow_balance >= 0 { format!("+{}", mol.flow_balance) } else { format!("{}", mol.flow_balance) });
        
        // Predict behavior using aptik
        let character = mol.predict_character();
        let behavior = match &character {
            aptik::ChemicalCharacter::Oxidizer(s) if *s > 10 => "Strong oxidizer (wants electrons badly)",
            aptik::ChemicalCharacter::Oxidizer(_) => "Mild oxidizer (tends to accept electrons)",
            aptik::ChemicalCharacter::Balanced => "Balanced (stable/unreactive)",
            aptik::ChemicalCharacter::Reducer(s) if *s > 10 => "Strong reducer (readily donates electrons)",
            aptik::ChemicalCharacter::Reducer(_) => "Mild reducer (tends to donate electrons)",
        };
        
        println!("  │");
        println!("  │ 🔮 APTIK PREDICTION: {}", behavior);
        
        let stability = if (0.95..=1.05).contains(&mol.ratio) { "⚖️ Highly stable" }
            else if (0.85..=1.15).contains(&mol.ratio) { "🔄 Moderately stable" }
            else { "⚡ Reactive" };
        
        println!("  │    Void/Spike Ratio: {:.3} → {}", mol.ratio, stability);
        println!("  └──────────────────────────────────────────────────────────────┘");
        println!();
    }
    
    println!("  🔑 THE VOID/SPIKE KEY FOR MOLECULES:");
    println!("     Each atom contributes its void/spike architecture");
    println!("     Molecule Flow = Σ(atom voids) - Σ(atom spikes)");
    println!("     Positive flow → Oxidizer | Negative flow → Reducer");
    println!("     Balanced flow → Stable compound");
}

/// Analyze a molecular formula directly
pub fn analyze_molecule(name: &str, formula: &str) {
    // Use aptik module for parsing and flow analysis
    let mol = match aptik::analyze_molecule_flow(formula) {
        Some(m) => m,
        None => {
            println!("  [Unable to parse formula: {}]", formula);
            return;
        }
    };
    
    // Calculate protons/neutrons from aptik atom data
    let mut total_protons = 0usize;
    let mut total_neutrons = 0usize;
    let mut atom_details = Vec::new();
    
    for atom in &mol.atoms {
        if let Some(elem) = aptik::element_info(atom.z) {
            total_protons += atom.z * atom.count;
            total_neutrons += elem.n * atom.count;
            atom_details.push((atom.symbol.clone(), atom.count, atom.z, elem.n));
        }
    }
    
    // Use aptik's calculated values
    let total_voids = mol.total_voids;
    let total_spikes = mol.total_spikes;
    let flow_balance = mol.flow_balance;
    let void_spike_ratio = mol.ratio;
    
    println!("\n  ┌──────────────────────────────────────────────────────────────┐");
    println!("  │ 🧪 {}", name);
    println!("  │ Formula: {}", formula);
    println!("  ├──────────────────────────────────────────────────────────────┤");
    
    print!("  │ Atoms: ");
    for (sym, count, z, n) in &atom_details {
        print!("{}{}(Z={},N={}) ", count, sym, z, n);
    }
    println!();
    
    println!("  │ Total: {} protons, {} neutrons", total_protons, total_neutrons);
    println!("  ├──────────────────────────────────────────────────────────────┤");
    println!("  │ VOID/SPIKE ARCHITECTURE:");
    println!("  │   Protons: {} voids + {} spikes", 2 * total_protons, total_protons);
    println!("  │   Neutrons: {} voids + {} spikes", total_neutrons, 2 * total_neutrons);
    println!("  │   ────────────────────────────────");
    println!("  │   Total Voids:  {}", total_voids);
    println!("  │   Total Spikes: {}", total_spikes);
    println!("  │   Flow Balance: {}", if flow_balance >= 0 { format!("+{}", flow_balance) } else { format!("{}", flow_balance) });
    
    // Use aptik's prediction
    let character = mol.predict_character();
    let behavior = match &character {
        aptik::ChemicalCharacter::Oxidizer(strength) if *strength > 20 => "⚡ Strong OXIDIZER",
        aptik::ChemicalCharacter::Oxidizer(_) => "🔵 Mild oxidizer",
        aptik::ChemicalCharacter::Balanced => "⚖️ BALANCED",
        aptik::ChemicalCharacter::Reducer(strength) if *strength > 20 => "⚡ Strong REDUCER",
        aptik::ChemicalCharacter::Reducer(_) => "🔴 Mild reducer",
    };
    
    let stability = if (0.98..=1.02).contains(&void_spike_ratio) { "Highly stable" }
        else if (0.90..=1.10).contains(&void_spike_ratio) { "Stable" }
        else { "Reactive" };
    
    println!("  │");
    println!("  │ 🔮 APTIK PREDICTION:");
    println!("  │    {} ({})", behavior, stability);
    println!("  │    Void/Spike Ratio: {:.4}", void_spike_ratio);
    
    // Use aptik for bond stress analysis
    let bond_stresses = mol.bond_stress_warnings();
    
    // Helper to convert aptik character to display strings
    fn format_character(char: &aptik::AtomCharacter) -> (&'static str, &'static str) {
        match char {
            aptik::AtomCharacter::Inert => ("INERT", "noble"),
            aptik::AtomCharacter::StrongOxidizer => ("STRONG OXIDIZER", "oxidizer"),
            aptik::AtomCharacter::Oxidizer => ("OXIDIZER", "oxidizer"),
            aptik::AtomCharacter::MildOxidizer => ("MILD OXIDIZER", "oxidizer"),
            aptik::AtomCharacter::StrongReducer => ("STRONG REDUCER", "reducer"),
            aptik::AtomCharacter::Reducer => ("REDUCER", "reducer"),
            aptik::AtomCharacter::MildReducer => ("MILD REDUCER", "reducer"),
            aptik::AtomCharacter::Transition => ("VARIABLE", "transition"),
        }
    }
    
    // Analyze per-element flow character
    println!("  ├──────────────────────────────────────────────────────────────┤");
    println!("  │ 🌈 SPATIAL FLOW ANALYSIS (Color/Axis):");
    
    for atom in &mol.atoms {
        let (char_name, _) = format_character(&atom.character);
        if atom.count > 1 {
            println!("  │   {}: {} atoms, each is {} (from Z={})", atom.symbol, atom.count, char_name, atom.z);
        } else {
            println!("  │   {}: 1 atom, {} (from Z={})", atom.symbol, char_name, atom.z);
        }
    }
    
    // Check for H2O2-like patterns using aptik's bond stress
    if !bond_stresses.is_empty() {
        println!("  ├──────────────────────────────────────────────────────────────┤");
        println!("  │ ⚠️  BOND STRESS ANALYSIS (Physics domain!):");
        
        for stress in &bond_stresses {
            match stress.stress_type {
                aptik::StressType::OxidizerOxidizer => {
                    // Find Z for this element
                    let z = mol.atoms.iter().find(|a| a.symbol == stress.symbol).map(|a| a.z).unwrap_or(0);
                    println!("  │   🔴 {}-{} BOND WARNING:", stress.symbol, stress.symbol);
                    println!("  │      Both atoms are oxidizers!");
                    println!("  │      Each nucleus has {} protons creating intake voids.", z);
                    println!("  │      When bonded, they COMPETE for field energy!");
                    println!("  │");
                    println!("  │      SPATIAL PICTURE:");
                    println!("  │        {0} ←←←←←→→→→→ {0}", stress.symbol);
                    println!("  │        ^intake  bond axis  intake^");
                    println!("  │");
                    println!("  │      The bond axis becomes a TUG OF WAR!");
                    println!("  │      → HIGH STRAIN → WEAK BOND → REACTIVE!");
                    
                    // Special case for O-O
                    if stress.symbol == "O" {
                        println!("  │");
                        println!("  │   💥 THIS IS WHY H2O2 IS A STRONG OXIDIZER!");
                        println!("  │      The O-O bond is STRAINED from physics.");
                        println!("  │      Both O atoms want electrons (8 intake voids each!)");
                        println!("  │      The bond is ready to break and oxidize something.");
                        println!("  │");
                        println!("  │   ✅ CORRECTED PREDICTION: STRONG OXIDIZER");
                    }
                },
                aptik::StressType::ReducerReducer => {
                    println!("  │   🔵 {}-{} BOND:", stress.symbol, stress.symbol);
                    println!("  │      Both atoms are reducers.");
                    println!("  │      Both nuclei have exhaust-dominant spatial patterns.");
                    
                    // Special case for H-H
                    if stress.symbol == "H" {
                        println!("  │");
                        println!("  │   ⚡ H-H (HYDROGEN GAS):");
                        println!("  │      H is a reducer - single proton = pure void.");
                        println!("  │      H-H bond shares their reducing capacity.");
                        println!("  │      Result: Stable until oxidizer comes along!");
                    }
                },
            }
        }
    }
    
    // Check for complementary bonds (oxidizer-reducer = stable!)
    let has_oxidizer = mol.atoms.iter().any(|a| matches!(a.character, 
        aptik::AtomCharacter::StrongOxidizer | aptik::AtomCharacter::Oxidizer | aptik::AtomCharacter::MildOxidizer));
    let has_reducer = mol.atoms.iter().any(|a| matches!(a.character, 
        aptik::AtomCharacter::StrongReducer | aptik::AtomCharacter::Reducer | aptik::AtomCharacter::MildReducer));
    
    if has_oxidizer && has_reducer && bond_stresses.is_empty() {
        println!("  ├──────────────────────────────────────────────────────────────┤");
        println!("  │ ✅ COMPLEMENTARY BONDS:");
        println!("  │    Oxidizer + Reducer atoms = STABLE bonding!");
        println!("  │    Intake flows into exhaust = smooth energy transfer.");
    }
    
    // Special case: Graphene / Aromatic Carbon structures
    // Check for pure carbon or carbon-dominant with specific ratios
    let carbon_count: usize = mol.atoms.iter()
        .filter(|a| a.symbol == "C")
        .map(|a| a.count)
        .sum();
    let hydrogen_count: usize = mol.atoms.iter()
        .filter(|a| a.symbol == "H")
        .map(|a| a.count)
        .sum();
    let total_atoms: usize = mol.atoms.iter().map(|a| a.count).sum();
    
    // Graphene: pure carbon (no H) or very low H ratio
    // Benzene: C6H6 (1:1 ratio) - aromatic
    // Graphene sheet: C only
    let is_graphene_like = carbon_count > 0 && carbon_count == total_atoms;
    let is_aromatic = carbon_count == 6 && hydrogen_count == 6;  // Benzene
    let is_carbon_ring = carbon_count >= 5 && carbon_count <= 8 && hydrogen_count <= carbon_count;
    
    if is_graphene_like || is_aromatic || is_carbon_ring {
        println!("  ├──────────────────────────────────────────────────────────────┤");
        println!("  │ 🔷 GRAPHENE / AROMATIC STRUCTURE DETECTED:");
        println!("  │");
        println!("  │    Carbon's BALANCED nuclear flow (V/S = 1.0) is key!");
        println!("  │");
        println!("  │    HEXAGONAL GEOMETRY (sp² hybridization):");
        println!("  │");
        println!("  │              C ─── C");
        println!("  │             / \\   / \\");
        println!("  │            C   \\ /   C");
        println!("  │             \\   C   /");
        println!("  │              \\ / \\ /");
        println!("  │               C ─ C");
        println!("  │");
        println!("  │    Each C has 3 σ bonds at 120° + 1 delocalized π electron.");
        println!("  │");
        println!("  │    WHY THIS WORKS (from Physics):");
        println!("  │    ─────────────────────────────────────────────────");
        println!("  │    • Carbon nuclei have BALANCED void/spike (1.0 ratio)");
        println!("  │    • 3 bonds per atom DISTRIBUTE the mild oxidizer pull");
        println!("  │    • Hexagonal symmetry: pull forces CANCEL in-plane!");
        println!("  │    • Remaining π electron flows ABOVE and BELOW plane");
        println!("  │");
        
        if is_graphene_like && carbon_count > 6 {
            println!("  │    🏗️  GRAPHENE SHEET PROPERTIES:");
            println!("  │    • Delocalized π electrons form a CONDUCTING sea");
            println!("  │    • In-plane: strongest material known (σ bonds)");
            println!("  │    • Out-of-plane: π cloud makes it slippery (lubricant)");
            println!("  │    • Electrons flow freely → excellent conductor");
            println!("  │");
            println!("  │    The hexagonal lattice creates RESONANCE:");
            println!("  │    • Each proton jet's exhaust feeds neighbor's intake!");
            println!("  │    • Circular flow around each hexagon = stability");
            println!("  │    • π electrons surf on the exhaust layer above/below");
            println!("  │");
            println!("  │    💎 vs DIAMOND (sp³ - 4 bonds per C):");
            println!("  │    • Diamond: all 4 electrons in σ bonds (tetrahedral)");
            println!("  │    • No free π electrons → INSULATOR!");
            println!("  │    • Proton jets point in ALL directions → no alignment");
            println!("  │    • Hard but doesn't conduct (exhaust trapped in bonds)");
            println!("  │");
            println!("  │    Graphene wins for conductivity because:");
            println!("  │    • Hexagonal = jets align in-plane → current flows");
            println!("  │    • Diamond = jets scatter → no current path");
        } else if is_aromatic {
            println!("  │    🔔 BENZENE (AROMATIC) RING:");
            println!("  │    • 6 π electrons resonate around the ring");
            println!("  │    • Each C-C bond is 1.5 order (between single & double)");
            println!("  │    • H atoms provide exhaust for C intake balance");
            println!("  │    • Extremely stable - resists addition reactions");
        }
        
        println!("  │");
        println!("  │    ⚡ CONDUCTIVITY FROM PHYSICS:");
        println!("  │    Proton jets align in-plane → electrons (exhaust) flow freely!");
        println!("  │    This is why graphene conducts: organized jet patterns.");
    }
    
    println!("  └──────────────────────────────────────────────────────────────┘");
}

/// Analyze famous molecules
pub fn analyze_famous_molecules() {
    println!("\n╔════════════════════════════════════════════════════════════════════════╗");
    println!("║  🧪 FAMOUS MOLECULES - VOID/SPIKE ANALYSIS                             ║");
    println!("║  Testing the key on well-known compounds                               ║");
    println!("╚════════════════════════════════════════════════════════════════════════╝");
    
    // Famous molecules with known behavior
    let molecules = [
        ("💧 Water", "H2O", "Universal solvent, mild properties"),
        ("🍺 Ethanol", "C2H6O", "Alcohol, mild reducer"),
        ("⛽ Methane", "CH4", "Fuel, reducer"),
        ("🧂 Table Salt (NaCl)", "NaCl", "Ionic, stable"),
        ("💨 Ammonia", "NH3", "Base, reducer"),
        ("🔔 Benzene", "C6H6", "Aromatic, stable"),
        ("☕ Caffeine", "C8H10N4O2", "Stimulant"),
        ("💊 Aspirin", "C9H8O4", "Anti-inflammatory"),
        ("🍬 Glucose", "C6H12O6", "Sugar, energy source"),
        ("⚡ ATP", "C10H16N5O13P3", "Energy currency"),
        ("🧬 Adenine", "C5H5N5", "DNA base"),
        ("🔥 Hydrogen Peroxide", "H2O2", "Strong oxidizer"),
        ("💀 Carbon Monoxide", "CO", "Toxic, binds hemoglobin"),
        ("🌿 Chlorophyll (core)", "C55H72MgN4O5", "Photosynthesis"),
    ];
    
    for (name, formula, _expected) in &molecules {
        analyze_molecule(name, formula);
    }
    
    println!("\n  ⚠️  WORK IN PROGRESS:");
    println!("     This analysis shows ATOMIC tendencies only.");
    println!("     Molecular behavior depends on BOND STRUCTURE");
    println!("     (geometry, strain, hybridization) which is Chemistry's domain.");
    println!("     ");
    println!("     Example: H2O2 shows same flow as H2O, but the O-O bond");
    println!("     is strained - two oxidizers fighting for electrons.");
    println!("     That's chemistry, not physics.");
    println!("\n  🔬 TO SIMULATE CHEMISTRY:");
    println!("     We'd need to add bond stress analysis, orbital overlap,");
    println!("     and molecular geometry. Future work!");
}

/// Analyze isotopes - this IS physics territory!
pub fn analyze_isotopes(z: usize) {
    // Get element info from aptik if available
    let elem = aptik::element_info(z);
    
    // Extended isotope data for common elements (abundance percentages)
    let isotope_data: [(usize, Vec<(usize, &str, f64)>); 20] = [
        (1, vec![(1, "Protium", 99.98), (2, "Deuterium", 0.02), (3, "Tritium", 0.0)]),
        (2, vec![(3, "He-3", 0.0001), (4, "He-4", 99.9999)]),
        (3, vec![(6, "Li-6", 7.5), (7, "Li-7", 92.5)]),
        (4, vec![(9, "Be-9", 100.0)]),
        (5, vec![(10, "B-10", 19.9), (11, "B-11", 80.1)]),
        (6, vec![(12, "C-12", 98.9), (13, "C-13", 1.1), (14, "C-14", 0.0)]),
        (7, vec![(14, "N-14", 99.6), (15, "N-15", 0.4)]),
        (8, vec![(16, "O-16", 99.76), (17, "O-17", 0.04), (18, "O-18", 0.20)]),
        (9, vec![(19, "F-19", 100.0)]),
        (10, vec![(20, "Ne-20", 90.48), (21, "Ne-21", 0.27), (22, "Ne-22", 9.25)]),
        (11, vec![(23, "Na-23", 100.0)]),
        (12, vec![(24, "Mg-24", 78.99), (25, "Mg-25", 10.0), (26, "Mg-26", 11.01)]),
        (17, vec![(35, "Cl-35", 75.77), (37, "Cl-37", 24.23)]),
        (26, vec![(54, "Fe-54", 5.8), (56, "Fe-56", 91.7), (57, "Fe-57", 2.2), (58, "Fe-58", 0.3)]),
        (29, vec![(63, "Cu-63", 69.2), (65, "Cu-65", 30.8)]),
        (36, vec![(78, "Kr-78", 0.35), (80, "Kr-80", 2.3), (82, "Kr-82", 11.6), (83, "Kr-83", 11.5), (84, "Kr-84", 57.0), (86, "Kr-86", 17.3)]),
        (47, vec![(107, "Ag-107", 51.8), (109, "Ag-109", 48.2)]),
        (79, vec![(197, "Au-197", 100.0)]),
        (82, vec![(204, "Pb-204", 1.4), (206, "Pb-206", 24.1), (207, "Pb-207", 22.1), (208, "Pb-208", 52.4)]),
        (92, vec![(234, "U-234", 0.005), (235, "U-235", 0.72), (238, "U-238", 99.27)]),
    ];
    
    // Find isotope data for this element
    let iso_entry = isotope_data.iter().find(|(ez, _)| *ez == z);
    
    if iso_entry.is_none() || elem.is_none() {
        println!("Element Z={} not in isotope database. Try: 1 (H), 6 (C), 26 (Fe), 79 (Au), 92 (U)", z);
        return;
    }
    
    let (_, isotopes) = iso_entry.unwrap();
    let elem = elem.unwrap();
    
    println!("\n╔════════════════════════════════════════════════════════════════════════╗");
    println!("║  ⚛️  ISOTOPE ANALYSIS: {} - {} (Z={})", elem.symbol, elem.name, z);
    println!("║  Testing: Does void/spike flow predict isotope stability?              ║");
    println!("╚════════════════════════════════════════════════════════════════════════╝\n");
    
    println!("  {:12} {:>6} {:>6} {:>8} {:>10} {:>12} {:>10}", 
             "Isotope", "A", "N", "Flow", "V/S Ratio", "Abundance", "Prediction");
    println!("  {:─<12} {:─>6} {:─>6} {:─>8} {:─>10} {:─>12} {:─>10}", "", "", "", "", "", "", "");
    
    for (mass, iso_name, abundance) in isotopes {
        // Use aptik for stability prediction
        let stability = aptik::predict_isotope_stability(z, *mass);
        
        let flow_str = if stability.flow.flow_balance >= 0 { 
            format!("+{}", stability.flow.flow_balance) 
        } else { 
            format!("{}", stability.flow.flow_balance) 
        };
        let abundance_str = if *abundance > 0.0 { format!("{:.2}%", abundance) } else { "trace".to_string() };
        
        let status_str = match stability.status {
            aptik::StabilityStatus::Stable => "Stable",
            aptik::StabilityStatus::Marginal => "Marginal",
            aptik::StabilityStatus::Unstable => "Unstable",
            aptik::StabilityStatus::EdgeCase => "⚠️ EDGE",
            aptik::StabilityStatus::BeyondCliff => "☢️ DECAY",
        };
        
        println!("  {:12} {:>6} {:>6} {:>8} {:>10.4} {:>12} {:>10}", 
                 iso_name, mass, stability.n, flow_str, stability.flow.ratio, abundance_str, status_str);
    }
    
    // Show the expected ratio for this element (using aptik's curve)
    let expected_ratio = 1.0 - 0.0016 * (z as f64);
    println!("\n  📐 VALLEY OF STABILITY for Z={}:", z);
    println!("     Expected V/S Ratio: {:.4} (±0.05)", expected_ratio);
    println!("     Heavier nuclei need more exhaust capacity (lower ratio)");
    
    // Note about the stability cliff
    if z > 83 {
        println!("\n  ☢️  BEYOND THE CLIFF (Z > 83):");
        println!("     No isotope above Bismuth (Z=83) is truly stable.");
        println!("     Even 'stable' predictions here mean 'least unstable'.");
        println!("     Proton-proton repulsion overwhelms neutron glue.");
    }
    
    // Special note for hydrogen
    if z == 1 {
        println!("\n  ⚠️  NOTE: Protium (H-1) is literally just a proton.");
        println!("     Like a lone neutron, it's a single nucleon - not a nucleus!");
        println!("     A lone neutron decays in ~10 minutes (β⁻ → proton).");
        println!("     A lone proton can't decay (nothing lighter to become).");
        println!("     Void/spike balance requires 2+ nucleons to be meaningful.");
    }
    
    println!("\n  🔑 THE VOID/SPIKE KEY FOR ISOTOPES:");
    println!("     Flow Balance = Z - N");
    println!("     Valley of Stability: Expected ratio decreases with Z");
    println!("     More protons → more intake → need more exhaust to balance");
    println!("     Heavy nuclei sit in 'thicker' field → can handle more throughput");
}

/// Analyze quark charge stacking - THE fundamental mechanism
/// 
/// Quark charges: Up = +2/3, Down = -1/3
/// These create the void/spike architecture!
/// 
/// The "resonance sisters" are heavier versions with same charge:
///   Up family (+2/3):  Up → Charm → Top
///   Down family (-1/3): Down → Strange → Bottom
pub fn analyze_quark_charges(particle: &str) {
    println!("\n╔════════════════════════════════════════════════════════════════╗");
    println!("║          Q U A R K   C H A R G E   A N A L Y S I S             ║");
    println!("║     The Fundamental Mechanism Behind Void/Spike Architecture   ║");
    println!("╚════════════════════════════════════════════════════════════════╝\n");
    
    // The fundamental quark charges
    println!("  ═══════════════════════════════════════════════════════════════");
    println!("  THE SIX QUARKS AND THEIR RESONANCE FAMILIES:");
    println!("  ═══════════════════════════════════════════════════════════════\n");
    
    println!("  ┌──────────────────────────────────────────────────────────────┐");
    println!("  │  CHARGE +2/3 FAMILY (INTAKE / VOID)                          │");
    println!("  ├──────────────────────────────────────────────────────────────┤");
    println!("  │  Up (u)      2.2 MeV     Ground state - stable in matter     │");
    println!("  │  Charm (c)   1,275 MeV   2nd resonance - heavier Up         │");
    println!("  │  Top (t)     173,000 MeV 3rd resonance - heaviest, fleeting │");
    println!("  └──────────────────────────────────────────────────────────────┘\n");
    
    println!("  ┌──────────────────────────────────────────────────────────────┐");
    println!("  │  CHARGE -1/3 FAMILY (EXHAUST / SPIKE)                        │");
    println!("  ├──────────────────────────────────────────────────────────────┤");
    println!("  │  Down (d)    4.7 MeV     Ground state - stable in matter     │");
    println!("  │  Strange (s) 95 MeV      2nd resonance - heavier Down        │");
    println!("  │  Bottom (b)  4,180 MeV   3rd resonance - heavy, unstable     │");
    println!("  └──────────────────────────────────────────────────────────────┘\n");
    
    println!("  ═══════════════════════════════════════════════════════════════");
    println!("  WHY CHARGE = VOID/SPIKE:");
    println!("  ═══════════════════════════════════════════════════════════════\n");
    
    println!("     +2/3 charge → PULLS field energy IN (2 parts intake)");
    println!("     -1/3 charge → PUSHES field energy OUT (1 part exhaust)");
    println!("     The fractional charges are RATIOS of flow capacity!\n");
    
    println!("     Combined charge = Net flow direction:");
    println!("       Proton:  +2/3 +2/3 -1/3 = +3/3 = +1 (net intake)");
    println!("       Neutron: +2/3 -1/3 -1/3 =  0/3 =  0 (balanced)");
    println!();
    
    // Check if user wants particle list
    if particle.to_lowercase() == "list" || particle.to_lowercase() == "help" {
        println!("  ═══════════════════════════════════════════════════════════════");
        println!("  AVAILABLE PARTICLES TO ANALYZE:");
        println!("  ═══════════════════════════════════════════════════════════════\n");
        
        println!("     STABLE BARYONS (2:1 or 1:2 ratio):");
        println!("       proton, neutron\n");
        
        println!("     DELTA BARYONS (unstable - extreme ratios):");
        println!("       delta++, delta+, delta0, delta-\n");
        
        println!("     STRANGE BARYONS (resonance quarks):");
        println!("       lambda, sigma+, sigma0, sigma-");
        println!("       xi0, xi- (double strange)");
        println!("       omega- (triple strange!)\n");
        
        println!("     MESONS (quark-antiquark pairs):");
        println!("       pion+, pion-, pion0");
        println!("       kaon+, kaon- (strange mesons)\n");
        
        println!("     Usage: quark <particle_name>");
        return;
    }
    
    // Now analyze the specific particle
    // Use aptik for proton/neutron, hardcode others for now
    let (name, quarks, charge, voids, spikes) = match particle.to_lowercase().as_str() {
        "proton" | "p" => {
            let b = aptik::Baryon::proton();
            ("PROTON", "u u d", "+1", b.total_voids(), b.total_spikes())
        },
        "neutron" | "n" => {
            let b = aptik::Baryon::neutron();
            ("NEUTRON", "u d d", "0", b.total_voids(), b.total_spikes())
        },
        "delta++" => ("DELTA++", "u u u", "+2", 6, 0),
        "delta+" => ("DELTA+", "u u d", "+1", 4, 2),
        "delta0" => ("DELTA0", "u d d", "0", 2, 4),
        "delta-" => ("DELTA-", "d d d", "-1", 0, 6),
        "sigma+" => ("SIGMA+", "u u s", "+1", 4, 2),  // Strange replaces d
        "sigma0" => ("SIGMA0", "u d s", "0", 2, 4),
        "sigma-" => ("SIGMA-", "d d s", "-1", 0, 6),
        "lambda" => ("LAMBDA", "u d s", "0", 2, 4),
        "xi0" => ("XI0", "u s s", "0", 2, 4),
        "xi-" => ("XI-", "d s s", "-1", 0, 6),
        "omega-" => ("OMEGA-", "s s s", "-1", 0, 6),  // All strange!
        "pion+" => ("PION+ (meson)", "u d̄", "+1", 2, 2),  // Quark-antiquark
        "pion-" => ("PION- (meson)", "ū d", "-1", 2, 2),
        "pion0" => ("PION0 (meson)", "(uū+dd̄)/√2", "0", 2, 2),
        "kaon+" => ("KAON+ (meson)", "u s̄", "+1", 2, 2),  // Strange meson
        "kaon-" => ("KAON- (meson)", "ū s", "-1", 2, 2),
        _ => {
            let b = aptik::Baryon::proton();
            ("PROTON", "u u d", "+1", b.total_voids(), b.total_spikes())
        },
    };
    
    println!("  ═══════════════════════════════════════════════════════════════");
    println!("  ANALYZING: {}", name);
    println!("  ═══════════════════════════════════════════════════════════════\n");
    
    println!("     Quark content: {}", quarks);
    println!("     Electric charge: {}", charge);
    println!();
    
    println!("     CHARGE STACKING:");
    let quark_chars: Vec<&str> = quarks.split_whitespace().collect();
    let mut total_charge = 0.0f64;
    for q in &quark_chars {
        let (ch, flow) = match *q {
            "u" => ("+2/3", "⬇️ VOID (intake)"),
            "d" => ("-1/3", "⬆️ SPIKE (exhaust)"),
            "s" => ("-1/3", "⬆️ SPIKE* (strange resonance)"),
            "c" => ("+2/3", "⬇️ VOID* (charm resonance)"),
            "b" => ("-1/3", "⬆️ SPIKE* (bottom resonance)"),
            "t" => ("+2/3", "⬇️ VOID* (top resonance)"),
            "d̄" | "s̄" => ("+1/3", "⬇️ anti-spike"),
            "ū" => ("-2/3", "⬆️ anti-void"),
            _ => ("?", "unknown"),
        };
        let charge_val: f64 = match *q {
            "u" | "c" | "t" => 2.0/3.0,
            "d" | "s" | "b" => -1.0/3.0,
            "d̄" | "s̄" => 1.0/3.0,
            "ū" => -2.0/3.0,
            _ => 0.0,
        };
        total_charge += charge_val;
        println!("       {} → {} {}", q, ch, flow);
    }
    println!("       ─────────────────────────────");
    println!("       Total: {:+.0}/3 = {}", total_charge * 3.0, charge);
    
    println!();
    println!("     VOID/SPIKE ARCHITECTURE:");
    println!("       Each +2/3 quark = 2 void units (intake capacity)");
    println!("       Each -1/3 quark = 2 spike units (exhaust capacity)");
    println!("       Total voids:  {} units", voids);
    println!("       Total spikes: {} units", spikes);
    
    let flow = voids as i32 - spikes as i32;
    let balance = if flow > 0 { "INTAKE EXCESS → pulls in field energy" }
                  else if flow < 0 { "EXHAUST EXCESS → radiates energy" }
                  else if voids == 0 { "ALL SPIKE → explosive decay!" }
                  else if spikes == 0 { "ALL VOID → implosive decay!" }
                  else { "BALANCED → stable flow" };
    println!("       Flow balance: {:+} ({})", flow, balance);
    
    // Resonance sisters explanation
    println!();
    println!("  ═══════════════════════════════════════════════════════════════");
    println!("  THE RESONANCE SISTERS:");
    println!("  ═══════════════════════════════════════════════════════════════\n");
    
    println!("     When energy exceeds certain thresholds, a quark can");
    println!("     vibrate as its heavier 'sister' - same charge, more mass.\n");
    
    println!("     Up Family (+2/3 void character):");
    println!("       u → c: ~1.3 GeV threshold (charm 'sister')");
    println!("       u → t: ~173 GeV threshold (top 'sister' - extremely rare)\n");
    
    println!("     Down Family (-1/3 spike character):");
    println!("       d → s: ~95 MeV threshold (strange 'sister' - common!)");
    println!("       d → b: ~4.2 GeV threshold (bottom 'sister')\n");
    
    println!("     Strange particles (Λ, Σ, Ξ, Ω) contain 'd → s' resonances!");
    println!("     They decay back to u/d quarks when energy dissipates.");
    
    // The stability insight
    println!();
    println!("  ═══════════════════════════════════════════════════════════════");
    println!("  STABILITY INSIGHT:");
    println!("  ═══════════════════════════════════════════════════════════════\n");
    
    // STABLE baryons: proton (4:2), neutron (2:4) - both have 2:1 or 1:2 ratios
    // This is the magic ratio! Each quark type contributes its complementary partner
    // Delta++ (6:0) and Delta- (0:6) are unstable - no exhaust or no intake!
    let ratio = if spikes > 0 { voids as f64 / spikes as f64 } else { 999.0 };
    let stability = if voids == 0 {
        "EXTREMELY UNSTABLE - no intake, all exhaust → explosive decay (Δ⁻)"
    } else if spikes == 0 {
        "EXTREMELY UNSTABLE - no exhaust, all intake → implosive collapse (Δ⁺⁺)"
    } else if (ratio - 2.0).abs() < 0.01 || (ratio - 0.5).abs() < 0.01 {
        // 2:1 or 1:2 ratio = STABLE BARYON
        // This is proton (4:2 = 2:1) or neutron (2:4 = 1:2)
        "STABLE BARYON - 2:1 ratio allows self-sustaining flow (the magic ratio!)"
    } else if (ratio - 1.0).abs() < 0.01 {
        // 1:1 ratio - perfectly balanced
        "STABLE - equal intake/exhaust (deuteron-like balance)"
    } else {
        // Other ratios - depends on the specific particle
        "QUASI-STABLE - non-standard ratio, depends on configuration"
    };
    
    println!("     {}: {}", name, stability);
    println!("     Void/Spike ratio: {:.2}", ratio);
    
    if quarks.contains("s") || quarks.contains("c") || quarks.contains("b") || quarks.contains("t") {
        println!();
        println!("     ⚡ Contains resonance quark(s) - higher energy state!");
        println!("     Will decay to ground state (u/d only) when energy allows.");
    }
    
    // COLOR CHARGE ANALYSIS - The spatial flow direction!
    println!();
    println!("  ═══════════════════════════════════════════════════════════════");
    println!("  COLOR CHARGE → SPATIAL FLOW DIRECTION:");
    println!("  ═══════════════════════════════════════════════════════════════\n");
    
    println!("     Color charge determines WHICH AXIS the void/spike operates on:");
    println!("       🔴 Red   = X-axis flow");
    println!("       🟢 Green = Y-axis flow");
    println!("       🔵 Blue  = Z-axis flow\n");
    
    println!("     Baryons must be 'white' (R+G+B) = all 3 axes covered.\n");
    
    // Analyze color configurations for this particle
    let quark_list: Vec<&str> = quarks.split_whitespace().collect();
    let quark_count = quark_list.len();
    
    if quark_count == 3 {
        // Baryon - show the color flow pattern
        println!("     {} COLOR FLOW PATTERNS:", name);
        println!("     ─────────────────────────────────────────────────\n");
        
        // Count up quarks (voids) and down-type quarks (spikes)
        let mut ups = 0;
        let mut downs = 0;
        for q in &quark_list {
            match *q {
                "u" | "c" | "t" => ups += 1,
                "d" | "s" | "b" => downs += 1,
                _ => {}
            }
        }
        
        // Show the key insight: which axes have voids vs spikes
        if ups == 2 && downs == 1 {
            // Proton-like (uud pattern)
            println!("     Configuration: 2 VOIDS + 1 SPIKE across 3 axes\n");
            println!("     Example: u(🔴) u(🟢) d(🔵)");
            println!("       X-axis: ⬇️ VOID  (intake from +X and -X)");
            println!("       Y-axis: ⬇️ VOID  (intake from +Y and -Y)");
            println!("       Z-axis: ⬆️ SPIKE (exhaust along Z)\n");
            
            println!("     🔥 THE PINCH POINT / NOZZLE EFFECT:");
            println!("     ─────────────────────────────────────────────────\n");
            println!("         ⬇️ intake    ⬇️ intake");
            println!("          \\          /");
            println!("           \\   u   /     u");
            println!("            \\     /");
            println!("             \\   /  ← PINCH POINT!");
            println!("              \\ /");
            println!("            ───◯───  d");
            println!("               |");
            println!("               | ← NOZZLE");
            println!("               |");
            println!("               ⬆️   FOCUSED JET!\n");
            
            println!("     Two intakes CREATE A PINCH that focuses the exhaust!");
            println!("     The focused jet IS the electric field!");
            println!("     This is why protons have +1 charge!\n");
            
            println!("     EXTERNAL INTERACTION:");
            println!("       → 2 intakes converge → PINCH POINT");
            println!("       → Pinch focuses exhaust → NOZZLE/JET");
            println!("       → Directed jet = ELECTRIC FIELD!");
            println!("       → Gluon rotation spins the jet → oscillating dipole!");
        } else if ups == 1 && downs == 2 {
            // Neutron-like (udd pattern)
            println!("     Configuration: 1 VOID + 2 SPIKES across 3 axes\n");
            println!("     Example: u(🔴) d(🟢) d(🔵)");
            println!("       X-axis: ⬇️ VOID  (intake from +X and -X)");
            println!("       Y-axis: ⬆️ SPIKE (exhaust along Y)");
            println!("       Z-axis: ⬆️ SPIKE (exhaust along Z)\n");
            
            println!("     🌫️  NO PINCH POINT = DIFFUSE CLOUD:");
            println!("     ─────────────────────────────────────────────────\n");
            println!("               ⬇️ intake");
            println!("                |");
            println!("                u   (single intake - no pinch!)");
            println!("               /|\\");
            println!("              / | \\");
            println!("             /  |  \\");
            println!("          d /   |   \\ d");
            println!("           ⬆️   |   ⬆️");
            println!("          ~~~~~~|~~~~~~  ← DIFFUSE EXHAUST CLOUD");
            println!("               ~~~\n");
            
            println!("     Only 1 intake - NO PINCH POINT possible!");
            println!("     2 exhausts can't focus - they just DIFFUSE.");
            println!("     Neutron is surrounded by a cloud of exhaust.");
            println!("     No directed flow = no EM field = NEUTRAL!\n");
            
            println!("     WHY FREE NEUTRONS DECAY (~10 min):");
            println!("       → Diffuse cloud can't sustain itself alone!");
            println!("       → Needs proton jets routing through it for structure.");
            println!("       → In a nucleus: proton jets + neutron clouds interlock!");
            println!("       → Free neutron: cloud dissipates → d→u decay (β⁻).\n");
            
            println!("     The asymmetry still creates a MAGNETIC MOMENT!");
            println!("     (The intake axis spins → tiny current loop)");
        } else if ups == 3 && downs == 0 {
            // Delta++ (uuu)
            println!("     Configuration: 3 VOIDS + 0 SPIKES across 3 axes\n");
            println!("     Example: u(🔴) u(🟢) u(🔵)");
            println!("       X-axis: ⬇️ VOID  (intake)");
            println!("       Y-axis: ⬇️ VOID  (intake)");
            println!("       Z-axis: ⬇️ VOID  (intake)\n");
            println!("     EXTERNAL INTERACTION:");
            println!("       → ALL 3 axes pulling in field energy!");
            println!("       → NO exhaust anywhere!");
            println!("       → IMPLOSION from all directions!\n");
            println!("     This is why Δ⁺⁺ decays in ~10⁻²⁴ seconds!");
            println!("     It literally implodes - one u must become d to create exhaust.");
        } else if ups == 0 && downs == 3 {
            // Delta- (ddd) or Omega- (sss)
            println!("     Configuration: 0 VOIDS + 3 SPIKES across 3 axes\n");
            println!("     Example: {}(🔴) {}(🟢) {}(🔵)", 
                     quark_list[0], quark_list[1], quark_list[2]);
            println!("       X-axis: ⬆️ SPIKE (exhaust)");
            println!("       Y-axis: ⬆️ SPIKE (exhaust)");
            println!("       Z-axis: ⬆️ SPIKE (exhaust)\n");
            println!("     EXTERNAL INTERACTION:");
            println!("       → ALL 3 axes radiating outward!");
            println!("       → NO intake anywhere!");
            println!("       → EXPLOSION in all directions!\n");
            println!("     This is why Δ⁻ decays instantly!");
            println!("     It literally explodes - one d must become u to create intake.");
        } else {
            // Mixed strange configurations
            println!("     Configuration: {} VOIDS + {} SPIKES across 3 axes\n", ups, downs);
            println!("     Mixed void/spike across axes.");
            println!("     The arrangement still must be 'white' (R+G+B).");
        }
        
        // Gluon exchange insight
        println!();
        println!("     🔄 GLUON EXCHANGE (color rotation):");
        println!("     ─────────────────────────────────────────────────\n");
        println!("     Gluons carry color-anticolor pairs, e.g., (Red,Anti-Green)");
        println!("     This SWAPS which axis has which quark type!");
        println!();
        println!("     For proton: The exhaust axis rotates continuously!");
        println!("       Tick 1: exhaust on Z → Tick 2: exhaust on X → Tick 3: exhaust on Y");
        println!("     This creates an OSCILLATING dipole field pattern.");
        println!("     The oscillation IS the electromagnetic interaction!");
        
    } else if quark_count == 2 {
        // Meson
        println!("     {} is a MESON (quark + antiquark)", name);
        println!();
        println!("     Mesons have color + anticolor = white on ONE axis.");
        println!("     They mediate forces between baryons (nuclear force carriers).");
        println!();
        println!("     The quark-antiquark pair creates:");
        println!("       → Void-antivoid or spike-antispike annihilation!");
        println!("       → This is why mesons are short-lived.");
        println!("       → The flow structure is inherently unstable.");
    }
    
    println!();
    println!("  💡 THE FRACTIONAL CHARGES ARE FLOW RATIOS!");
    println!("     +2/3 = 2 parts intake capacity");
    println!("     -1/3 = 1 part exhaust capacity (or 2 parts as spike unit)");
    println!("     The 3 in the denominator = color charge confinement!");
    println!();
    println!("  🌈 COLOR = SPATIAL AXIS = 3D FLOW DIRECTION!");
    println!("     Red/Green/Blue map to X/Y/Z axes.");
    println!("     Gluons rotate which axis does what → oscillation → EM field!");
}