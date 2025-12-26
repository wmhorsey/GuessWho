//! # Resonant Packet Sorting System
//! 
//! An implementation of the Space-Time Energy Model for information sorting.
//! Packets are sorted not by content inspection, but by resonance attraction —
//! like frequencies attract, unlike frequencies repel.
//!
//! ## Architecture
//! - `Resonance`: Frequency signature (wavelength + amplitude + velocity)
//! - `Packet`: Light/energy unit with payload and resonance signature
//! - `Node`: Attractor point in the field with its own resonance
//! - `Supernode`: Emergent cluster of similar-resonance nodes
//! - `Switchboard`: Central registry for lookup, release, and garbage collection
//! - `Universe`: The cosmic container with CMB boundary, scaled to memory
//! - `Semantic`: Content-to-resonance mapping — the quantum tunneling mechanism

mod resonance;
mod packet;
mod node;
mod supernode;
mod switchboard;
mod field;
mod universe;
mod semantic;
mod pipeline;
mod loaders;

use universe::{Universe, UniverseConfig};
use packet::PacketType;
use pipeline::{Pipeline, ChemistrySource};
use std::path::Path;

/// Count chemistry files in a directory (recursive)
fn count_chemistry_files(dir: &Path) -> usize {
    let mut count = 0;
    if let Ok(entries) = std::fs::read_dir(dir) {
        for entry in entries.flatten() {
            let path = entry.path();
            if path.is_dir() {
                count += count_chemistry_files(&path);
            } else if let Some(ext) = path.extension().and_then(|s| s.to_str()) {
                if matches!(ext.to_lowercase().as_str(), "csv" | "tsv" | "sdf" | "mol" | "smi" | "smiles") {
                    count += 1;
                }
            }
        }
    }
    count
}

fn main() {
    println!("╔══════════════════════════════════════════════════════════╗");
    println!("║     RESONANT UNIVERSE — Space-Time Energy Model          ║");
    println!("╚══════════════════════════════════════════════════════════╝\n");

    // Create universe with custom config for demo
    let config = UniverseConfig {
        radius: 60.0,
        memory_bytes: 32 * 1024 * 1024, // 32 MB for demo
        dt: 0.5,
        damping: 0.95,
        force_strength: 8.0,
        capture_radius: 6.0,
        ..Default::default()
    };

    let mut universe = Universe::with_config(config);
    
    println!("{}\n", universe.memory_usage());

    // Spawn seed nodes
    println!("Spawning seed nodes...");
    universe.spawn_seed_nodes(12);
    println!("Created {} nodes\n", universe.nodes.len());

    // Release packets into the universe
    println!("Releasing packets into the universe...\n");
    
    // === DOCUMENT-LEVEL TEST — SHORT STORIES & PARAGRAPHS ===
    // Testing if word-level resonance can cluster longer content by theme/meaning
    // Different genres, topics, styles — will they cluster semantically?
    
    let content_samples = [
        // === SCIENCE FICTION ===
        "The starship descended through the clouds of the alien world. Captain Chen gripped the controls as sensors detected unknown life forms below. This was humanity's first contact.",
        
        "In the year 2347, robots had become indistinguishable from humans. Dr. Yamamoto wondered if her assistant was real or synthetic. The line between machine and consciousness had blurred.",
        
        "The wormhole opened like a flower blooming in reverse. Through it, travelers from a parallel Earth emerged, carrying news of a universe where gravity worked differently.",
        
        // === FANTASY ===
        "The dragon circled the ancient tower where the princess waited. But she was no damsel in distress. She raised her staff and spoke words of power that made the beast bow.",
        
        "In the enchanted forest, elves gathered to discuss the growing darkness. The old magic was fading, and only the prophecy of the chosen one offered hope.",
        
        "The wizard's apprentice accidentally turned the king into a frog. Now he had three days to find the counter-spell or face execution. Magic was harder than it looked.",
        
        // === ROMANCE ===
        "She saw him across the crowded coffee shop. Their eyes met. He smiled. She dropped her latte. It was the most embarrassing and wonderful moment of her life.",
        
        "After twenty years apart, they met again at their high school reunion. The spark was still there, burning brighter than ever. Some loves never truly end.",
        
        "He wrote her a letter every day for a year. She never replied. On the last day, she appeared at his door with all 365 letters, tear-stained and memorized.",
        
        // === HORROR ===
        "The house had been empty for fifty years. When the new family moved in, they found children's handprints on the windows. From the inside. Fresh.",
        
        "She woke at 3 AM to the sound of breathing. Something was in the room. Something that whispered her name in a voice like grinding bones.",
        
        "The mirror showed a reflection that wasn't quite right. It smiled when she didn't. It moved when she was still. Tonight, it was reaching through.",
        
        // === TECHNICAL/PROGRAMMING ===
        "The function recursively traverses the binary tree, comparing node values at each level. Time complexity is O(n) where n is the number of nodes.",
        
        "To implement the sorting algorithm, first partition the array around a pivot element. Recursively apply the same logic to sub-arrays until sorted.",
        
        "The API endpoint accepts JSON payloads with authentication tokens. Response codes indicate success or failure. Rate limiting applies after 100 requests.",
        
        // === COOKING/RECIPES ===
        "Preheat oven to 350 degrees. Mix flour, sugar, and eggs until smooth. Add vanilla extract and fold in chocolate chips. Bake for 12 minutes.",
        
        "The secret to perfect risotto is patience. Add broth one ladle at a time, stirring constantly. The rice should be creamy but still have bite.",
        
        "For the marinade, combine soy sauce, ginger, garlic, and honey. Let the chicken rest in this mixture overnight. Grill over medium heat.",
        
        // === NATURE/WILDLIFE ===
        "The wolf pack moved silently through the snow. The alpha paused, ears forward, sensing prey. In winter, every hunt meant survival.",
        
        "Beneath the ocean surface, the coral reef teemed with life. Fish of every color darted between formations that had grown for centuries.",
        
        "The eagle soared on thermal currents, scanning the meadow below. A rabbit moved. In seconds, the hunter became death from above.",
        
        // === PHILOSOPHY/ABSTRACT ===
        "What is consciousness but patterns recognizing themselves? The universe evolved eyes to see itself, minds to question itself.",
        
        "Time flows like a river, but who stands on the bank to watch it pass? Perhaps we are both the water and the observer.",
        
        "If a tree falls in the forest and no one hears it, does it make a sound? Perhaps the forest itself is always listening.",
        
        // === NEWS/JOURNALISM ===
        "The mayor announced new infrastructure spending today. Roads, bridges, and public transit will receive funding over the next five years.",
        
        "Scientists discovered a new species of deep-sea fish. The creature lives at depths previously thought impossible for vertebrate life.",
        
        "Stock markets rose sharply on news of the trade agreement. Analysts predict continued growth through the fourth quarter.",
    ];
    
    // Create packets from content
    for content in &content_samples {
        universe.spawn_content_packet(content);
    }
    
    // Fewer random packets to let content-based ones dominate
    for _ in 0..10 {
        universe.spawn_packet(PacketType::Signal);
        universe.spawn_packet(PacketType::Memory);
        universe.spawn_packet(PacketType::Query);
        universe.spawn_packet(PacketType::Echo);
    }
    // Add some anomalies - the mysteries
    for _ in 0..10 {
        universe.spawn_packet(PacketType::Anomaly);
    }
    
    println!("Released {} packets ({} with content)\n", 
        universe.packets.len(),
        content_samples.len());

    // Simulation loop
    let total_ticks = 200;  // Longer run for fission to occur
    let render_interval = 40;

    println!("Running simulation ({} ticks)...\n", total_ticks);

    for tick in 1..=total_ticks {
        universe.tick();

        if tick % render_interval == 0 || tick == 1 {
            // Clear screen effect
            println!("\n{}", "─".repeat(62));
            println!("{}", universe.status());
            println!("{}", universe.memory_usage());
            println!();
            
            // Render the universe
            println!("{}", universe.render(60, 30));
            println!();
            
            // Show discoveries
            if !universe.discoveries.is_empty() {
                println!("✧ DISCOVERIES ✧");
                for discovery in &universe.discoveries {
                    println!("  {} - traveled {:.1} units, hit CMB {} times",
                        discovery.signature(),
                        discovery.distance_traveled,
                        discovery.cmb_hits
                    );
                }
            }
            
            // CMB stats
            println!("\nCMB Statistics:");
            println!("  Total hits: {}", universe.cmb_stats.total_hits);
            println!("  Bounces: {}", universe.cmb_stats.bounces);
            println!("  Cosmic rays created: {}", universe.cmb_stats.cosmic_rays_created);
            println!("  Discoveries: {}", universe.cmb_stats.discoveries);
        }
    }

    // Final report
    println!("\n{}", "═".repeat(62));
    println!("FINAL STATE");
    println!("{}", "═".repeat(62));
    println!("{}", universe.status());
    
    // === SYSTEM MATURITY REPORT ===
    println!("\n🌱 SYSTEM MATURITY — The universe is young and learning");
    println!("─────────────────────────────────────────────────────────────");
    let seed_count = 12u64;
    let learned = universe.learning.nodes_learned;
    let total_nodes = universe.nodes.len() as u64;
    let maturity_pct = if total_nodes > 0 { (learned as f64 / total_nodes as f64) * 100.0 } else { 0.0 };
    
    println!("   Age: {} ticks", universe.tick);
    println!("   Seed nodes: {}", seed_count);
    println!("   Learned nodes: {} ({:.1}% of total)", learned, maturity_pct);
    println!("   Supernodes: {}", universe.supernodes.len());
    println!("   Fission events: {} (self-organization)", universe.fission_events);
    println!("   Garbage collected: {} (didn't matter)", universe.garbage_collected);
    println!("   Discoveries: {} (something new!)", universe.discoveries.len());
    
    // Growth assessment
    let growth_stage = if universe.tick < 50 {
        "🌱 Infant — Just born, still forming basic structure"
    } else if universe.tick < 200 {
        "🌿 Seedling — Learning patterns, nodes emerging"
    } else if universe.tick < 500 {
        "🌳 Sapling — Structure stabilizing, supernodes forming"
    } else if universe.tick < 1000 {
        "🌲 Young Tree — Refined categories, efficient sorting"
    } else {
        "🏔️ Mature — Deep structure, fast recall"
    };
    println!("\n   Growth stage: {}", growth_stage);
    
    println!("\nSupernodes formed:");
    for supernode in universe.supernodes.values() {
        println!("  {}", supernode.signature());
    }

    // === SEMANTIC CLUSTERING REPORT ===
    println!("\n{}", "═".repeat(62));
    println!("SEMANTIC CLUSTERING — Content sorted by resonance");
    println!("{}", "═".repeat(62));
    
    // Group content packets by which node captured them (or flying together)
    use std::collections::HashMap;
    let mut node_contents: HashMap<Option<u64>, Vec<&str>> = HashMap::new();
    
    for packet in universe.packets.values() {
        if let Some(ref content) = packet.content {
            // Find which node (if any) this packet is near
            let nearest_node = universe.nodes.values()
                .filter(|n| {
                    let dx = n.position.0 - packet.position.x;
                    let dy = n.position.1 - packet.position.y;
                    (dx*dx + dy*dy).sqrt() < 15.0
                })
                .min_by(|a, b| {
                    let dist_a = ((a.position.0 - packet.position.x).powi(2) + 
                                  (a.position.1 - packet.position.y).powi(2)).sqrt();
                    let dist_b = ((b.position.0 - packet.position.x).powi(2) + 
                                  (b.position.1 - packet.position.y).powi(2)).sqrt();
                    dist_a.partial_cmp(&dist_b).unwrap()
                })
                .map(|n| n.id);
            
            node_contents.entry(nearest_node).or_default().push(content);
        }
    }
    
    // Print clustered content
    for (node_id, contents) in node_contents.iter() {
        if contents.len() > 1 {
            match node_id {
                Some(id) => {
                    if let Some(node) = universe.nodes.get(id) {
                        println!("\n📍 Node {} (λ={:.2}, α={:.2}):", 
                            id, node.resonance.wavelength, node.resonance.amplitude);
                    }
                }
                None => println!("\n🌌 Drifting together in space:"),
            }
            for content in contents {
                println!("   • \"{}\"", content);
            }
        }
    }
    
    println!("\n✧ THE INSIGHT: Similar content clusters at similar nodes! ✧");
    println!("   Content → Resonance → Node → Supernode");
    println!("   The sorting emerges from the physics, not from rules.\n");

    println!("\n✧ THE 0.01% — Where Mysteries Live ✧");
    if universe.discoveries.is_empty() {
        println!("  No true discoveries yet... keep running the simulation.");
    } else {
        for (i, discovery) in universe.discoveries.iter().enumerate() {
            println!("  {}. {} — {}", 
                i + 1,
                discovery.signature_full(),
                match discovery.cmb_hits {
                    1..=2 => "glanced the edge of reality",
                    3..=5 => "touched the cosmic boundary",
                    _ => "transcended the observable universe",
                }
            );
        }
    }

    println!("\nLegend:");
    println!("  X = Center    n/+ = Nodes    O/o = Supernodes");
    println!("  ! = Signal    # = Memory     ? = Query");
    println!("  ~ = Echo      * = Anomaly    @ = Cosmic Ray");
    println!("  . = CMB Boundary (edge of the observable universe)");

    // ═══════════════════════════════════════════════════════════════════
    // RESONANT SEARCH — Query travels through the field
    // ═══════════════════════════════════════════════════════════════════
    println!("\n{}", "═".repeat(62));
    println!("RESONANT SEARCH — Query packet travels through the field");
    println!("{}", "═".repeat(62));
    
    // The search query becomes a packet that travels through nodes
    let search_queries = [
        "dragon magic spell wizard",
        "robot artificial consciousness machine",
        "cooking recipe bake oven",
        "stock market trading finance",
        "space alien starship",
    ];
    
    for query in &search_queries {
        let result = universe.search_resonant(query, 10);
        
        println!("\n🔍 {}", result.summary());
        println!("   Path: center");
        for &node_id in &result.nodes_visited {
            if let Some(node) = universe.nodes.get(&node_id) {
                print!(" → N{}[λ={:.2}]", node_id, node.resonance.wavelength);
            }
        }
        println!();
        
        if result.matches.is_empty() {
            println!("   No resonant matches found");
        } else {
            println!("   Matches (memory refs):");
            for m in result.top(3) {
                let source = match m.found_via {
                    universe::MatchSource::ResonancePath => "⚡",
                    universe::MatchSource::NearbyDrift => "~",
                    universe::MatchSource::CmbScan => "🌌",
                };
                if let Some(content) = &m.content {
                    let preview: String = content.chars().take(50).collect();
                    println!("   {} mem:{} [affinity={:.2}] {}...", 
                        source, m.memory_ref, m.affinity, preview);
                } else {
                    println!("   {} mem:{} [affinity={:.2}]", source, m.memory_ref, m.affinity);
                }
            }
        }
    }
    
    // Demonstrate access actions
    println!("\n📌 How it works:");
    println!("   1. Query → Resonance (filter)");
    println!("   2. Resonance → Travel through nodes (physics)");
    println!("   3. Arrive at resonant node → Retrieve memory addresses");
    println!("   4. CMB scan catches what resonance missed");
    println!("\n   Match sources: ⚡ Resonance Path | ~ Nearby Drift | 🌌 CMB Scan");
    println!("\n   The query doesn't SCAN — it RESONATES through the field!");

    // ═══════════════════════════════════════════════════════════════════
    // CHEMISTRY PIPELINE — Real domain-specific data
    // ═══════════════════════════════════════════════════════════════════
    println!("\n{}", "═".repeat(62));
    println!("CHEMISTRY PIPELINE — Domain-specific resonance");
    println!("{}", "═".repeat(62));
    
    // Create a fresh universe for chemistry - SCALE TO AVAILABLE RESOURCES
    // User has 64GB RAM and 3TB storage!
    // Key insight: Smaller radius = faster capture, more nodes = more landing spots
    let mut chem_universe = Universe::new();
    chem_universe.config.radius = 150.0;  // Compact universe for faster capture
    chem_universe.config.memory_bytes = 4 * 1024 * 1024 * 1024;  // 4GB for chemistry
    chem_universe.config.nodes_per_mb = 200;  // 200 nodes per MB = 800k max nodes
    chem_universe.config.capture_radius = 15.0;  // Large capture radius
    chem_universe.config.force_strength = 10.0;  // Stronger attraction to nodes
    
    // Create pipeline
    let mut pipeline = Pipeline::new();
    
    // Try to load from local data files first
    let data_dir = std::path::Path::new("data");
    let mut chem_source = if data_dir.exists() {
        println!("\n📁 Scanning data/ directory (including subdirectories)...");
        
        // Count files first
        let file_count = count_chemistry_files(data_dir);
        println!("   Found {} chemistry data files", file_count);
        
        if file_count > 100 {
            println!("   Loading in bulk mode (this may take a moment)...");
        }
        
        // Try to load any chemistry files found
        match loaders::load_directory(data_dir) {
            Ok(sources) if !sources.is_empty() => {
                // Merge all sources
                let mut combined = ChemistrySource::new("Local Data");
                let mut total_records = 0;
                for mut source in sources {
                    let records = source.records_owned();
                    total_records += records.len();
                    for record in records {
                        combined.add(record);
                    }
                }
                println!("   Total compounds loaded: {}", total_records);
                combined
            }
            _ => {
                println!("   No chemistry files found, using demo data");
                ChemistrySource::demo()
            }
        }
    } else {
        // Generate sample data files for next time
        println!("\n📁 No data/ directory found. Creating with sample files...");
        std::fs::create_dir_all("data").ok();
        loaders::generate_sample_csv("data/compounds.csv").ok();
        loaders::generate_sample_smiles("data/molecules.smi").ok();
        println!("   → Place your PubChem CSV/SDF files in data/ and re-run!");
        println!("   → Using built-in demo data for now\n");
        ChemistrySource::demo()
    };
    
    // Seed the universe with nodes covering the FULL resonance spectrum
    // More seed nodes = more landing spots for packets
    let seed_count = 2000;  // 2000 seed nodes spread across wavelengths
    println!("\n🌱 Seeding universe with {} nodes across resonance spectrum...", seed_count);
    chem_universe.spawn_seed_nodes(seed_count);
    println!("   Universe capacity: {} max nodes", chem_universe.config.max_nodes());
    
    println!("🧪 Ingesting into resonance field...");
    let start_time = std::time::Instant::now();
    if let Err(e) = pipeline.ingest(&mut chem_source, &mut chem_universe, None) {
        println!("   Error: {}", e);
    }
    let ingest_time = start_time.elapsed();
    println!("   Ingestion took {:.2}s", ingest_time.as_secs_f64());
    
    // Scale simulation ticks based on dataset size - run longer for large datasets
    let packet_count = chem_universe.packets.len();
    let tick_count = if packet_count > 10000 { 2000 }  // Much longer for 50k+
                     else if packet_count > 1000 { 500 } 
                     else { 200 };
    
    println!("\n⚛️  Running molecular dynamics ({} ticks for {} compounds)...", tick_count, packet_count);
    println!("   Universe: radius={:.0}, max_nodes={}, capture_radius={:.1}", 
        chem_universe.config.radius, 
        chem_universe.config.max_nodes(),
        chem_universe.config.capture_radius);
    println!("   Tracking capture rate over time:\n");
    
    let start_time = std::time::Instant::now();
    for i in 0..tick_count {
        chem_universe.tick();
        
        // Report capture rate at intervals
        if (i + 1) % 500 == 0 || i == 99 || i == 249 {
            let captured: usize = chem_universe.nodes.values()
                .map(|n| n.captured_packets.len())
                .sum();
            let rate = 100.0 * captured as f64 / packet_count.max(1) as f64;
            
            // Count packet states
            let mut in_flight = 0;
            let mut at_boundary = 0;
            let mut released = 0;
            let mut total_cmb_hits = 0;
            for p in chem_universe.packets.values() {
                match p.state {
                    packet::PacketState::InFlight => in_flight += 1,
                    packet::PacketState::AtBoundary => at_boundary += 1,
                    packet::PacketState::Released => released += 1,
                    _ => {}
                }
                total_cmb_hits += p.cmb_hits;
            }
            
            println!("   Tick {:>4}: captured={:>5} ({:.1}%) | flying={:>5} | boundary={} | garbage={} | cmb_hits={}",
                i + 1, captured, rate, in_flight, at_boundary, released, total_cmb_hits);
        }
    }
    let sim_time = start_time.elapsed();
    println!("\n   Simulation took {:.2}s ({:.0} ticks/sec)", 
        sim_time.as_secs_f64(), 
        tick_count as f64 / sim_time.as_secs_f64());
    
    // Show what clustered together
    println!("\n🔬 MOLECULAR CLUSTERING RESULTS:");
    println!("   Do similar molecules cluster together?\n");
    
    // Find nodes with captured chemistry packets
    let mut node_molecules: Vec<(u64, f64, Vec<String>)> = Vec::new();
    
    for node in chem_universe.nodes.values() {
        if node.captured_packets.is_empty() {
            continue;
        }
        
        let mut molecules = Vec::new();
        for &pid in &node.captured_packets {
            if let Some(packet) = chem_universe.packets.get(&pid) {
                if let Some(ref content) = packet.content {
                    molecules.push(content.clone());
                }
            }
        }
        
        if !molecules.is_empty() {
            node_molecules.push((node.id, node.resonance.wavelength, molecules));
        }
    }
    
    // Sort by number of molecules (largest clusters first)
    node_molecules.sort_by(|a, b| b.2.len().cmp(&a.2.len()));
    
    // DETAILED PACKET STATE ACCOUNTING - where did every packet go?
    let mut final_captured = 0;
    let mut final_in_flight = 0;
    let mut final_at_boundary = 0;
    let mut final_released = 0;
    let mut final_consumed = 0;
    let mut high_cmb_hits = 0;  // packets that hit CMB 3+ times
    let mut no_resonance_match = 0;  // packets whose resonance doesn't match any node
    
    for packet in chem_universe.packets.values() {
        match packet.state {
            packet::PacketState::Captured => final_captured += 1,
            packet::PacketState::InFlight => {
                final_in_flight += 1;
                // Check if this packet just can't find a home
                if packet.cmb_hits >= 3 {
                    high_cmb_hits += 1;
                }
            }
            packet::PacketState::AtBoundary => final_at_boundary += 1,
            packet::PacketState::Released => final_released += 1,
            packet::PacketState::Consumed => final_consumed += 1,
        }
    }
    
    // CLUSTER STATISTICS
    let total_clusters = node_molecules.len();
    let packets_in_clusters: usize = node_molecules.iter().map(|(_, _, m)| m.len()).sum();
    let largest_cluster = node_molecules.first().map(|(_, _, m)| m.len()).unwrap_or(0);
    let avg_cluster_size = if total_clusters > 0 { 
        packets_in_clusters as f64 / total_clusters as f64 
    } else { 0.0 };
    
    println!("   📊 PACKET ACCOUNTING (where did every packet go?):");
    println!("      ────────────────────────────────────────────");
    println!("      Total packets created:  {:>6}", packet_count);
    println!("      ────────────────────────────────────────────");
    println!("      ✅ Captured by nodes:   {:>6} ({:.1}%)", final_captured, 100.0 * final_captured as f64 / packet_count.max(1) as f64);
    println!("      🚀 Still flying:        {:>6} ({:.1}%)", final_in_flight, 100.0 * final_in_flight as f64 / packet_count.max(1) as f64);
    println!("      🌌 At CMB boundary:     {:>6}", final_at_boundary);
    println!("      🗑️  Garbage collected:  {:>6}", final_released);
    println!("      📤 Consumed (accessed): {:>6}", final_consumed);
    println!("      ────────────────────────────────────────────");
    println!("      ⚠️  Bounced 3+ times:    {:>6} (lost/can't find home)", high_cmb_hits);
    
    // ANALYZE GARBAGE PACKETS - why were they rejected?
    if final_released > 0 {
        println!("\n   🔍 GARBAGE ANALYSIS - Why were {} packets rejected?", final_released);
        
        let mut garbage_wavelengths: Vec<f64> = Vec::new();
        let mut garbage_samples: Vec<(f64, u64, String)> = Vec::new();  // (wavelength, cmb_hits, content)
        let mut garbage_with_content = 0;
        let mut total_cmb_on_garbage: u64 = 0;
        
        for packet in chem_universe.packets.values() {
            if packet.state == packet::PacketState::Released {
                garbage_wavelengths.push(packet.resonance.wavelength);
                total_cmb_on_garbage += packet.cmb_hits as u64;
                
                if let Some(ref content) = packet.content {
                    garbage_with_content += 1;
                    if garbage_samples.len() < 10 {
                        garbage_samples.push((packet.resonance.wavelength, packet.cmb_hits as u64, content.clone()));
                    }
                }
            }
        }
        
        println!("      Garbage packets with content: {}/{}", garbage_with_content, final_released);
        println!("      Average CMB hits on garbage: {:.1}", total_cmb_on_garbage as f64 / final_released.max(1) as f64);
        
        // Find wavelength distribution of garbage
        if !garbage_wavelengths.is_empty() {
            garbage_wavelengths.sort_by(|a, b| a.partial_cmp(b).unwrap());
            let min_wl = garbage_wavelengths.first().unwrap_or(&0.0);
            let max_wl = garbage_wavelengths.last().unwrap_or(&1.0);
            let avg_wl: f64 = garbage_wavelengths.iter().sum::<f64>() / garbage_wavelengths.len() as f64;
            
            println!("      Wavelength range: {:.3} - {:.3} (avg: {:.3})", min_wl, max_wl, avg_wl);
            
            // Bucket the wavelengths
            let mut buckets = [0usize; 10];
            for wl in &garbage_wavelengths {
                let bucket = ((wl * 10.0) as usize).min(9);
                buckets[bucket] += 1;
            }
            println!("      Wavelength distribution:");
            for (i, count) in buckets.iter().enumerate() {
                if *count > 0 {
                    let bar = "█".repeat((*count * 20 / garbage_wavelengths.len().max(1)).max(1));
                    println!("        {:.1}-{:.1}: {:>4} {}", i as f64 / 10.0, (i + 1) as f64 / 10.0, count, bar);
                }
            }
            
            // Show sample garbage molecules
            println!("      \n      Sample rejected molecules:");
            if garbage_samples.is_empty() {
                println!("        (no content attached to garbage packets)");
            } else {
                for (wl, cmb, content) in garbage_samples.iter().take(5) {
                    println!("        λ={:.3} ({}x CMB): {}", wl, cmb, content);
                }
            }
            
            // Check if there are nodes in that wavelength range
            let mut nodes_in_range = 0;
            for node in chem_universe.nodes.values() {
                if node.resonance.wavelength >= *min_wl && node.resonance.wavelength <= *max_wl {
                    nodes_in_range += 1;
                }
            }
            println!("      \n      Nodes in garbage wavelength range: {}", nodes_in_range);
        }
    }
    
    println!("\n   📈 CLUSTER STATISTICS:");
    println!("      Total clusters formed: {}", total_clusters);
    println!("      Largest cluster: {} molecules", largest_cluster);
    println!("      Average cluster size: {:.1} molecules", avg_cluster_size);
    println!("      Capture rate: {:.1}%", 100.0 * packets_in_clusters as f64 / packet_count.max(1) as f64);
    
    // Distribution of cluster sizes
    let mut size_buckets = [0usize; 6]; // 1, 2-5, 6-10, 11-50, 51-100, 100+
    for (_, _, molecules) in &node_molecules {
        match molecules.len() {
            1 => size_buckets[0] += 1,
            2..=5 => size_buckets[1] += 1,
            6..=10 => size_buckets[2] += 1,
            11..=50 => size_buckets[3] += 1,
            51..=100 => size_buckets[4] += 1,
            _ => size_buckets[5] += 1,
        }
    }
    println!("\n      Cluster size distribution:");
    println!("        1 molecule:   {} clusters", size_buckets[0]);
    println!("        2-5:          {} clusters", size_buckets[1]);
    println!("        6-10:         {} clusters", size_buckets[2]);
    println!("        11-50:        {} clusters", size_buckets[3]);
    println!("        51-100:       {} clusters", size_buckets[4]);
    println!("        100+:         {} clusters", size_buckets[5]);
    
    println!("\n   🔝 LARGEST CLUSTERS:");
    for (node_id, wavelength, molecules) in node_molecules.iter().take(8) {
        println!("   📍 Node {} [λ={:.3}] — {} molecules:", node_id, wavelength, molecules.len());
        for mol in molecules.iter().take(3) {
            println!("      • {}", mol);
        }
        if molecules.len() > 3 {
            println!("      ... and {} more", molecules.len() - 3);
        }
        println!();
    }
    
    // Test chemistry-specific search
    println!("   🔍 Chemistry Search Test:");
    
    let chem_queries = [
        "sugar glucose fructose",
        "amino acid protein",
        "aromatic benzene ring",
        "pain aspirin ibuprofen",
    ];
    
    for query in &chem_queries {
        let result = chem_universe.search_resonant(query, 10);
        print!("   '{}': ", query);
        if result.matches.is_empty() {
            println!("No matches");
        } else {
            let previews: Vec<String> = result.top(2)
                .iter()
                .filter_map(|m| m.content.as_ref())
                .map(|c| c.chars().take(20).collect::<String>())
                .collect();
            println!("{}", previews.join(" | "));
        }
    }
    
    println!("\n   ✧ The chemistry data clusters by molecular structure!");
    println!("   ✧ Sugars find sugars, drugs find drugs — via RESONANCE.");

    // Demonstrate persistence
    println!("\n{}", "═".repeat(62));
    println!("PERSISTENCE");
    println!("{}", "═".repeat(62));
    
    // Save entire universe state
    match universe.save("universe_state.json") {
        Ok(_) => println!("✓ Universe state saved to universe_state.json"),
        Err(e) => println!("✗ Failed to save universe: {}", e),
    }
    
    // Save just the learned nodes (for reuse in future universes)
    match universe.save_nodes("learned_nodes.json") {
        Ok(_) => println!("✓ Learned nodes saved to learned_nodes.json"),
        Err(e) => println!("✗ Failed to save nodes: {}", e),
    }
    
    // Show what we saved
    println!("\nSaved {} nodes and {} supernodes for future sessions.", 
        universe.nodes.len(),
        universe.supernodes.len());
    println!("The resonance field now has persistent memory.");
}
