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

use universe::{Universe, UniverseConfig};
use packet::{Packet, PacketType};

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
    // EVENT-DRIVEN ACCESS DEMONSTRATION
    // ═══════════════════════════════════════════════════════════════════
    println!("\n{}", "═".repeat(62));
    println!("EVENT-DRIVEN ACCESS — Search triggers consumption");
    println!("{}", "═".repeat(62));
    
    // Search by content - this is how real access would work
    let search_queries = [
        "dragon magic spell",
        "robot artificial consciousness",
        "cooking recipe bake",
        "stock market trading",
    ];
    
    for query in &search_queries {
        println!("\n🔍 Searching: \"{}\"", query);
        let results = universe.search_by_content(query);
        
        if results.is_empty() {
            println!("   No resonant matches found");
        } else {
            println!("   Found {} resonant packets:", results.len().min(3));
            for (i, packet) in results.iter().take(3).enumerate() {
                if let Some(content) = &packet.content {
                    // Truncate long content
                    let preview: String = content.chars().take(60).collect();
                    println!("   {}. [λ={:.3}] {}...", i+1, packet.resonance.wavelength, preview);
                }
            }
        }
    }
    
    // Demonstrate access actions
    println!("\n📌 Access Actions:");
    println!("   • Lock   → Data becomes permanent, node strengthens");
    println!("   • Boost  → Node grows (captures more similar content)");
    println!("   • Reject → Data drifts away, node weakens");
    println!("   • Consume → Data used and released");

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
