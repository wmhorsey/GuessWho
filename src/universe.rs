//! Universe — The cosmic container for the resonant field
//!
//! The universe has a finite radius (the CMB boundary) scaled to available memory.
//! Packets that hit the boundary either bounce back as cosmic rays or escape
//! as discoveries. The universe renders itself for visualization.

use std::collections::HashMap;
use std::path::Path;
use std::fs;
use serde::{Serialize, Deserialize};
use crate::packet::{Packet, PacketId, PacketState, PacketType, Position};
use crate::node::{Node, NodeId};
use crate::supernode::{Supernode, SupernodeId};
use crate::resonance::Resonance;

/// Configuration for the universe
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct UniverseConfig {
    /// Radius of the observable universe (CMB boundary)
    pub radius: f64,
    /// Available memory in bytes (determines capacity)
    pub memory_bytes: usize,
    /// Packets per MB of memory
    pub packets_per_mb: usize,
    /// Nodes per MB of memory  
    pub nodes_per_mb: usize,
    /// Time step for physics simulation
    pub dt: f64,
    /// Damping factor for velocity
    pub damping: f64,
    /// Force strength multiplier
    pub force_strength: f64,
    /// Capture radius for nodes
    pub capture_radius: f64,
}

impl Default for UniverseConfig {
    fn default() -> Self {
        Self {
            radius: 100.0,
            memory_bytes: 64 * 1024 * 1024, // 64 MB default
            packets_per_mb: 1000,
            nodes_per_mb: 10,
            dt: 1.0,
            damping: 0.98,
            force_strength: 5.0,
            capture_radius: 5.0,
        }
    }
}

impl UniverseConfig {
    /// Create config scaled to system memory
    pub fn from_available_memory() -> Self {
        // Try to detect available memory (fallback to 64MB)
        let memory_bytes = detect_available_memory().unwrap_or(64 * 1024 * 1024);
        
        // Scale radius with memory (sqrt scaling for 2D area feel)
        let mb = memory_bytes as f64 / (1024.0 * 1024.0);
        let radius = 50.0 + (mb.sqrt() * 10.0);
        
        Self {
            radius,
            memory_bytes,
            ..Default::default()
        }
    }

    pub fn max_packets(&self) -> usize {
        (self.memory_bytes / (1024 * 1024)) * self.packets_per_mb
    }

    pub fn max_nodes(&self) -> usize {
        (self.memory_bytes / (1024 * 1024)) * self.nodes_per_mb
    }
}

/// Statistics about CMB interactions
#[derive(Default, Debug, Clone, Serialize, Deserialize)]
pub struct CMBStats {
    /// Total packets that hit the CMB
    pub total_hits: u64,
    /// Packets that bounced back
    pub bounces: u64,
    /// Packets that became cosmic rays
    pub cosmic_rays_created: u64,
    /// True discoveries (packets that defied all classification)
    pub discoveries: u64,
}

/// Actions that can be taken when accessing data
#[derive(Debug, Clone, Copy)]
pub enum AccessAction {
    /// Lock the data permanently at its node
    Lock,
    /// Consume and release the data
    Consume,
    /// Boost the node (data was useful)
    Boost,
    /// Reject - data was wrong, let it drift
    Reject,
}

/// Result of an access operation
#[derive(Debug)]
pub enum AccessResult {
    /// Data locked at a node
    Locked { node_id: Option<NodeId> },
    /// Data consumed
    Consumed,
    /// Node boosted
    Boosted { node_id: Option<NodeId> },
    /// Data rejected, drifting away
    Rejected,
}

/// A single match from a resonant search
#[derive(Debug, Clone)]
pub struct SearchMatch {
    /// The packet that matched
    pub packet_id: PacketId,
    /// Memory reference (the actual data address)
    pub memory_ref: u64,
    /// How strongly it resonated with the query
    pub affinity: f64,
    /// Which node it was found at (0 if drifting)
    pub node_id: NodeId,
    /// Content preview (if available)
    pub content: Option<String>,
    /// How the match was found
    pub found_via: MatchSource,
}

/// How a search match was discovered
#[derive(Debug, Clone, PartialEq)]
pub enum MatchSource {
    /// Found via resonance path (query traveled to this node)
    ResonancePath,
    /// Found nearby while traveling
    NearbyDrift,
    /// Found via CMB boundary scan (keyword match)
    CmbScan,
}

/// Result of a resonant search
#[derive(Debug)]
pub struct SearchResult {
    /// The original query
    pub query: String,
    /// The query's resonance signature
    pub query_resonance: Resonance,
    /// How many hops the query took through the field
    pub hops: usize,
    /// Nodes visited during the search
    pub nodes_visited: Vec<NodeId>,
    /// Matches found, ordered by affinity
    pub matches: Vec<SearchMatch>,
}

impl SearchResult {
    /// Get the top N matches
    pub fn top(&self, n: usize) -> &[SearchMatch] {
        &self.matches[..self.matches.len().min(n)]
    }
    
    /// Get memory references for all matches
    pub fn memory_refs(&self) -> Vec<u64> {
        self.matches.iter().map(|m| m.memory_ref).collect()
    }
    
    /// Display a summary
    pub fn summary(&self) -> String {
        format!(
            "Query '{}' [λ={:.3}] → {} hops, {} nodes, {} matches",
            self.query,
            self.query_resonance.wavelength,
            self.hops,
            self.nodes_visited.len(),
            self.matches.len()
        )
    }
}

/// The observable universe
#[derive(Serialize, Deserialize)]
pub struct Universe {
    /// Configuration
    pub config: UniverseConfig,
    /// All packets in the universe
    pub packets: HashMap<PacketId, Packet>,
    /// All nodes in the universe
    pub nodes: HashMap<NodeId, Node>,
    /// Supernodes (emergent clusters)
    pub supernodes: HashMap<SupernodeId, Supernode>,
    /// Current tick
    pub tick: u64,
    /// CMB statistics
    pub cmb_stats: CMBStats,
    /// Discoveries (packets that escaped or became anomalies)
    pub discoveries: Vec<Packet>,
    /// Fission events counter
    pub fission_events: u64,
    /// Garbage collected (packets that didn't matter)
    pub garbage_collected: u64,
    /// Learning metrics - how the system grows
    pub learning: LearningMetrics,
    /// Next IDs
    next_packet_id: PacketId,
    next_node_id: NodeId,
    next_supernode_id: SupernodeId,
}

/// Metrics tracking how the system learns and grows over time
#[derive(Default, Debug, Clone, Serialize, Deserialize)]
pub struct LearningMetrics {
    /// Nodes spawned from observed patterns
    pub nodes_learned: u64,
    /// Total accesses (searches, system calls)
    pub total_accesses: u64,
    /// Successful accesses (data was useful)
    pub successful_accesses: u64,
    /// Rejections (data was wrong)
    pub rejections: u64,
    /// Boosts (data was very useful)
    pub boosts: u64,
    /// Locks (data became permanent)
    pub locks: u64,
    /// Age of the universe in ticks
    pub age: u64,
    /// Seed nodes (original, not learned)
    pub seed_nodes: u64,
}

impl Universe {
    /// Create a new universe with default config
    pub fn new() -> Self {
        Self::with_config(UniverseConfig::default())
    }

    /// Create a universe with custom config
    pub fn with_config(config: UniverseConfig) -> Self {
        Self {
            config,
            packets: HashMap::new(),
            nodes: HashMap::new(),
            supernodes: HashMap::new(),
            tick: 0,
            cmb_stats: CMBStats::default(),
            discoveries: Vec::new(),
            fission_events: 0,
            garbage_collected: 0,
            learning: LearningMetrics::default(),
            next_packet_id: 1,
            next_node_id: 1,
            next_supernode_id: 1,
        }
    }

    /// Create a universe scaled to available system memory
    pub fn from_system_memory() -> Self {
        Self::with_config(UniverseConfig::from_available_memory())
    }

    /// Spawn a new packet into the universe
    pub fn spawn_packet(&mut self, packet_type: PacketType) -> Option<PacketId> {
        if self.packets.len() >= self.config.max_packets() {
            return None; // Universe at capacity
        }

        let id = self.next_packet_id;
        self.next_packet_id += 1;

        let memory_ref = id * 1000 + self.tick;
        let mut packet = Packet::new(id, packet_type, memory_ref, self.config.radius);
        packet.born_at = self.tick;

        self.packets.insert(id, packet);
        Some(id)
    }
    
    /// Insert a pre-built packet into the universe
    pub fn insert_packet(&mut self, mut packet: Packet) -> Option<PacketId> {
        if self.packets.len() >= self.config.max_packets() {
            return None;
        }
        
        packet.born_at = self.tick;
        let id = packet.id;
        if id >= self.next_packet_id {
            self.next_packet_id = id + 1;
        }
        
        self.packets.insert(id, packet);
        Some(id)
    }

    /// Spawn a packet from actual content — THE QUANTUM TUNNELING MECHANISM
    /// 
    /// The resonance is DERIVED from the content itself. This means:
    /// - Similar content → similar resonance → clusters at same nodes
    /// - The "meaning" is encoded in the resonance signature
    /// - Sorting by resonance IS sorting by semantic similarity
    pub fn spawn_content_packet(&mut self, content: &str) -> Option<PacketId> {
        if self.packets.len() >= self.config.max_packets() {
            return None; // Universe at capacity
        }

        let id = self.next_packet_id;
        self.next_packet_id += 1;

        let memory_ref = id * 1000 + self.tick;
        let mut packet = Packet::from_content(id, content, memory_ref, self.config.radius);
        packet.born_at = self.tick;

        self.packets.insert(id, packet);
        Some(id)
    }

    /// Spawn a node at a position with given resonance
    pub fn spawn_node(&mut self, resonance: Resonance, position: Position) -> Option<NodeId> {
        if self.nodes.len() >= self.config.max_nodes() {
            return None;
        }

        let id = self.next_node_id;
        self.next_node_id += 1;

        let mut node = Node::new(id, resonance, (position.x, position.y));
        node.capture_radius = self.config.capture_radius;

        self.nodes.insert(id, node);
        Some(id)
    }

    /// Spawn seed nodes distributed around the universe
    pub fn spawn_seed_nodes(&mut self, count: usize) {
        use rand::Rng;
        let mut rng = rand::thread_rng();

        for i in 0..count {
            // Create clusters of similar resonances
            let base_wl = (i / 2) as f64 / ((count / 2) as f64).max(1.0);
            let variation = rng.gen_range(-0.05..0.05);
            let wavelength = (base_wl + variation).clamp(0.01, 0.99);
            let amplitude = rng.gen_range(0.5..1.0);

            let resonance = Resonance::full(
                wavelength,
                amplitude,
                crate::resonance::Velocity::random(),
            );

            // Position on a shell inside the universe
            let angle = (i as f64) * std::f64::consts::PI * 2.0 / (count as f64);
            let radius = self.config.radius * 0.5 + rng.gen_range(-10.0..10.0);
            let position = Position::new(
                angle.cos() * radius,
                angle.sin() * radius,
                rng.gen_range(-10.0..10.0),
            );

            self.spawn_node(resonance, position);
        }
    }

    /// Advance the universe by one tick
    /// 
    /// The tick only handles PHYSICS - movement, forces, clustering.
    /// Consumption is EVENT-DRIVEN (search, access, system calls).
    pub fn tick(&mut self) {
        self.tick += 1;

        // Phase 1: Calculate and apply forces
        self.apply_resonance_forces();

        // Phase 2: Update positions and check CMB
        self.update_positions();

        // Phase 3: Process captures (temporary - packets orbit nodes)
        self.process_captures();

        // Phase 4: Release stale captures (packets that don't fit drift away)
        if self.tick % 5 == 0 {
            self.release_poor_fits();
        }

        // Phase 5: Form supernodes
        if self.tick % 5 == 0 {
            self.form_supernodes();
        }

        // Phase 6: Learn new nodes from uncaptured packet clusters
        if self.tick % 10 == 0 {
            self.learn_new_nodes();
        }

        // Phase 7: Fission - split nodes/supernodes with too much internal tension
        if self.tick % 15 == 0 {
            self.process_fission();
        }
        
        // Phase 8: Garbage collection - clean up old drifting packets
        if self.tick % 20 == 0 {
            self.garbage_collect();
        }
    }

    /// Apply resonance-based forces to all packets
    fn apply_resonance_forces(&mut self) {
        let node_data: Vec<(NodeId, Resonance, (f64, f64), f64)> = self.nodes.values()
            .map(|n| (n.id, n.resonance.clone(), n.position, n.strength))
            .collect();

        let supernode_data: Vec<(Resonance, (f64, f64), f64)> = self.supernodes.values()
            .map(|sn| (sn.resonance.clone(), sn.center, sn.strength()))
            .collect();

        for packet in self.packets.values_mut() {
            if packet.state != PacketState::InFlight {
                continue;
            }

            // Forces from nodes
            for (_, node_res, node_pos, strength) in &node_data {
                let affinity = node_res.affinity(&packet.resonance);
                let target = Position::new(node_pos.0, node_pos.1, 0.0);
                let force = affinity * strength * self.config.force_strength;
                packet.apply_force(&target, force);
            }

            // Forces from supernodes (stronger)
            for (sn_res, sn_pos, strength) in &supernode_data {
                let affinity = sn_res.affinity(&packet.resonance);
                let target = Position::new(sn_pos.0, sn_pos.1, 0.0);
                let force = affinity * strength * self.config.force_strength * 2.0;
                packet.apply_force(&target, force);
            }

            // Apply damping
            packet.apply_damping(self.config.damping);
        }
    }

    /// Update all packet positions and check CMB boundary
    fn update_positions(&mut self) {
        let radius = self.config.radius;
        let mut cmb_packets: Vec<PacketId> = Vec::new();

        for packet in self.packets.values_mut() {
            if packet.state != PacketState::InFlight {
                continue;
            }

            packet.update_position(self.config.dt);

            if packet.check_cmb(radius) {
                cmb_packets.push(packet.id);
                self.cmb_stats.total_hits += 1;
            }
        }

        // Process CMB hits
        for packet_id in cmb_packets {
            if let Some(packet) = self.packets.get_mut(&packet_id) {
                // Anomalies that hit CMB multiple times become discoveries
                if packet.packet_type == PacketType::Anomaly && packet.cmb_hits >= 3 {
                    self.cmb_stats.discoveries += 1;
                    self.discoveries.push(packet.clone());
                    packet.state = PacketState::Released;
                    continue;
                }

                // Convert to cosmic ray after first hit (if anomaly)
                if packet.packet_type == PacketType::Anomaly && packet.cmb_hits == 1 {
                    packet.packet_type = PacketType::CosmicRay;
                    self.cmb_stats.cosmic_rays_created += 1;
                }

                // Bounce back
                packet.bounce_from_cmb(self.config.radius);
                self.cmb_stats.bounces += 1;
            }
        }
    }

    /// Process node captures
    fn process_captures(&mut self) {
        let mut captures: Vec<(NodeId, PacketId)> = Vec::new();

        for node in self.nodes.values() {
            for packet in self.packets.values() {
                if packet.state != PacketState::InFlight {
                    continue;
                }

                let dist = ((packet.position.x - node.position.0).powi(2)
                          + (packet.position.y - node.position.1).powi(2)).sqrt();

                if dist < node.capture_radius {
                    let affinity = node.resonance.affinity(&packet.resonance);
                    if affinity > 0.1 {
                        captures.push((node.id, packet.id));
                    }
                }
            }
        }

        for (node_id, packet_id) in captures {
            if let Some(node) = self.nodes.get_mut(&node_id) {
                if let Some(packet) = self.packets.get_mut(&packet_id) {
                    if packet.state == PacketState::InFlight {
                        packet.capture();
                        node.captured_packets.push(packet_id);
                        node.total_captures += 1;
                        
                        // Record wavelength for fission detection
                        node.record_wavelength(packet.resonance.wavelength);

                        // Record in supernode
                        if let Some(sn_id) = node.supernode_id {
                            if let Some(supernode) = self.supernodes.get_mut(&sn_id) {
                                supernode.record_sort();
                            }
                        }
                    }
                }
            }
        }
    }

    /// Release packets that are poor fits for their capturing node
    /// These packets will drift away and potentially find better homes
    fn release_poor_fits(&mut self) {
        let mut to_release: Vec<(NodeId, PacketId)> = Vec::new();

        for node in self.nodes.values() {
            for &packet_id in &node.captured_packets {
                if let Some(packet) = self.packets.get(&packet_id) {
                    let affinity = node.resonance.affinity(&packet.resonance);
                    // If affinity dropped below threshold, release
                    if affinity < 0.3 {
                        to_release.push((node.id, packet_id));
                    }
                }
            }
        }

        for (node_id, packet_id) in to_release {
            if let Some(node) = self.nodes.get_mut(&node_id) {
                node.release_packet(packet_id);
            }
            if let Some(packet) = self.packets.get_mut(&packet_id) {
                packet.state = PacketState::InFlight; // Back to drifting
            }
        }
    }

    /// Garbage collection - remove packets that have been drifting too long
    /// without finding a home. These are flagged as trash.
    /// 
    /// NOTE: We're now more lenient. A packet needs to:
    /// 1. Be drifting for 500+ ticks (not 100)
    /// 2. Hit the CMB 5+ times (not 2)
    /// This gives packets more time to find resonant nodes.
    fn garbage_collect(&mut self) {
        let current_tick = self.tick;
        let mut trash: Vec<PacketId> = Vec::new();

        for packet in self.packets.values() {
            // Packets drifting for 500+ ticks without capture AND 
            // hit boundary 5+ times = truly lost
            if packet.state == PacketState::InFlight 
               && current_tick - packet.born_at > 500 
               && packet.cmb_hits >= 5  // Hit boundary many times = really lost
            {
                trash.push(packet.id);
            }
        }

        let trash_count = trash.len() as u64;
        for packet_id in trash {
            if let Some(packet) = self.packets.get_mut(&packet_id) {
                packet.state = PacketState::Released; // Mark as garbage
            }
        }
        self.garbage_collected += trash_count;
    }

    // ═══════════════════════════════════════════════════════════════════
    // EVENT-DRIVEN ACCESS — Called by external searches/system calls
    // ═══════════════════════════════════════════════════════════════════

    /// Search by sending a query packet through the field
    /// 
    /// The search query becomes a PACKET that travels through the universe.
    /// It's attracted to nodes with similar resonance. When it arrives at
    /// a node, we retrieve the memory addresses of captured packets there.
    /// 
    /// This is the TRUE search — the query resonates through the field.
    pub fn search_resonant(&mut self, query: &str, max_hops: usize) -> SearchResult {
        use crate::semantic::words_to_resonance;
        use crate::packet::Position;
        
        let query_resonance = words_to_resonance(query);
        
        // Start from center
        let mut position = Position::new(0.0, 0.0, 0.0);
        let mut visited_nodes: Vec<NodeId> = Vec::new();
        let mut matches: Vec<SearchMatch> = Vec::new();
        let mut hops = 0;
        
        // The query packet travels through the field
        while hops < max_hops {
            hops += 1;
            
            // Find the most resonant node from current position
            let mut best_node: Option<(NodeId, f64, (f64, f64))> = None;
            
            for node in self.nodes.values() {
                if visited_nodes.contains(&node.id) {
                    continue;
                }
                
                let affinity = query_resonance.affinity(&node.resonance);
                if affinity > 0.2 {
                    // Calculate distance from current position
                    let dx = node.position.0 - position.x;
                    let dy = node.position.1 - position.y;
                    let dist = (dx * dx + dy * dy).sqrt();
                    
                    // Score = affinity / distance (closer & more resonant = better)
                    let score = affinity / (dist + 1.0);
                    
                    if best_node.is_none() || score > best_node.unwrap().1 {
                        best_node = Some((node.id, score, node.position));
                    }
                }
            }
            
            match best_node {
                Some((node_id, _score, node_pos)) => {
                    // Move to this node
                    position = Position::new(node_pos.0, node_pos.1, 0.0);
                    visited_nodes.push(node_id);
                    
                    // Retrieve matches from this node
                    if let Some(node) = self.nodes.get(&node_id) {
                        for &packet_id in &node.captured_packets {
                            if let Some(packet) = self.packets.get(&packet_id) {
                                // Calculate affinity from resonance
                                let resonance_affinity = query_resonance.affinity(&packet.resonance);
                                
                                // BONUS: Check content for keyword matches (semantic boost)
                                let content_bonus = if let Some(ref content) = packet.content {
                                    let query_words: Vec<&str> = query.split_whitespace()
                                        .map(|w| w.trim_matches(|c: char| !c.is_alphanumeric()))
                                        .filter(|w| w.len() > 2)
                                        .collect();
                                    let content_lower = content.to_lowercase();
                                    let matching_words = query_words.iter()
                                        .filter(|w| content_lower.contains(&w.to_lowercase()))
                                        .count();
                                    if !query_words.is_empty() {
                                        (matching_words as f64 / query_words.len() as f64) * 0.5
                                    } else {
                                        0.0
                                    }
                                } else {
                                    0.0
                                };
                                
                                let total_affinity = (resonance_affinity + content_bonus).min(1.0);
                                
                                if total_affinity > 0.3 {
                                    matches.push(SearchMatch {
                                        packet_id,
                                        memory_ref: packet.memory_ref,
                                        affinity: total_affinity,
                                        node_id,
                                        content: packet.content.clone(),
                                        found_via: MatchSource::ResonancePath,
                                    });
                                }
                            }
                        }
                    }
                    
                    // Also check nearby packets (not yet captured)
                    for packet in self.packets.values() {
                        if packet.state != PacketState::InFlight {
                            continue;
                        }
                        let dx = packet.position.x - position.x;
                        let dy = packet.position.y - position.y;
                        let dist = (dx * dx + dy * dy).sqrt();
                        
                        if dist < 10.0 {
                            let resonance_affinity = query_resonance.affinity(&packet.resonance);
                            
                            // BONUS: Check content for keyword matches
                            let content_bonus = if let Some(ref content) = packet.content {
                                let query_words: Vec<&str> = query.split_whitespace()
                                    .map(|w| w.trim_matches(|c: char| !c.is_alphanumeric()))
                                    .filter(|w| w.len() > 2)
                                    .collect();
                                let content_lower = content.to_lowercase();
                                let matching_words = query_words.iter()
                                    .filter(|w| content_lower.contains(&w.to_lowercase()))
                                    .count();
                                if !query_words.is_empty() {
                                    (matching_words as f64 / query_words.len() as f64) * 0.5
                                } else {
                                    0.0
                                }
                            } else {
                                0.0
                            };
                            
                            let total_affinity = (resonance_affinity + content_bonus).min(1.0);
                            
                            if total_affinity > 0.3 {
                                matches.push(SearchMatch {
                                    packet_id: packet.id,
                                    memory_ref: packet.memory_ref,
                                    affinity: total_affinity,
                                    node_id: 0, // Not at a node
                                    content: packet.content.clone(),
                                    found_via: MatchSource::NearbyDrift,
                                });
                            }
                        }
                    }
                }
                None => {
                    // No more resonant nodes to visit
                    break;
                }
            }
        }
        
        // FINAL PASS: CMB boundary check - scan ALL content for keyword matches
        // This is the "final scan at the end to add more context" that catches
        // what resonance alone might have missed
        for packet in self.packets.values() {
            if packet.state == PacketState::Released {
                continue;
            }
            if matches.iter().any(|m| m.packet_id == packet.id) {
                continue; // Already found
            }
            
            if let Some(ref content) = packet.content {
                let query_words: Vec<&str> = query.split_whitespace()
                    .map(|w| w.trim_matches(|c: char| !c.is_alphanumeric()))
                    .filter(|w| w.len() > 2)
                    .collect();
                let content_lower = content.to_lowercase();
                let matching_words = query_words.iter()
                    .filter(|w| content_lower.contains(&w.to_lowercase()))
                    .count();
                
                // If any keywords match, include with lower affinity (found via CMB scan)
                if matching_words > 0 {
                    let content_affinity = matching_words as f64 / query_words.len() as f64;
                    matches.push(SearchMatch {
                        packet_id: packet.id,
                        memory_ref: packet.memory_ref,
                        affinity: content_affinity * 0.8, // CMB scan penalty
                        node_id: 0,
                        content: packet.content.clone(),
                        found_via: MatchSource::CmbScan,
                    });
                }
            }
        }
        
        // Sort by affinity
        matches.sort_by(|a, b| b.affinity.partial_cmp(&a.affinity).unwrap());
        matches.dedup_by(|a, b| a.packet_id == b.packet_id);
        
        self.learning.total_accesses += 1;
        
        SearchResult {
            query: query.to_string(),
            query_resonance,
            hops,
            nodes_visited: visited_nodes,
            matches,
        }
    }

    /// Simple search (legacy) — scans all packets without traveling
    pub fn search_by_resonance(&self, query_resonance: &Resonance) -> Vec<&Packet> {
        let mut results: Vec<(&Packet, f64)> = self.packets.values()
            .filter(|p| p.state != PacketState::Released)
            .map(|p| (p, query_resonance.affinity(&p.resonance)))
            .filter(|(_, affinity)| *affinity > 0.3)
            .collect();

        results.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap());
        results.into_iter().map(|(p, _)| p).collect()
    }

    /// Search by content string — converts content to resonance, then searches
    pub fn search_by_content(&self, query: &str) -> Vec<&Packet> {
        use crate::semantic::words_to_resonance;
        let query_resonance = words_to_resonance(query);
        self.search_by_resonance(&query_resonance)
    }

    /// Access a packet — this is what happens when data is actually USED
    /// 
    /// Access can result in:
    /// - LOCK: Packet becomes permanent at its node
    /// - CONSUME: Packet is used up and released
    /// - BOOST: Packet's node gets stronger (more relevant)
    pub fn access_packet(&mut self, packet_id: PacketId, action: AccessAction) -> Option<AccessResult> {
        let _packet = self.packets.get(&packet_id)?;
        
        // Track access
        self.learning.total_accesses += 1;
        
        // Find which node (if any) has this packet
        let capturing_node_id = self.nodes.values()
            .find(|n| n.captured_packets.contains(&packet_id))
            .map(|n| n.id);

        match action {
            AccessAction::Lock => {
                // Packet becomes permanent - strengthen the node
                self.learning.locks += 1;
                self.learning.successful_accesses += 1;
                if let Some(node_id) = capturing_node_id {
                    if let Some(node) = self.nodes.get_mut(&node_id) {
                        node.strength += 0.1; // Node grows stronger
                    }
                }
                Some(AccessResult::Locked { node_id: capturing_node_id })
            }
            AccessAction::Consume => {
                // Use the data and release
                self.learning.successful_accesses += 1;
                if let Some(node_id) = capturing_node_id {
                    if let Some(node) = self.nodes.get_mut(&node_id) {
                        node.release_packet(packet_id);
                    }
                }
                if let Some(packet) = self.packets.get_mut(&packet_id) {
                    packet.consume();
                }
                Some(AccessResult::Consumed)
            }
            AccessAction::Boost => {
                // Data was useful - boost the node
                self.learning.boosts += 1;
                self.learning.successful_accesses += 1;
                if let Some(node_id) = capturing_node_id {
                    if let Some(node) = self.nodes.get_mut(&node_id) {
                        node.strength += 0.2;
                        node.capture_radius += 0.5; // Grows to capture more similar content
                    }
                }
                Some(AccessResult::Boosted { node_id: capturing_node_id })
            }
            AccessAction::Reject => {
                // Data was wrong - weaken association, let it drift
                self.learning.rejections += 1;
                if let Some(node_id) = capturing_node_id {
                    if let Some(node) = self.nodes.get_mut(&node_id) {
                        node.release_packet(packet_id);
                        node.strength = (node.strength - 0.1).max(0.5);
                    }
                }
                if let Some(packet) = self.packets.get_mut(&packet_id) {
                    packet.state = PacketState::InFlight; // Drift away
                }
                Some(AccessResult::Rejected)
            }
        }
    }

    /// Form supernodes from clustered nodes
    fn form_supernodes(&mut self) {
        use crate::node::node_affinity;

        let unassigned: Vec<NodeId> = self.nodes.values()
            .filter(|n| n.supernode_id.is_none())
            .map(|n| n.id)
            .collect();

        for node_id in unassigned {
            let mut best_match: Option<SupernodeId> = None;
            let mut best_affinity = 0.0;

            if let Some(node) = self.nodes.get(&node_id) {
                for supernode in self.supernodes.values() {
                    if supernode.should_include(node) {
                        let affinity = supernode.resonance.affinity(&node.resonance);
                        if affinity > best_affinity {
                            best_affinity = affinity;
                            best_match = Some(supernode.id);
                        }
                    }
                }
            }

            if let Some(sn_id) = best_match {
                if let Some(node) = self.nodes.get_mut(&node_id) {
                    if let Some(supernode) = self.supernodes.get_mut(&sn_id) {
                        supernode.add_member(node);
                    }
                }
            } else {
                // Try to form new supernode with another unassigned node
                let mut partner_id: Option<NodeId> = None;
                
                if let Some(node) = self.nodes.get(&node_id) {
                    for other_node in self.nodes.values() {
                        if other_node.id != node_id && other_node.supernode_id.is_none() {
                            if node_affinity(node, other_node) > 0.5 {
                                partner_id = Some(other_node.id);
                                break;
                            }
                        }
                    }
                }

                if let Some(partner) = partner_id {
                    let sn_id = self.next_supernode_id;
                    self.next_supernode_id += 1;

                    if let Some(n1) = self.nodes.get(&node_id) {
                        let mut supernode = Supernode::new(sn_id, n1);
                        
                        if let Some(n1_mut) = self.nodes.get_mut(&node_id) {
                            n1_mut.supernode_id = Some(sn_id);
                        }
                        if let Some(n2) = self.nodes.get_mut(&partner) {
                            supernode.add_member(n2);
                        }

                        self.supernodes.insert(sn_id, supernode);
                    }
                }
            }
        }

        // Update supernode centers
        let node_list: Vec<Node> = self.nodes.values().cloned().collect();
        for supernode in self.supernodes.values_mut() {
            supernode.update_center(&node_list);
        }
    }

    /// Learn new nodes from clusters of uncaptured packets
    fn learn_new_nodes(&mut self) {
        // Find packets that have been in flight for a while without being captured
        // Collect the data we need without borrowing self
        let tick = self.tick;
        let long_flight_data: Vec<(f64, f64, f64, f64)> = self.packets.values()
            .filter(|p| p.state == PacketState::InFlight && tick - p.born_at > 20)
            .map(|p| (p.position.x, p.position.y, p.resonance.wavelength, p.resonance.amplitude))
            .collect();

        if long_flight_data.len() < 3 {
            return; // Not enough uncaptured packets to learn from
        }

        // Group packets by resonance proximity (using wavelength as key similarity)
        let mut clusters: Vec<Vec<(f64, f64, f64, f64)>> = Vec::new();
        
        'outer: for packet_data in long_flight_data {
            for cluster in &mut clusters {
                if let Some(first) = cluster.first() {
                    // Simple wavelength similarity check
                    let wl_diff = (first.2 - packet_data.2).abs();
                    if wl_diff < 0.15 {
                        cluster.push(packet_data);
                        continue 'outer;
                    }
                }
            }
            clusters.push(vec![packet_data]);
        }

        // Create new node for clusters of 3+ packets
        for cluster in clusters {
            if cluster.len() >= 3 {
                // Calculate average position and resonance
                let mut avg_x = 0.0;
                let mut avg_y = 0.0;
                let mut avg_wl = 0.0;
                let mut avg_amp = 0.0;

                for (x, y, wl, amp) in &cluster {
                    avg_x += x;
                    avg_y += y;
                    avg_wl += wl;
                    avg_amp += amp;
                }

                let count = cluster.len() as f64;
                avg_x /= count;
                avg_y /= count;
                avg_wl /= count;
                avg_amp /= count;

                let resonance = Resonance::full(
                    avg_wl,
                    avg_amp,
                    crate::resonance::Velocity::random(),
                );

                let position = Position::new(avg_x, avg_y, 0.0);

                if self.spawn_node(resonance, position).is_some() {
                    // Node learned! It will start attracting these packets
                    self.learning.nodes_learned += 1;
                }
            }
        }
    }

    /// Process fission events - split nodes and supernodes with too much internal tension
    fn process_fission(&mut self) {
        self.process_node_fission();
        self.process_supernode_fission();
    }

    /// Split nodes that have captured packets with too-diverse resonances
    /// Split nodes that have processed packets with too-diverse resonances
    fn process_node_fission(&mut self) {
        use crate::node::NodeFission;
        
        // Find nodes that should undergo fission based on wavelength history
        let fission_candidates: Vec<(NodeId, NodeFission)> = self.nodes.values()
            .filter(|n| n.should_fission_from_history())
            .filter_map(|n| n.compute_fission_from_history().map(|f| (n.id, f)))
            .collect();

        for (original_node_id, fission) in fission_candidates {
            // Create new node for the lower wavelength range
            let new_node_id = self.next_node_id;
            self.next_node_id += 1;

            let new_node = Node::new(
                new_node_id,
                Resonance::new(fission.new_wavelength),
                fission.new_position,
            );

            // Update original node's resonance to the upper range
            if let Some(original) = self.nodes.get_mut(&original_node_id) {
                original.resonance = Resonance::new(fission.updated_wavelength);
                // Clear history after fission - start fresh
                original.wavelength_history.clear();
                
                // Remove from supernode if present (it may no longer belong)
                if let Some(sn_id) = original.supernode_id {
                    original.supernode_id = None;
                    if let Some(supernode) = self.supernodes.get_mut(&sn_id) {
                        supernode.member_ids.retain(|&id| id != original_node_id);
                    }
                }
            }

            self.nodes.insert(new_node_id, new_node);
            self.fission_events += 1;
        }
    }

    /// Split supernodes that have member nodes with too-diverse resonances
    fn process_supernode_fission(&mut self) {
        use crate::supernode::SupernodeFission;
        
        let node_list: Vec<Node> = self.nodes.values().cloned().collect();
        
        // Find supernodes that should undergo fission
        let fission_candidates: Vec<(SupernodeId, SupernodeFission)> = self.supernodes.values()
            .filter(|sn| sn.should_fission(&node_list))
            .filter_map(|sn| sn.compute_fission(&node_list).map(|f| (sn.id, f)))
            .collect();

        for (original_sn_id, fission) in fission_candidates {
            // Create new supernode for leaving members
            let new_sn_id = self.next_supernode_id;
            self.next_supernode_id += 1;

            // Find a seed node for the new supernode
            if let Some(&seed_id) = fission.leave_members.first() {
                if let Some(seed_node) = self.nodes.get(&seed_id) {
                    let mut new_supernode = Supernode::new(new_sn_id, seed_node);
                    new_supernode.resonance = fission.new_resonance;
                    
                    // Update nodes to point to new supernode
                    for node_id in &fission.leave_members {
                        if let Some(node) = self.nodes.get_mut(node_id) {
                            node.supernode_id = Some(new_sn_id);
                        }
                        if !new_supernode.member_ids.contains(node_id) {
                            new_supernode.member_ids.push(*node_id);
                        }
                    }

                    // Update original supernode - keep only staying members
                    if let Some(original) = self.supernodes.get_mut(&original_sn_id) {
                        original.member_ids.retain(|id| fission.stay_members.contains(id));
                        
                        // Recalculate original's resonance
                        let remaining_wl: f64 = fission.stay_members.iter()
                            .filter_map(|id| self.nodes.get(id))
                            .map(|n| n.resonance.wavelength)
                            .sum::<f64>() / fission.stay_members.len().max(1) as f64;
                        
                        original.resonance = Resonance::new(remaining_wl);
                    }

                    self.supernodes.insert(new_sn_id, new_supernode);
                    self.fission_events += 1;
                }
            }
        }
    }

    /// Render the universe as ASCII art
    pub fn render(&self, width: usize, height: usize) -> String {
        let mut canvas: Vec<Vec<char>> = vec![vec![' '; width]; height];
        let scale = self.config.radius / (width.min(height) as f64 / 2.0);

        let cx = width / 2;
        let cy = height / 2;

        // Draw CMB boundary (circle)
        let cmb_radius = (self.config.radius / scale) as usize;
        for angle in 0..360 {
            let rad = (angle as f64).to_radians();
            let x = cx as i32 + (cmb_radius as f64 * rad.cos()) as i32;
            let y = cy as i32 + (cmb_radius as f64 * rad.sin()) as i32;
            if x >= 0 && x < width as i32 && y >= 0 && y < height as i32 {
                canvas[y as usize][x as usize] = '.';
            }
        }

        // Draw supernodes (large circles)
        for supernode in self.supernodes.values() {
            let x = cx as i32 + (supernode.center.0 / scale) as i32;
            let y = cy as i32 + (supernode.center.1 / scale) as i32;
            if x >= 1 && x < (width - 1) as i32 && y >= 1 && y < (height - 1) as i32 {
                canvas[y as usize][x as usize] = 'O';
                // Draw a small ring around supernode
                for dx in -1..=1 {
                    for dy in -1..=1 {
                        if dx != 0 || dy != 0 {
                            let nx = (x + dx) as usize;
                            let ny = (y + dy) as usize;
                            if canvas[ny][nx] == ' ' {
                                canvas[ny][nx] = 'o';
                            }
                        }
                    }
                }
            }
        }

        // Draw nodes
        for node in self.nodes.values() {
            let x = cx as i32 + (node.position.0 / scale) as i32;
            let y = cy as i32 + (node.position.1 / scale) as i32;
            if x >= 0 && x < width as i32 && y >= 0 && y < height as i32 {
                let symbol = if node.supernode_id.is_some() { '+' } else { 'n' };
                canvas[y as usize][x as usize] = symbol;
            }
        }

        // Draw packets
        for packet in self.packets.values() {
            let (px, py) = packet.position.project_2d();
            let x = cx as i32 + (px / scale) as i32;
            let y = cy as i32 + (py / scale) as i32;
            if x >= 0 && x < width as i32 && y >= 0 && y < height as i32 {
                let symbol = packet.packet_type.symbol();
                canvas[y as usize][x as usize] = symbol;
            }
        }

        // Draw center marker
        canvas[cy][cx] = 'X';

        // Build string
        let mut result = String::new();
        result.push_str(&format!("╔{}╗\n", "═".repeat(width)));
        for row in &canvas {
            result.push('║');
            result.extend(row);
            result.push_str("║\n");
        }
        result.push_str(&format!("╚{}╝", "═".repeat(width)));

        result
    }

    /// Get a status report
    pub fn status(&self) -> String {
        let in_flight = self.packets.values().filter(|p| p.state == PacketState::InFlight).count();
        let captured = self.packets.values().filter(|p| p.state == PacketState::Captured).count();
        let nodes_in_sn = self.nodes.values().filter(|n| n.supernode_id.is_some()).count();

        format!(
            "Tick: {} | Packets: {} ({} flying, {} captured) | Nodes: {} ({} in supernodes) | Supernodes: {} | CMB hits: {} | Fissions: {} | Discoveries: {}",
            self.tick,
            self.packets.len(),
            in_flight,
            captured,
            self.nodes.len(),
            nodes_in_sn,
            self.supernodes.len(),
            self.cmb_stats.total_hits,
            self.fission_events,
            self.discoveries.len()
        )
    }

    /// Memory usage estimate
    pub fn memory_usage(&self) -> String {
        let packet_size = std::mem::size_of::<Packet>();
        let node_size = std::mem::size_of::<Node>();
        let used = self.packets.len() * packet_size + self.nodes.len() * node_size;
        let max = self.config.memory_bytes;
        let pct = (used as f64 / max as f64) * 100.0;

        format!(
            "Memory: {:.2} MB / {:.2} MB ({:.1}%)",
            used as f64 / (1024.0 * 1024.0),
            max as f64 / (1024.0 * 1024.0),
            pct
        )
    }

    // =========================================
    // PERSISTENCE
    // =========================================

    /// Save the universe state to a JSON file
    pub fn save<P: AsRef<Path>>(&self, path: P) -> Result<(), String> {
        let json = serde_json::to_string_pretty(self)
            .map_err(|e| format!("Serialization error: {}", e))?;
        
        fs::write(path, json)
            .map_err(|e| format!("File write error: {}", e))?;
        
        Ok(())
    }

    /// Load a universe state from a JSON file
    pub fn load<P: AsRef<Path>>(path: P) -> Result<Self, String> {
        let json = fs::read_to_string(path)
            .map_err(|e| format!("File read error: {}", e))?;
        
        serde_json::from_str(&json)
            .map_err(|e| format!("Deserialization error: {}", e))
    }

    /// Save just the learned node configurations (lightweight)
    pub fn save_nodes<P: AsRef<Path>>(&self, path: P) -> Result<(), String> {
        let nodes: Vec<&Node> = self.nodes.values().collect();
        let json = serde_json::to_string_pretty(&nodes)
            .map_err(|e| format!("Serialization error: {}", e))?;
        
        fs::write(path, json)
            .map_err(|e| format!("File write error: {}", e))?;
        
        Ok(())
    }

    /// Load and merge nodes from a file (additive learning)
    pub fn load_nodes<P: AsRef<Path>>(&mut self, path: P) -> Result<usize, String> {
        let json = fs::read_to_string(path)
            .map_err(|e| format!("File read error: {}", e))?;
        
        let nodes: Vec<Node> = serde_json::from_str(&json)
            .map_err(|e| format!("Deserialization error: {}", e))?;
        
        let mut added = 0;
        for mut node in nodes {
            // Assign new ID to avoid conflicts
            let new_id = self.next_node_id;
            self.next_node_id += 1;
            node.id = new_id;
            node.supernode_id = None; // Will be reassigned by form_supernodes
            node.captured_packets.clear();
            
            self.nodes.insert(new_id, node);
            added += 1;
        }
        
        Ok(added)
    }
}

impl Default for Universe {
    fn default() -> Self {
        Self::new()
    }
}

/// Try to detect available system memory using sysinfo
fn detect_available_memory() -> Option<usize> {
    use sysinfo::System;
    
    let sys = System::new_all();
    let available = sys.available_memory() as usize;
    
    if available > 0 {
        Some(available)
    } else {
        None
    }
}

/// Get total system memory
pub fn get_system_memory_info() -> (usize, usize) {
    use sysinfo::System;
    
    let sys = System::new_all();
    let total = sys.total_memory() as usize;
    let available = sys.available_memory() as usize;
    
    (total, available)
}
