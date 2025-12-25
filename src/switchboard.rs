//! Switchboard — Central registry for packet lifecycle management
//!
//! When packets are consumed or destroyed, they "snap back" to the switchboard.
//! The switchboard handles:
//! - Lookup: Packet was used, trust it will re-enter with new address
//! - Release: Consolidation candidates (50k entries for "red" → 1 canonical)
//! - Garbage collection: True waste that gets recycled

use std::collections::HashMap;
use crate::packet::{Packet, PacketId, PacketState, PacketType, MemoryRef};
use crate::resonance::Resonance;

/// What happens when a packet snaps back
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SnapbackAction {
    /// Re-enter the system with a new address
    Lookup,
    /// Consolidate with canonical reference
    Release,
    /// True garbage — recycle completely
    Garbage,
}

/// A record of a canonical reference (for deduplication)
#[derive(Clone, Debug)]
pub struct CanonicalRef {
    /// The canonical memory reference
    pub memory_ref: MemoryRef,
    /// The resonance signature of this canonical
    pub resonance: Resonance,
    /// How many packets point to this canonical
    pub reference_count: u64,
    /// The packet type this canonical represents
    pub packet_type: PacketType,
}

/// The central switchboard managing packet lifecycle
#[derive(Debug)]
pub struct Switchboard {
    /// Next available packet ID
    next_packet_id: PacketId,
    /// Next available memory reference
    next_memory_ref: MemoryRef,
    /// Packets awaiting processing (snapped back)
    pending: Vec<Packet>,
    /// Canonical references for deduplication
    canonicals: HashMap<PacketType, Vec<CanonicalRef>>,
    /// Statistics
    pub total_lookups: u64,
    pub total_releases: u64,
    pub total_garbage: u64,
    pub total_created: u64,
    /// Threshold for considering something a duplicate
    duplicate_threshold: f64,
    /// Maximum references before consolidation triggers
    consolidation_threshold: u64,
}

impl Switchboard {
    /// Create a new switchboard
    pub fn new() -> Self {
        Self {
            next_packet_id: 1,
            next_memory_ref: 1000,
            pending: Vec::new(),
            canonicals: HashMap::new(),
            total_lookups: 0,
            total_releases: 0,
            total_garbage: 0,
            total_created: 0,
            duplicate_threshold: 0.9,
            consolidation_threshold: 10, // Lower for demo purposes
        }
    }

    /// Create a new packet and assign it an ID and memory reference
    pub fn create_packet(&mut self, packet_type: PacketType, universe_radius: f64) -> Packet {
        let id = self.next_packet_id;
        self.next_packet_id += 1;
        
        let memory_ref = self.next_memory_ref;
        self.next_memory_ref += 1;
        
        self.total_created += 1;
        
        let mut packet = Packet::new(id, packet_type, memory_ref, universe_radius);
        packet.born_at = self.total_created;
        
        packet
    }

    /// Receive a packet that has snapped back
    pub fn receive_snapback(&mut self, packet: Packet) {
        self.pending.push(packet);
    }

    /// Process all pending packets
    pub fn process_pending(&mut self) -> Vec<(Packet, SnapbackAction)> {
        let packets = std::mem::take(&mut self.pending);
        let mut results = Vec::new();

        for packet in packets {
            let action = self.determine_action(&packet);
            
            match action {
                SnapbackAction::Lookup => {
                    self.total_lookups += 1;
                }
                SnapbackAction::Release => {
                    self.total_releases += 1;
                    self.record_canonical(&packet);
                }
                SnapbackAction::Garbage => {
                    self.total_garbage += 1;
                }
            }

            results.push((packet, action));
        }

        results
    }

    /// Determine what action to take for a snapped-back packet
    fn determine_action(&self, packet: &Packet) -> SnapbackAction {
        match packet.state {
            PacketState::Consumed => {
                // Was used — likely a lookup
                SnapbackAction::Lookup
            }
            PacketState::Released => {
                // Check if this is a duplicate that should be consolidated
                if self.is_duplicate(packet) {
                    SnapbackAction::Release
                } else {
                    SnapbackAction::Garbage
                }
            }
            _ => {
                // Anomaly or error state — garbage
                SnapbackAction::Garbage
            }
        }
    }

    /// Check if this packet is a duplicate of an existing canonical
    fn is_duplicate(&self, packet: &Packet) -> bool {
        if let Some(canonicals) = self.canonicals.get(&packet.packet_type) {
            for canonical in canonicals {
                if canonical.resonance.affinity(&packet.resonance) > self.duplicate_threshold {
                    return true;
                }
            }
        }
        false
    }

    /// Record a packet as a canonical reference
    fn record_canonical(&mut self, packet: &Packet) {
        let canonicals = self.canonicals
            .entry(packet.packet_type)
            .or_insert_with(Vec::new);

        // Check if we should merge with existing canonical
        for canonical in canonicals.iter_mut() {
            if canonical.resonance.affinity(&packet.resonance) > self.duplicate_threshold {
                canonical.reference_count += 1;
                return;
            }
        }

        // Create new canonical
        canonicals.push(CanonicalRef {
            memory_ref: packet.memory_ref,
            resonance: packet.resonance.clone(),
            reference_count: 1,
            packet_type: packet.packet_type,
        });
    }

    /// Check if we need to consolidate any canonicals
    pub fn needs_consolidation(&self) -> Vec<PacketType> {
        let mut result = Vec::new();
        
        for (packet_type, canonicals) in &self.canonicals {
            let total_refs: u64 = canonicals.iter()
                .map(|c| c.reference_count)
                .sum();
            
            if total_refs > self.consolidation_threshold {
                result.push(*packet_type);
            }
        }

        result
    }

    /// Perform consolidation for a packet type
    pub fn consolidate(&mut self, packet_type: PacketType) -> Option<CanonicalRef> {
        if let Some(canonicals) = self.canonicals.get_mut(&packet_type) {
            if canonicals.is_empty() {
                return None;
            }

            // Find the canonical with highest reference count
            let winner_idx = canonicals.iter()
                .enumerate()
                .max_by_key(|(_, c)| c.reference_count)
                .map(|(i, _)| i)?;

            // Merge all into the winner
            let mut winner = canonicals[winner_idx].clone();
            for (i, canonical) in canonicals.iter().enumerate() {
                if i != winner_idx {
                    winner.reference_count += canonical.reference_count;
                }
            }

            // Replace with just the winner
            *canonicals = vec![winner.clone()];

            Some(winner)
        } else {
            None
        }
    }

    /// Get statistics report
    pub fn report(&self) -> String {
        let total_processed = self.total_lookups + self.total_releases + self.total_garbage;
        let total = total_processed.max(1) as f64;
        
        format!(
            "Switchboard Stats:\n\
             - Created: {}\n\
             - Lookups: {} ({:.1}%)\n\
             - Releases: {} ({:.1}%)\n\
             - Garbage: {} ({:.1}%)\n\
             - Pending: {}",
            self.total_created,
            self.total_lookups, (self.total_lookups as f64 / total) * 100.0,
            self.total_releases, (self.total_releases as f64 / total) * 100.0,
            self.total_garbage, (self.total_garbage as f64 / total) * 100.0,
            self.pending.len()
        )
    }
}

impl Default for Switchboard {
    fn default() -> Self {
        Self::new()
    }
}
