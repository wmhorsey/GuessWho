//! Node — Resonant attractors in the field
//!
//! Nodes are points in the field with their own resonance signature.
//! They attract packets with similar resonance and repel those with
//! dissimilar resonance. They don't analyze packets — they just resonate.

use crate::resonance::Resonance;
use crate::packet::{Packet, PacketId, PacketState};
use serde::{Serialize, Deserialize};

/// Unique identifier for a node
pub type NodeId = u64;

/// A resonant node in the field
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct Node {
    /// Unique identifier
    pub id: NodeId,
    /// The node's resonance signature
    pub resonance: Resonance,
    /// Position in the field
    pub position: (f64, f64),
    /// Capture radius — how close a packet must be to get captured
    pub capture_radius: f64,
    /// Strength of the node's attractive/repulsive force
    pub strength: f64,
    /// Packets currently captured by this node
    pub captured_packets: Vec<PacketId>,
    /// Total packets this node has ever captured
    pub total_captures: u64,
    /// Is this node part of a supernode?
    pub supernode_id: Option<u64>,
    /// History of captured packet wavelengths (for fission detection)
    /// Tracks the last N wavelengths that passed through
    #[serde(default)]
    pub wavelength_history: Vec<f64>,
}

impl Node {
    /// Create a new node with the given resonance
    pub fn new(id: NodeId, resonance: Resonance, position: (f64, f64)) -> Self {
        Self {
            id,
            resonance,
            position,
            capture_radius: 5.0,  // Wider capture range
            strength: 3.0,        // Stronger attraction
            captured_packets: Vec::new(),
            total_captures: 0,
            supernode_id: None,
            wavelength_history: Vec::new(),
        }
    }
    
    /// Record a wavelength when a packet is captured
    pub fn record_wavelength(&mut self, wavelength: f64) {
        self.wavelength_history.push(wavelength);
        // Keep only the last 20 wavelengths
        if self.wavelength_history.len() > 20 {
            self.wavelength_history.remove(0);
        }
    }
    /// Positive = attraction, negative = repulsion
    pub fn force_on(&self, packet: &Packet) -> f64 {
        let affinity = self.resonance.affinity(&packet.resonance);
        affinity * self.strength
    }

    /// Check if a packet is within capture range
    pub fn can_capture(&self, packet: &Packet) -> bool {
        if packet.state != PacketState::InFlight {
            return false;
        }
        
        let dx = packet.position.x - self.position.0;
        let dy = packet.position.y - self.position.1;
        let distance = (dx * dx + dy * dy).sqrt();
        
        // Must be close AND have positive affinity
        distance < self.capture_radius && self.force_on(packet) > 0.1
    }

    /// Attempt to capture a packet
    pub fn capture(&mut self, packet: &mut Packet) -> bool {
        if self.can_capture(packet) {
            packet.capture();
            self.captured_packets.push(packet.id);
            self.total_captures += 1;
            true
        } else {
            false
        }
    }

    /// Release a captured packet (for consumption or garbage collection)
    pub fn release_packet(&mut self, packet_id: PacketId) -> bool {
        if let Some(pos) = self.captured_packets.iter().position(|&id| id == packet_id) {
            self.captured_packets.remove(pos);
            true
        } else {
            false
        }
    }

    /// Count of currently captured packets
    pub fn captured_count(&self) -> usize {
        self.captured_packets.len()
    }

    /// Drift the node's position slightly toward other similar nodes
    pub fn drift_toward(&mut self, target: (f64, f64), strength: f64) {
        let dx = target.0 - self.position.0;
        let dy = target.1 - self.position.1;
        self.position.0 += dx * strength * 0.1;
        self.position.1 += dy * strength * 0.1;
    }

    /// Get display signature
    pub fn signature(&self) -> String {
        format!(
            "N{}[{}|cap:{}]",
            self.id,
            self.resonance.signature(),
            self.captured_packets.len()
        )
    }
}

/// Calculate affinity between two nodes (for supernode formation)
pub fn node_affinity(a: &Node, b: &Node) -> f64 {
    a.resonance.affinity(&b.resonance)
}

/// Result of a node fission event - splits the node's resonance range
#[derive(Debug)]
pub struct NodeFission {
    /// Wavelength for the new node (covers lower range)
    pub new_wavelength: f64,
    /// Suggested position for the new node
    pub new_position: (f64, f64),
    /// Updated wavelength for original node (covers upper range)
    pub updated_wavelength: f64,
}

impl Node {
    /// Calculate internal tension from wavelength history
    /// High tension = historically very different wavelengths have passed through
    pub fn internal_tension_from_history(&self) -> f64 {
        if self.wavelength_history.len() < 4 {
            return 0.0;
        }

        // Calculate variance in wavelengths
        let mean: f64 = self.wavelength_history.iter().sum::<f64>() / self.wavelength_history.len() as f64;
        let variance: f64 = self.wavelength_history.iter()
            .map(|wl| (wl - mean).powi(2))
            .sum::<f64>() / self.wavelength_history.len() as f64;

        variance.sqrt() // Standard deviation as tension metric
    }

    /// Check if this node should undergo fission based on wavelength history
    pub fn should_fission_from_history(&self) -> bool {
        // Need enough history to detect pattern
        if self.wavelength_history.len() < 6 {
            return false;
        }
        
        let tension = self.internal_tension_from_history();
        tension > 0.12 // Threshold for fission
    }

    /// Compute fission based on wavelength history
    pub fn compute_fission_from_history(&self) -> Option<NodeFission> {
        if self.wavelength_history.len() < 6 {
            return None;
        }

        // Sort wavelengths to find the split point
        let mut sorted_wls: Vec<f64> = self.wavelength_history.clone();
        sorted_wls.sort_by(|a, b| a.partial_cmp(b).unwrap());

        // Find the largest gap to split on
        let mut max_gap = 0.0;
        let mut split_idx = sorted_wls.len() / 2;

        for i in 1..sorted_wls.len() {
            let gap = sorted_wls[i] - sorted_wls[i-1];
            if gap > max_gap {
                max_gap = gap;
                split_idx = i;
            }
        }

        // Only split if gap is significant
        if max_gap < 0.08 {
            return None;
        }

        // Calculate wavelengths for each cluster
        let lower_wls = &sorted_wls[..split_idx];
        let upper_wls = &sorted_wls[split_idx..];

        if lower_wls.is_empty() || upper_wls.is_empty() {
            return None;
        }

        let lower_wl: f64 = lower_wls.iter().sum::<f64>() / lower_wls.len() as f64;
        let upper_wl: f64 = upper_wls.iter().sum::<f64>() / upper_wls.len() as f64;

        // Offset position slightly
        let offset = 5.0;
        let new_position = (
            self.position.0 + offset,
            self.position.1 - offset,
        );

        Some(NodeFission {
            new_wavelength: lower_wl,
            new_position,
            updated_wavelength: upper_wl,
        })
    }

    /// Legacy method - use history-based fission instead
    #[allow(dead_code)]
    pub fn internal_tension(&self, packets: &std::collections::HashMap<PacketId, crate::packet::Packet>) -> f64 {
        if self.captured_packets.len() < 2 {
            return 0.0;
        }

        let wavelengths: Vec<f64> = self.captured_packets.iter()
            .filter_map(|id| packets.get(id))
            .map(|p| p.resonance.wavelength)
            .collect();

        if wavelengths.len() < 2 {
            return 0.0;
        }

        let mean: f64 = wavelengths.iter().sum::<f64>() / wavelengths.len() as f64;
        let variance: f64 = wavelengths.iter()
            .map(|wl| (wl - mean).powi(2))
            .sum::<f64>() / wavelengths.len() as f64;

        variance.sqrt()
    }

    #[allow(dead_code)]
    pub fn should_fission(&self, _packets: &std::collections::HashMap<PacketId, crate::packet::Packet>) -> bool {
        self.should_fission_from_history()
    }

    #[allow(dead_code)]
    pub fn compute_fission(&self, _packets: &std::collections::HashMap<PacketId, crate::packet::Packet>) -> Option<NodeFission> {
        self.compute_fission_from_history()
    }
}
