//! Supernode — Emergent clusters of similar-resonance nodes
//!
//! Supernodes form when nodes with similar resonance congregate.
//! They act as categorical sorters, channeling packets to their
//! constituent nodes without ever inspecting the packet contents.

use crate::resonance::Resonance;
use crate::node::{Node, NodeId};
use serde::{Serialize, Deserialize};

/// Unique identifier for a supernode
pub type SupernodeId = u64;

/// An emergent cluster of similar-resonance nodes
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct Supernode {
    /// Unique identifier
    pub id: SupernodeId,
    /// The aggregate resonance of this supernode (average of members)
    pub resonance: Resonance,
    /// Member node IDs
    pub member_ids: Vec<NodeId>,
    /// Center of mass position
    pub center: (f64, f64),
    /// The resonance threshold for membership
    pub affinity_threshold: f64,
    /// Total packets sorted through this supernode
    pub packets_sorted: u64,
}

impl Supernode {
    /// Create a new supernode from a seed node
    pub fn new(id: SupernodeId, seed_node: &Node) -> Self {
        Self {
            id,
            resonance: seed_node.resonance.clone(),
            member_ids: vec![seed_node.id],
            center: seed_node.position,
            affinity_threshold: 0.6, // Nodes must have >0.6 affinity to join
            packets_sorted: 0,
        }
    }

    /// Check if a node should belong to this supernode
    pub fn should_include(&self, node: &Node) -> bool {
        self.resonance.affinity(&node.resonance) > self.affinity_threshold
    }

    /// Add a node to this supernode
    pub fn add_member(&mut self, node: &mut Node) {
        if !self.member_ids.contains(&node.id) {
            self.member_ids.push(node.id);
            node.supernode_id = Some(self.id);
            self.recalculate_resonance_with(node);
        }
    }

    /// Remove a node from this supernode
    pub fn remove_member(&mut self, node_id: NodeId) {
        self.member_ids.retain(|&id| id != node_id);
    }

    /// Update the aggregate resonance when adding a new node
    fn recalculate_resonance_with(&mut self, new_node: &Node) {
        let count = self.member_ids.len() as f64;
        if count <= 1.0 {
            self.resonance = new_node.resonance.clone();
        } else {
            // Weighted average of wavelengths
            let old_weight = (count - 1.0) / count;
            let new_weight = 1.0 / count;
            
            let new_wl = self.resonance.wavelength * old_weight 
                         + new_node.resonance.wavelength * new_weight;
            
            self.resonance = Resonance::new(new_wl);
        }
    }

    /// Update center of mass based on member positions
    pub fn update_center(&mut self, nodes: &[Node]) {
        if self.member_ids.is_empty() {
            return;
        }

        let mut sum_x = 0.0;
        let mut sum_y = 0.0;
        let mut count = 0.0;

        for node in nodes {
            if self.member_ids.contains(&node.id) {
                sum_x += node.position.0;
                sum_y += node.position.1;
                count += 1.0;
            }
        }

        if count > 0.0 {
            self.center = (sum_x / count, sum_y / count);
        }
    }

    /// Get the strength of this supernode (based on member count)
    pub fn strength(&self) -> f64 {
        (self.member_ids.len() as f64).sqrt()
    }

    /// Record that a packet was sorted through this supernode
    pub fn record_sort(&mut self) {
        self.packets_sorted += 1;
    }

    /// Check if this supernode is still viable (has enough members)
    pub fn is_viable(&self) -> bool {
        self.member_ids.len() >= 2
    }

    /// Get display signature
    pub fn signature(&self) -> String {
        format!(
            "SN{}[{}|members:{}|sorted:{}]",
            self.id,
            self.resonance.signature(),
            self.member_ids.len(),
            self.packets_sorted
        )
    }

    /// Calculate internal tension from member node resonances
    /// High tension = member nodes have very different resonances
    pub fn internal_tension(&self, nodes: &[Node]) -> f64 {
        if self.member_ids.len() < 2 {
            return 0.0;
        }

        let wavelengths: Vec<f64> = nodes.iter()
            .filter(|n| self.member_ids.contains(&n.id))
            .map(|n| n.resonance.wavelength)
            .collect();

        if wavelengths.len() < 2 {
            return 0.0;
        }

        // Calculate variance in wavelengths
        let mean: f64 = wavelengths.iter().sum::<f64>() / wavelengths.len() as f64;
        let variance: f64 = wavelengths.iter()
            .map(|wl| (wl - mean).powi(2))
            .sum::<f64>() / wavelengths.len() as f64;

        variance.sqrt()
    }

    /// Check if this supernode should undergo fission
    pub fn should_fission(&self, nodes: &[Node]) -> bool {
        // Need at least 4 members to consider splitting
        if self.member_ids.len() < 4 {
            return false;
        }
        
        let tension = self.internal_tension(nodes);
        tension > 0.12 // Lower threshold than nodes - supernodes should be more cohesive
    }

    /// Compute fission - which nodes should split off into a new supernode
    pub fn compute_fission(&self, nodes: &[Node]) -> Option<SupernodeFission> {
        if self.member_ids.len() < 4 {
            return None;
        }

        // Gather node wavelengths
        let mut node_wls: Vec<(NodeId, f64)> = nodes.iter()
            .filter(|n| self.member_ids.contains(&n.id))
            .map(|n| (n.id, n.resonance.wavelength))
            .collect();

        if node_wls.len() < 4 {
            return None;
        }

        // Sort by wavelength
        node_wls.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap());

        // Find the largest gap to split on
        let mut max_gap = 0.0;
        let mut split_idx = node_wls.len() / 2;

        for i in 1..node_wls.len() {
            let gap = node_wls[i].1 - node_wls[i-1].1;
            if gap > max_gap {
                max_gap = gap;
                split_idx = i;
            }
        }

        // Only split if gap is significant
        if max_gap < 0.08 {
            return None;
        }

        let stay_members: Vec<NodeId> = node_wls[..split_idx].iter().map(|(id, _)| *id).collect();
        let leave_members: Vec<NodeId> = node_wls[split_idx..].iter().map(|(id, _)| *id).collect();

        if stay_members.is_empty() || leave_members.is_empty() {
            return None;
        }

        // Calculate new resonance for leaving cluster
        let leave_wl: f64 = leave_members.iter()
            .filter_map(|id| nodes.iter().find(|n| n.id == *id))
            .map(|n| n.resonance.wavelength)
            .sum::<f64>() / leave_members.len() as f64;

        Some(SupernodeFission {
            stay_members,
            leave_members,
            new_resonance: Resonance::new(leave_wl),
        })
    }
}

/// Result of a supernode fission event
#[derive(Debug)]
pub struct SupernodeFission {
    /// Nodes that stay with the original supernode
    pub stay_members: Vec<NodeId>,
    /// Nodes that leave to form a new supernode
    pub leave_members: Vec<NodeId>,
    /// Suggested resonance for the new supernode
    pub new_resonance: Resonance,
}

/// Attempt to merge two supernodes if their resonances are close enough
pub fn try_merge(a: &Supernode, b: &Supernode) -> bool {
    a.resonance.affinity(&b.resonance) > 0.8
}
