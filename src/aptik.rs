// aptik.rs - THE KEY: Physics Engine for Void/Spike Architecture
//
// This module contains the core physics model:
// - Quark charge stacking (+2/3, -1/3) → void/spike architecture  
// - Pinch point theory → electric charge
// - Nuclear flow predictions → chemical behavior
//
// Three sub-systems:
// 1. Quarks - fundamental particle analysis
// 2. Elements - nuclear void/spike calculation
// 3. Molecules - bond stress and flow prediction

// ============================================================================
// CORE PHYSICS: The Void/Spike Model
// ============================================================================

/// Calculate void/spike architecture for a nucleus
/// 
/// Proton (UUD) = 2 voids (intake) + 1 spike (exhaust)
/// Neutron (UDD) = 1 void (intake) + 2 spikes (exhaust)
pub fn nuclear_flow(z: usize, n: usize) -> NuclearFlow {
    let proton_voids = 2 * z;
    let proton_spikes = z;
    let neutron_voids = n;
    let neutron_spikes = 2 * n;
    
    let total_voids = proton_voids + neutron_voids;
    let total_spikes = proton_spikes + neutron_spikes;
    let flow_balance = total_voids as i32 - total_spikes as i32;
    
    NuclearFlow {
        z,
        n,
        voids: total_voids,
        spikes: total_spikes,
        flow_balance,
        ratio: total_voids as f64 / total_spikes.max(1) as f64,
    }
}

#[derive(Debug, Clone)]
pub struct NuclearFlow {
    pub z: usize,           // Proton count
    pub n: usize,           // Neutron count  
    pub voids: usize,       // Total intake capacity
    pub spikes: usize,      // Total exhaust capacity
    pub flow_balance: i32,  // Voids - Spikes (= Z - N)
    pub ratio: f64,         // Void/Spike ratio
}

impl NuclearFlow {
    /// Predict if this nuclear configuration is stable
    pub fn is_stable(&self) -> bool {
        if self.z == 0 || self.n == 0 && self.z > 1 {
            return false; // Need both nucleon types (except H-1)
        }
        
        // Valley of stability: expected ratio decreases with Z
        let expected_ratio = 1.0 - 0.0016 * (self.z as f64);
        let tolerance = 0.05;
        
        (self.ratio - expected_ratio).abs() <= tolerance
    }
    
    /// Predict chemical character from nuclear flow
    pub fn chemical_character(&self) -> ChemicalCharacter {
        if self.flow_balance > 0 {
            ChemicalCharacter::Oxidizer(self.flow_balance)
        } else if self.flow_balance < 0 {
            ChemicalCharacter::Reducer(self.flow_balance.abs())
        } else {
            ChemicalCharacter::Balanced
        }
    }
}

#[derive(Debug, Clone, PartialEq)]
pub enum ChemicalCharacter {
    Oxidizer(i32),   // Intake excess (wants electrons)
    Reducer(i32),    // Exhaust excess (gives electrons)
    Balanced,        // Equal flow
}

// ============================================================================
// QUARK MODEL: The Foundation
// ============================================================================

#[derive(Debug, Clone)]
pub struct Quark {
    pub name: &'static str,
    pub symbol: char,
    pub charge: f64,        // +2/3 or -1/3
    pub is_resonance: bool, // Strange, Charm, etc. are resonance states
    pub mass_mev: f64,
}

impl Quark {
    pub fn up() -> Self { Quark { name: "Up", symbol: 'u', charge: 2.0/3.0, is_resonance: false, mass_mev: 2.2 } }
    pub fn down() -> Self { Quark { name: "Down", symbol: 'd', charge: -1.0/3.0, is_resonance: false, mass_mev: 4.7 } }
    pub fn strange() -> Self { Quark { name: "Strange", symbol: 's', charge: -1.0/3.0, is_resonance: true, mass_mev: 95.0 } }
    pub fn charm() -> Self { Quark { name: "Charm", symbol: 'c', charge: 2.0/3.0, is_resonance: true, mass_mev: 1275.0 } }
    pub fn bottom() -> Self { Quark { name: "Bottom", symbol: 'b', charge: -1.0/3.0, is_resonance: true, mass_mev: 4180.0 } }
    pub fn top() -> Self { Quark { name: "Top", symbol: 't', charge: 2.0/3.0, is_resonance: true, mass_mev: 173000.0 } }
    
    /// Void units from this quark (intake capacity)
    pub fn voids(&self) -> usize {
        if self.charge > 0.0 { 2 } else { 0 }
    }
    
    /// Spike units from this quark (exhaust capacity)
    pub fn spikes(&self) -> usize {
        if self.charge < 0.0 { 2 } else { 0 }
    }
}

#[derive(Debug, Clone)]
pub struct Baryon {
    pub name: &'static str,
    pub quarks: [Quark; 3],
}

impl Baryon {
    pub fn proton() -> Self {
        Baryon { name: "Proton", quarks: [Quark::up(), Quark::up(), Quark::down()] }
    }
    
    pub fn neutron() -> Self {
        Baryon { name: "Neutron", quarks: [Quark::up(), Quark::down(), Quark::down()] }
    }
    
    pub fn total_voids(&self) -> usize {
        self.quarks.iter().map(|q| q.voids()).sum()
    }
    
    pub fn total_spikes(&self) -> usize {
        self.quarks.iter().map(|q| q.spikes()).sum()
    }
    
    pub fn electric_charge(&self) -> f64 {
        self.quarks.iter().map(|q| q.charge).sum()
    }
    
    /// Does this baryon have a pinch point (2+ intakes)?
    pub fn has_pinch_point(&self) -> bool {
        let up_count = self.quarks.iter().filter(|q| q.charge > 0.0).count();
        up_count >= 2
    }
    
    /// Predict stability from void/spike ratio
    pub fn stability(&self) -> BaryonStability {
        let v = self.total_voids();
        let s = self.total_spikes();
        
        if v == 0 {
            BaryonStability::Explosive  // All exhaust, no intake (Δ⁻)
        } else if s == 0 {
            BaryonStability::Implosive  // All intake, no exhaust (Δ⁺⁺)
        } else {
            let ratio = v as f64 / s as f64;
            if (ratio - 2.0).abs() < 0.1 || (ratio - 0.5).abs() < 0.1 {
                BaryonStability::Stable  // 2:1 or 1:2 ratio (p, n)
            } else if (ratio - 1.0).abs() < 0.1 {
                BaryonStability::Stable  // 1:1 balanced
            } else {
                BaryonStability::Unstable
            }
        }
    }
}

#[derive(Debug, Clone, PartialEq)]
pub enum BaryonStability {
    Stable,     // Self-sustaining flow
    Unstable,   // Will decay
    Implosive,  // All intake - collapses
    Explosive,  // All exhaust - explodes
}

// ============================================================================
// ELEMENT DATA
// ============================================================================

/// Get element info by atomic number
pub fn element_info(z: usize) -> Option<ElementInfo> {
    let data: [(usize, &str, &str, usize); 92] = [
        (1, "H", "Hydrogen", 1), (2, "He", "Helium", 4), (3, "Li", "Lithium", 7),
        (4, "Be", "Beryllium", 9), (5, "B", "Boron", 11), (6, "C", "Carbon", 12),
        (7, "N", "Nitrogen", 14), (8, "O", "Oxygen", 16), (9, "F", "Fluorine", 19),
        (10, "Ne", "Neon", 20), (11, "Na", "Sodium", 23), (12, "Mg", "Magnesium", 24),
        (13, "Al", "Aluminum", 27), (14, "Si", "Silicon", 28), (15, "P", "Phosphorus", 31),
        (16, "S", "Sulfur", 32), (17, "Cl", "Chlorine", 35), (18, "Ar", "Argon", 40),
        (19, "K", "Potassium", 39), (20, "Ca", "Calcium", 40), (21, "Sc", "Scandium", 45),
        (22, "Ti", "Titanium", 48), (23, "V", "Vanadium", 51), (24, "Cr", "Chromium", 52),
        (25, "Mn", "Manganese", 55), (26, "Fe", "Iron", 56), (27, "Co", "Cobalt", 59),
        (28, "Ni", "Nickel", 58), (29, "Cu", "Copper", 64), (30, "Zn", "Zinc", 65),
        (31, "Ga", "Gallium", 70), (32, "Ge", "Germanium", 73), (33, "As", "Arsenic", 75),
        (34, "Se", "Selenium", 79), (35, "Br", "Bromine", 80), (36, "Kr", "Krypton", 84),
        (37, "Rb", "Rubidium", 85), (38, "Sr", "Strontium", 88), (39, "Y", "Yttrium", 89),
        (40, "Zr", "Zirconium", 91), (41, "Nb", "Niobium", 93), (42, "Mo", "Molybdenum", 96),
        (43, "Tc", "Technetium", 98), (44, "Ru", "Ruthenium", 101), (45, "Rh", "Rhodium", 103),
        (46, "Pd", "Palladium", 106), (47, "Ag", "Silver", 108), (48, "Cd", "Cadmium", 112),
        (49, "In", "Indium", 115), (50, "Sn", "Tin", 119), (51, "Sb", "Antimony", 122),
        (52, "Te", "Tellurium", 128), (53, "I", "Iodine", 127), (54, "Xe", "Xenon", 131),
        (55, "Cs", "Cesium", 133), (56, "Ba", "Barium", 137), (57, "La", "Lanthanum", 139),
        (58, "Ce", "Cerium", 140), (59, "Pr", "Praseodymium", 141), (60, "Nd", "Neodymium", 144),
        (61, "Pm", "Promethium", 145), (62, "Sm", "Samarium", 150), (63, "Eu", "Europium", 152),
        (64, "Gd", "Gadolinium", 157), (65, "Tb", "Terbium", 159), (66, "Dy", "Dysprosium", 163),
        (67, "Ho", "Holmium", 165), (68, "Er", "Erbium", 167), (69, "Tm", "Thulium", 169),
        (70, "Yb", "Ytterbium", 173), (71, "Lu", "Lutetium", 175), (72, "Hf", "Hafnium", 178),
        (73, "Ta", "Tantalum", 181), (74, "W", "Tungsten", 184), (75, "Re", "Rhenium", 186),
        (76, "Os", "Osmium", 190), (77, "Ir", "Iridium", 192), (78, "Pt", "Platinum", 195),
        (79, "Au", "Gold", 197), (80, "Hg", "Mercury", 201), (81, "Tl", "Thallium", 204),
        (82, "Pb", "Lead", 207), (83, "Bi", "Bismuth", 209), (84, "Po", "Polonium", 209),
        (85, "At", "Astatine", 210), (86, "Rn", "Radon", 222), (87, "Fr", "Francium", 223),
        (88, "Ra", "Radium", 226), (89, "Ac", "Actinium", 227), (90, "Th", "Thorium", 232),
        (91, "Pa", "Protactinium", 231), (92, "U", "Uranium", 238),
    ];
    
    data.iter().find(|(elem_z, _, _, _)| *elem_z == z).map(|(_, symbol, name, mass)| {
        let n = mass - z;
        let flow = nuclear_flow(z, n);
        ElementInfo {
            z,
            n,
            symbol: symbol.to_string(),
            name: name.to_string(),
            mass: *mass,
            flow,
        }
    })
}

#[derive(Debug, Clone)]
pub struct ElementInfo {
    pub z: usize,
    pub n: usize,
    pub symbol: String,
    pub name: String,
    pub mass: usize,
    pub flow: NuclearFlow,
}

// ============================================================================
// MOLECULE ANALYSIS
// ============================================================================

/// Parse a molecular formula like "H2O" or "C6H12O6"
pub fn parse_formula(formula: &str) -> Vec<(String, usize)> {
    let mut atoms = Vec::new();
    let mut current_element = String::new();
    let mut current_count = String::new();
    
    for c in formula.chars() {
        if c.is_uppercase() {
            if !current_element.is_empty() {
                let count = if current_count.is_empty() { 1 } 
                           else { current_count.parse().unwrap_or(1) };
                atoms.push((current_element, count));
                current_count.clear();
            }
            current_element = c.to_string();
        } else if c.is_lowercase() {
            current_element.push(c);
        } else if c.is_numeric() {
            current_count.push(c);
        }
    }
    if !current_element.is_empty() {
        let count = if current_count.is_empty() { 1 } 
                   else { current_count.parse().unwrap_or(1) };
        atoms.push((current_element, count));
    }
    atoms
}

/// Analyze a molecule's void/spike architecture
pub fn analyze_molecule_flow(formula: &str) -> Option<MoleculeFlow> {
    let atoms = parse_formula(formula);
    
    let mut total_voids = 0usize;
    let mut total_spikes = 0usize;
    let mut atom_flows: Vec<AtomContribution> = Vec::new();
    
    for (symbol, count) in &atoms {
        // Look up element by symbol
        let elem = (1..=92).filter_map(element_info)
            .find(|e| e.symbol == *symbol)?;
        
        let atom_voids = elem.flow.voids * count;
        let atom_spikes = elem.flow.spikes * count;
        
        total_voids += atom_voids;
        total_spikes += atom_spikes;
        
        atom_flows.push(AtomContribution {
            symbol: symbol.clone(),
            count: *count,
            z: elem.z,
            voids_per_atom: elem.flow.voids,
            spikes_per_atom: elem.flow.spikes,
            character: atom_chemical_character(elem.z),
        });
    }
    
    let flow_balance = total_voids as i32 - total_spikes as i32;
    
    Some(MoleculeFlow {
        formula: formula.to_string(),
        atoms: atom_flows,
        total_voids,
        total_spikes,
        flow_balance,
        ratio: total_voids as f64 / total_spikes.max(1) as f64,
    })
}

/// Get chemical character based on Z (periodic table position)
pub fn atom_chemical_character(z: usize) -> AtomCharacter {
    match z {
        2 | 10 | 18 | 36 | 54 | 86 => AtomCharacter::Inert,
        8 | 9 | 17 | 35 => AtomCharacter::StrongOxidizer,
        7 | 16 | 34 => AtomCharacter::Oxidizer,
        6 | 15 | 33 | 5 => AtomCharacter::MildOxidizer,
        1 | 3 | 11 | 19 | 37 | 55 => AtomCharacter::StrongReducer,
        4 | 12 | 20 | 38 | 56 => AtomCharacter::Reducer,
        21..=30 | 39..=48 | 72..=80 => AtomCharacter::Transition,
        _ => AtomCharacter::MildReducer,
    }
}

#[derive(Debug, Clone)]
pub struct MoleculeFlow {
    pub formula: String,
    pub atoms: Vec<AtomContribution>,
    pub total_voids: usize,
    pub total_spikes: usize,
    pub flow_balance: i32,
    pub ratio: f64,
}

#[derive(Debug, Clone)]
pub struct AtomContribution {
    pub symbol: String,
    pub count: usize,
    pub z: usize,
    pub voids_per_atom: usize,
    pub spikes_per_atom: usize,
    pub character: AtomCharacter,
}

#[derive(Debug, Clone, PartialEq)]
pub enum AtomCharacter {
    StrongOxidizer,
    Oxidizer,
    MildOxidizer,
    Inert,
    MildReducer,
    Reducer,
    StrongReducer,
    Transition,
}

impl MoleculeFlow {
    /// Predict overall molecule character
    pub fn predict_character(&self) -> ChemicalCharacter {
        if self.flow_balance > 0 {
            ChemicalCharacter::Oxidizer(self.flow_balance)
        } else if self.flow_balance < 0 {
            ChemicalCharacter::Reducer(self.flow_balance.abs())
        } else {
            ChemicalCharacter::Balanced
        }
    }
    
    /// Check for bond stress (same-element bonds between oxidizers)
    pub fn bond_stress_warnings(&self) -> Vec<BondStress> {
        let mut warnings = Vec::new();
        
        for atom in &self.atoms {
            if atom.count >= 2 {
                match atom.character {
                    AtomCharacter::StrongOxidizer | AtomCharacter::Oxidizer => {
                        warnings.push(BondStress {
                            symbol: atom.symbol.clone(),
                            stress_type: StressType::OxidizerOxidizer,
                            severity: if atom.character == AtomCharacter::StrongOxidizer { 
                                "HIGH" 
                            } else { 
                                "MODERATE" 
                            }.to_string(),
                        });
                    },
                    AtomCharacter::StrongReducer | AtomCharacter::Reducer => {
                        warnings.push(BondStress {
                            symbol: atom.symbol.clone(),
                            stress_type: StressType::ReducerReducer,
                            severity: "LOW".to_string(),
                        });
                    },
                    _ => {}
                }
            }
        }
        
        warnings
    }
}

#[derive(Debug, Clone)]
pub struct BondStress {
    pub symbol: String,
    pub stress_type: StressType,
    pub severity: String,
}

#[derive(Debug, Clone, PartialEq)]
pub enum StressType {
    OxidizerOxidizer,  // Tug of war (like O-O in H2O2)
    ReducerReducer,    // Shared exhaust (like H-H)
}

// ============================================================================
// ISOTOPE STABILITY
// ============================================================================

/// Predict isotope stability using the valley of stability
pub fn predict_isotope_stability(z: usize, mass: usize) -> IsotopeStability {
    let n = mass.saturating_sub(z);
    let flow = nuclear_flow(z, n);
    
    // Valley of stability curve
    let expected_ratio = 1.0 - 0.0016 * (z as f64);
    let tolerance = 0.05;
    let deviation = (flow.ratio - expected_ratio).abs();
    
    let status = if mass == 1 && n == 0 {
        StabilityStatus::EdgeCase  // Protium - just a proton
    } else if z > 83 {
        StabilityStatus::BeyondCliff  // Above bismuth - all radioactive
    } else if deviation <= tolerance {
        StabilityStatus::Stable
    } else if deviation <= tolerance * 2.0 {
        StabilityStatus::Marginal
    } else {
        StabilityStatus::Unstable
    };
    
    IsotopeStability {
        z,
        n,
        mass,
        flow,
        expected_ratio,
        deviation,
        status,
    }
}

#[derive(Debug, Clone)]
pub struct IsotopeStability {
    pub z: usize,
    pub n: usize,
    pub mass: usize,
    pub flow: NuclearFlow,
    pub expected_ratio: f64,
    pub deviation: f64,
    pub status: StabilityStatus,
}

#[derive(Debug, Clone, PartialEq)]
pub enum StabilityStatus {
    Stable,
    Marginal,
    Unstable,
    EdgeCase,    // Single nucleon
    BeyondCliff, // Z > 83
}

// ============================================================================
// TESTS
// ============================================================================

#[cfg(test)]
mod tests {
    use super::*;
    
    #[test]
    fn test_proton_flow() {
        let p = Baryon::proton();
        assert_eq!(p.total_voids(), 4);
        assert_eq!(p.total_spikes(), 2);
        assert!(p.has_pinch_point());
        assert_eq!(p.stability(), BaryonStability::Stable);
    }
    
    #[test]
    fn test_neutron_flow() {
        let n = Baryon::neutron();
        assert_eq!(n.total_voids(), 2);
        assert_eq!(n.total_spikes(), 4);
        assert!(!n.has_pinch_point());
        assert_eq!(n.stability(), BaryonStability::Stable);
    }
    
    #[test]
    fn test_h2o2_stress() {
        let mol = analyze_molecule_flow("H2O2").unwrap();
        let stress = mol.bond_stress_warnings();
        assert!(stress.iter().any(|s| s.symbol == "O" && s.stress_type == StressType::OxidizerOxidizer));
    }
    
    #[test]
    fn test_carbon_balanced() {
        let elem = element_info(6).unwrap();
        assert_eq!(elem.flow.ratio, 1.0);
    }
}
