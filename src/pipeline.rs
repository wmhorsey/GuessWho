// ═══════════════════════════════════════════════════════════════════════════
// DATA PIPELINE — Consuming real datasets into the resonance field
// ═══════════════════════════════════════════════════════════════════════════

use crate::packet::{Packet, MemoryRef};
use crate::resonance::Resonance;
use crate::universe::Universe;
use std::path::Path;
use std::fs::File;
use std::io::{BufRead, BufReader};

/// A data source that can be consumed into the universe
pub trait DataSource {
    /// Name of this data source
    fn name(&self) -> &str;
    
    /// Iterator over records
    fn records(&mut self) -> Box<dyn Iterator<Item = DataRecord> + '_>;
    
    /// Convert a record to resonance (domain-specific)
    fn to_resonance(&self, record: &DataRecord) -> Resonance;
}

/// A single record from a data source
#[derive(Debug, Clone)]
pub struct DataRecord {
    /// Unique identifier in the source
    pub id: String,
    /// Primary content (for display/preview)
    pub content: String,
    /// Structured fields
    pub fields: Vec<(String, String)>,
    /// Raw bytes (for memory reference)
    pub raw_size: usize,
}

impl DataRecord {
    pub fn new(id: impl Into<String>, content: impl Into<String>) -> Self {
        Self {
            id: id.into(),
            content: content.into(),
            fields: Vec::new(),
            raw_size: 0,
        }
    }
    
    pub fn with_field(mut self, key: impl Into<String>, value: impl Into<String>) -> Self {
        self.fields.push((key.into(), value.into()));
        self
    }
    
    pub fn with_size(mut self, size: usize) -> Self {
        self.raw_size = size;
        self
    }
}

/// Pipeline statistics
#[derive(Debug, Default)]
pub struct PipelineStats {
    pub records_processed: u64,
    pub packets_created: u64,
    pub bytes_ingested: u64,
    pub errors: u64,
}

/// The main data pipeline
pub struct Pipeline {
    /// Memory reference counter
    next_memory_ref: MemoryRef,
    /// Statistics
    pub stats: PipelineStats,
}

impl Pipeline {
    pub fn new() -> Self {
        Self {
            next_memory_ref: 1_000_000, // Start after demo range
            stats: PipelineStats::default(),
        }
    }
    
    /// Ingest a data source into the universe
    pub fn ingest<S: DataSource>(
        &mut self, 
        source: &mut S, 
        universe: &mut Universe,
        limit: Option<usize>,
    ) -> Result<u64, String> {
        let mut count = 0;
        
        // Collect records first to avoid borrow issues
        let records: Vec<DataRecord> = source.records().collect();
        let total = limit.unwrap_or(records.len()).min(records.len());
        
        println!("📥 Ingesting {} records from '{}'...", total, source.name());
        
        let start = std::time::Instant::now();
        
        for record in &records {
            if let Some(max) = limit {
                if count >= max {
                    break;
                }
            }
            
            // Convert to resonance using domain-specific logic
            let resonance = source.to_resonance(record);
            
            // Create packet
            let memory_ref = self.next_memory_ref;
            self.next_memory_ref += 1;
            
            let packet = Packet::from_content(
                self.stats.packets_created + 1,
                &record.content,
                memory_ref,
                universe.config.radius,
            ).with_resonance(resonance);
            
            universe.insert_packet(packet);
            
            self.stats.records_processed += 1;
            self.stats.packets_created += 1;
            self.stats.bytes_ingested += record.raw_size as u64;
            count += 1;
            
            // Progress indicator for large datasets
            if count % 5000 == 0 {
                let elapsed = start.elapsed().as_secs_f64();
                let rate = count as f64 / elapsed;
                print!("\r   Progress: {}/{} ({:.0} records/sec)...", count, total, rate);
            }
        }
        
        let elapsed = start.elapsed().as_secs_f64();
        println!("\r   ✓ Ingested {} records in {:.2}s ({:.0}/sec)                    ", 
                 count, elapsed, count as f64 / elapsed);
        Ok(count as u64)
    }
    
    /// Get current stats
    pub fn stats(&self) -> &PipelineStats {
        &self.stats
    }
}

// ═══════════════════════════════════════════════════════════════════════════
// CSV DATA SOURCE — Generic CSV file reader
// ═══════════════════════════════════════════════════════════════════════════

pub struct CsvSource {
    name: String,
    path: String,
    content_column: usize,
    id_column: Option<usize>,
    records: Vec<DataRecord>,
    loaded: bool,
}

impl CsvSource {
    pub fn new(name: impl Into<String>, path: impl Into<String>) -> Self {
        Self {
            name: name.into(),
            path: path.into(),
            content_column: 0,
            id_column: None,
            records: Vec::new(),
            loaded: false,
        }
    }
    
    pub fn content_column(mut self, col: usize) -> Self {
        self.content_column = col;
        self
    }
    
    pub fn id_column(mut self, col: usize) -> Self {
        self.id_column = Some(col);
        self
    }
    
    fn load(&mut self) -> Result<(), String> {
        if self.loaded {
            return Ok(());
        }
        
        let file = File::open(&self.path)
            .map_err(|e| format!("Failed to open {}: {}", self.path, e))?;
        let reader = BufReader::new(file);
        
        let mut line_num = 0;
        for line in reader.lines() {
            let line = line.map_err(|e| format!("Read error: {}", e))?;
            line_num += 1;
            
            // Skip header
            if line_num == 1 {
                continue;
            }
            
            let fields: Vec<&str> = line.split(',').collect();
            
            if fields.len() > self.content_column {
                let content = fields[self.content_column].trim().trim_matches('"');
                let id = self.id_column
                    .and_then(|i| fields.get(i))
                    .map(|s| s.trim().to_string())
                    .unwrap_or_else(|| line_num.to_string());
                
                let record = DataRecord::new(id, content)
                    .with_size(line.len());
                
                // Add all fields
                let mut record = record;
                for (i, field) in fields.iter().enumerate() {
                    record = record.with_field(format!("col_{}", i), field.trim());
                }
                
                self.records.push(record);
            }
        }
        
        self.loaded = true;
        Ok(())
    }
}

impl DataSource for CsvSource {
    fn name(&self) -> &str {
        &self.name
    }
    
    fn records(&mut self) -> Box<dyn Iterator<Item = DataRecord> + '_> {
        if let Err(e) = self.load() {
            eprintln!("Failed to load CSV: {}", e);
            return Box::new(std::iter::empty());
        }
        Box::new(self.records.iter().cloned())
    }
    
    fn to_resonance(&self, record: &DataRecord) -> Resonance {
        // Default: use word-based resonance from content
        crate::semantic::words_to_resonance(&record.content)
    }
}

// ═══════════════════════════════════════════════════════════════════════════
// CHEMISTRY DATA SOURCE — PubChem/ChEMBL format
// ═══════════════════════════════════════════════════════════════════════════

/// Chemical compound record
#[derive(Debug, Clone)]
pub struct ChemicalRecord {
    pub id: String,              // CID, ChEMBL ID, etc.
    pub name: String,            // Common name
    pub formula: String,         // Molecular formula (C6H12O6)
    pub smiles: Option<String>,  // SMILES notation
    pub inchi: Option<String>,   // InChI key
    pub mass: Option<f64>,       // Molecular weight
    pub properties: Vec<(String, String)>,
}

/// Chemistry-specific data source
pub struct ChemistrySource {
    name: String,
    records: Vec<ChemicalRecord>,
}

impl ChemistrySource {
    pub fn new(name: impl Into<String>) -> Self {
        Self {
            name: name.into(),
            records: Vec::new(),
        }
    }
    
    /// Add a chemical record
    pub fn add(&mut self, record: ChemicalRecord) {
        self.records.push(record);
    }
    
    /// Take ownership of all records (for merging sources)
    pub fn records_owned(&mut self) -> Vec<ChemicalRecord> {
        std::mem::take(&mut self.records)
    }
    
    /// Load from PubChem SDF format (simplified)
    pub fn from_pubchem_csv(path: impl AsRef<Path>) -> Result<Self, String> {
        let file = File::open(path.as_ref())
            .map_err(|e| format!("Failed to open file: {}", e))?;
        let reader = BufReader::new(file);
        
        let mut source = Self::new("PubChem");
        let mut first_line = true;
        let mut headers: Vec<String> = Vec::new();
        
        for line in reader.lines() {
            let line = line.map_err(|e| format!("Read error: {}", e))?;
            
            if first_line {
                headers = line.split(',').map(|s| s.trim().to_lowercase()).collect();
                first_line = false;
                continue;
            }
            
            let fields: Vec<&str> = line.split(',').collect();
            
            // Find relevant columns
            let get_field = |name: &str| -> Option<String> {
                headers.iter().position(|h| h.contains(name))
                    .and_then(|i| fields.get(i))
                    .map(|s| s.trim().trim_matches('"').to_string())
                    .filter(|s| !s.is_empty())
            };
            
            let record = ChemicalRecord {
                id: get_field("cid").or_else(|| get_field("id")).unwrap_or_default(),
                name: get_field("name").or_else(|| get_field("title")).unwrap_or_default(),
                formula: get_field("formula").or_else(|| get_field("molecular")).unwrap_or_default(),
                smiles: get_field("smiles"),
                inchi: get_field("inchi"),
                mass: get_field("mass").or_else(|| get_field("weight"))
                    .and_then(|s| s.parse().ok()),
                properties: Vec::new(),
            };
            
            if !record.id.is_empty() || !record.name.is_empty() {
                source.add(record);
            }
        }
        
        Ok(source)
    }
    
    /// Convert molecular formula to resonance
    /// 
    /// The key insight: molecular formula encodes STRUCTURE
    /// C6H12O6 (glucose) should resonate similarly to C6H10O5 (starch unit)
    pub fn formula_to_resonance(formula: &str) -> Resonance {
        // Parse element counts
        let mut elements: Vec<(char, u32)> = Vec::new();
        let mut current_element = String::new();
        let mut current_count = String::new();
        
        for c in formula.chars() {
            if c.is_uppercase() {
                // Save previous element
                if !current_element.is_empty() {
                    let count: u32 = current_count.parse().unwrap_or(1);
                    if let Some(first) = current_element.chars().next() {
                        elements.push((first, count));
                    }
                }
                current_element = c.to_string();
                current_count = String::new();
            } else if c.is_lowercase() {
                current_element.push(c);
            } else if c.is_numeric() {
                current_count.push(c);
            }
        }
        // Don't forget last element
        if !current_element.is_empty() {
            let count: u32 = current_count.parse().unwrap_or(1);
            if let Some(first) = current_element.chars().next() {
                elements.push((first, count));
            }
        }
        
        // Calculate resonance from element composition
        // Each element contributes based on its atomic properties
        let element_wavelength = |e: char| -> f64 {
            match e {
                'H' => 0.01,  // Hydrogen - very light
                'C' => 0.12,  // Carbon - organic backbone
                'N' => 0.14,  // Nitrogen
                'O' => 0.16,  // Oxygen
                'S' => 0.32,  // Sulfur
                'P' => 0.31,  // Phosphorus
                'F' => 0.19,  // Fluorine
                'B' => 0.11,  // Boron
                'I' => 0.53,  // Iodine
                _ => 0.50,    // Other elements
            }
        };
        
        // Weighted average wavelength
        let mut total_weight = 0.0;
        let mut weighted_lambda = 0.0;
        
        for (element, count) in &elements {
            let w = element_wavelength(*element);
            weighted_lambda += w * (*count as f64);
            total_weight += *count as f64;
        }
        
        let wavelength = if total_weight > 0.0 {
            (weighted_lambda / total_weight).min(1.0)
        } else {
            0.5
        };
        
        // Amplitude based on molecular complexity (number of atoms)
        let atom_count: u32 = elements.iter().map(|(_, c)| c).sum();
        let amplitude = (atom_count as f64 / 50.0).min(1.0);
        
        // Phase based on element diversity
        let diversity = elements.len() as f64 / 10.0;
        let phase = (diversity * std::f64::consts::PI * 2.0) % (std::f64::consts::PI * 2.0);
        
        Resonance::new(wavelength)
            .with_amplitude(amplitude)
            .with_phase(phase)
    }
    
    /// Convert SMILES to resonance (more structural info)
    pub fn smiles_to_resonance(smiles: &str) -> Resonance {
        // SMILES encodes structure: rings, bonds, branches
        let ring_count = smiles.matches(|c: char| c.is_numeric()).count();
        let branch_count = smiles.matches('(').count();
        let double_bonds = smiles.matches('=').count();
        let triple_bonds = smiles.matches('#').count();
        let aromatic = smiles.chars().filter(|c| c.is_lowercase()).count();
        
        // Wavelength from structural complexity
        let complexity = (ring_count + branch_count * 2 + aromatic) as f64;
        let wavelength = (complexity / 20.0).min(1.0);
        
        // Amplitude from size
        let size = smiles.len() as f64;
        let amplitude = (size / 100.0).min(1.0);
        
        // Phase from bond types
        let bond_factor = (double_bonds * 2 + triple_bonds * 3) as f64;
        let phase = (bond_factor * 0.5) % (std::f64::consts::PI * 2.0);
        
        Resonance::new(wavelength)
            .with_amplitude(amplitude)
            .with_phase(phase)
    }
}

impl DataSource for ChemistrySource {
    fn name(&self) -> &str {
        &self.name
    }
    
    fn records(&mut self) -> Box<dyn Iterator<Item = DataRecord> + '_> {
        Box::new(self.records.iter().map(|chem| {
            let content = if !chem.name.is_empty() {
                format!("{} ({})", chem.name, chem.formula)
            } else {
                chem.formula.clone()
            };
            
            let mut record = DataRecord::new(&chem.id, content)
                .with_field("formula", &chem.formula)
                .with_field("name", &chem.name);
            
            if let Some(ref smiles) = chem.smiles {
                record = record.with_field("smiles", smiles);
            }
            if let Some(mass) = chem.mass {
                record = record.with_field("mass", format!("{:.2}", mass));
            }
            
            record
        }))
    }
    
    fn to_resonance(&self, record: &DataRecord) -> Resonance {
        // Prefer SMILES if available, then formula
        if let Some(smiles) = record.fields.iter().find(|(k, _)| k == "smiles") {
            return Self::smiles_to_resonance(&smiles.1);
        }
        
        if let Some(formula) = record.fields.iter().find(|(k, _)| k == "formula") {
            return Self::formula_to_resonance(&formula.1);
        }
        
        // Fallback to content-based
        crate::semantic::words_to_resonance(&record.content)
    }
}

// ═══════════════════════════════════════════════════════════════════════════
// DEMO DATA — Sample chemistry data for testing
// ═══════════════════════════════════════════════════════════════════════════

impl ChemistrySource {
    /// Create demo dataset with common molecules
    pub fn demo() -> Self {
        let mut source = Self::new("Chemistry Demo");
        
        // Sugars (should cluster together)
        source.add(ChemicalRecord {
            id: "5793".into(),
            name: "Glucose".into(),
            formula: "C6H12O6".into(),
            smiles: Some("OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O".into()),
            inchi: None,
            mass: Some(180.16),
            properties: vec![("category".into(), "sugar".into())],
        });
        
        source.add(ChemicalRecord {
            id: "5988".into(),
            name: "Sucrose".into(),
            formula: "C12H22O11".into(),
            smiles: Some("OC[C@H]1O[C@@](CO)(O[C@H]2O[C@H](CO)[C@@H](O)[C@H](O)[C@H]2O)[C@@H](O)[C@@H]1O".into()),
            inchi: None,
            mass: Some(342.30),
            properties: vec![("category".into(), "sugar".into())],
        });
        
        source.add(ChemicalRecord {
            id: "6134".into(),
            name: "Fructose".into(),
            formula: "C6H12O6".into(),
            smiles: Some("OC[C@H]1O[C@](O)(CO)[C@@H](O)[C@@H]1O".into()),
            inchi: None,
            mass: Some(180.16),
            properties: vec![("category".into(), "sugar".into())],
        });
        
        // Amino acids (should cluster together)
        source.add(ChemicalRecord {
            id: "5960".into(),
            name: "Alanine".into(),
            formula: "C3H7NO2".into(),
            smiles: Some("CC(N)C(O)=O".into()),
            inchi: None,
            mass: Some(89.09),
            properties: vec![("category".into(), "amino_acid".into())],
        });
        
        source.add(ChemicalRecord {
            id: "6106".into(),
            name: "Glycine".into(),
            formula: "C2H5NO2".into(),
            smiles: Some("NCC(O)=O".into()),
            inchi: None,
            mass: Some(75.07),
            properties: vec![("category".into(), "amino_acid".into())],
        });
        
        source.add(ChemicalRecord {
            id: "6137".into(),
            name: "Leucine".into(),
            formula: "C6H13NO2".into(),
            smiles: Some("CC(C)CC(N)C(O)=O".into()),
            inchi: None,
            mass: Some(131.17),
            properties: vec![("category".into(), "amino_acid".into())],
        });
        
        source.add(ChemicalRecord {
            id: "6306".into(),
            name: "Tryptophan".into(),
            formula: "C11H12N2O2".into(),
            smiles: Some("NC(CC1=CNC2=CC=CC=C12)C(O)=O".into()),
            inchi: None,
            mass: Some(204.23),
            properties: vec![("category".into(), "amino_acid".into())],
        });
        
        // Drugs (should cluster together)
        source.add(ChemicalRecord {
            id: "2244".into(),
            name: "Aspirin".into(),
            formula: "C9H8O4".into(),
            smiles: Some("CC(=O)OC1=CC=CC=C1C(O)=O".into()),
            inchi: None,
            mass: Some(180.16),
            properties: vec![("category".into(), "drug".into())],
        });
        
        source.add(ChemicalRecord {
            id: "1983".into(),
            name: "Acetaminophen".into(),
            formula: "C8H9NO2".into(),
            smiles: Some("CC(=O)NC1=CC=C(O)C=C1".into()),
            inchi: None,
            mass: Some(151.16),
            properties: vec![("category".into(), "drug".into())],
        });
        
        source.add(ChemicalRecord {
            id: "3672".into(),
            name: "Ibuprofen".into(),
            formula: "C13H18O2".into(),
            smiles: Some("CC(C)CC1=CC=C(C=C1)C(C)C(O)=O".into()),
            inchi: None,
            mass: Some(206.28),
            properties: vec![("category".into(), "drug".into())],
        });
        
        // Organic solvents
        source.add(ChemicalRecord {
            id: "702".into(),
            name: "Ethanol".into(),
            formula: "C2H6O".into(),
            smiles: Some("CCO".into()),
            inchi: None,
            mass: Some(46.07),
            properties: vec![("category".into(), "solvent".into())],
        });
        
        source.add(ChemicalRecord {
            id: "887".into(),
            name: "Methanol".into(),
            formula: "CH4O".into(),
            smiles: Some("CO".into()),
            inchi: None,
            mass: Some(32.04),
            properties: vec![("category".into(), "solvent".into())],
        });
        
        source.add(ChemicalRecord {
            id: "180".into(),
            name: "Acetone".into(),
            formula: "C3H6O".into(),
            smiles: Some("CC(C)=O".into()),
            inchi: None,
            mass: Some(58.08),
            properties: vec![("category".into(), "solvent".into())],
        });
        
        // Hydrocarbons
        source.add(ChemicalRecord {
            id: "8030".into(),
            name: "Benzene".into(),
            formula: "C6H6".into(),
            smiles: Some("c1ccccc1".into()),
            inchi: None,
            mass: Some(78.11),
            properties: vec![("category".into(), "hydrocarbon".into())],
        });
        
        source.add(ChemicalRecord {
            id: "1140".into(),
            name: "Toluene".into(),
            formula: "C7H8".into(),
            smiles: Some("Cc1ccccc1".into()),
            inchi: None,
            mass: Some(92.14),
            properties: vec![("category".into(), "hydrocarbon".into())],
        });
        
        source.add(ChemicalRecord {
            id: "7500".into(),
            name: "Naphthalene".into(),
            formula: "C10H8".into(),
            smiles: Some("c1ccc2ccccc2c1".into()),
            inchi: None,
            mass: Some(128.17),
            properties: vec![("category".into(), "hydrocarbon".into())],
        });
        
        // Vitamins
        source.add(ChemicalRecord {
            id: "54670067".into(),
            name: "Vitamin C".into(),
            formula: "C6H8O6".into(),
            smiles: Some("OC[C@H](O)[C@H]1OC(=O)C(O)=C1O".into()),
            inchi: None,
            mass: Some(176.12),
            properties: vec![("category".into(), "vitamin".into())],
        });
        
        source.add(ChemicalRecord {
            id: "1130".into(),
            name: "Vitamin B3".into(),
            formula: "C6H5NO2".into(),
            smiles: Some("OC(=O)c1cccnc1".into()),
            inchi: None,
            mass: Some(123.11),
            properties: vec![("category".into(), "vitamin".into())],
        });
        
        // Caffeine and related
        source.add(ChemicalRecord {
            id: "2519".into(),
            name: "Caffeine".into(),
            formula: "C8H10N4O2".into(),
            smiles: Some("Cn1cnc2c1c(=O)n(c(=O)n2C)C".into()),
            inchi: None,
            mass: Some(194.19),
            properties: vec![("category".into(), "stimulant".into())],
        });
        
        source.add(ChemicalRecord {
            id: "2153".into(),
            name: "Theobromine".into(),
            formula: "C7H8N4O2".into(),
            smiles: Some("Cn1cnc2c1c(=O)[nH]c(=O)n2C".into()),
            inchi: None,
            mass: Some(180.16),
            properties: vec![("category".into(), "stimulant".into())],
        });
        
        source
    }
}

// ═══════════════════════════════════════════════════════════════════════════
// PUBCHEM FETCHER — Download from PubChem API
// ═══════════════════════════════════════════════════════════════════════════

/// Fetch compounds from PubChem (requires network)
#[cfg(feature = "network")]
pub mod pubchem {
    use super::*;
    
    const PUBCHEM_API: &str = "https://pubchem.ncbi.nlm.nih.gov/rest/pug";
    
    /// Fetch a compound by CID
    pub fn fetch_compound(cid: u64) -> Result<ChemicalRecord, String> {
        // Would use reqwest or similar
        unimplemented!("Network fetch not implemented - use demo data or CSV")
    }
    
    /// Search compounds by name
    pub fn search_compounds(query: &str, limit: usize) -> Result<Vec<ChemicalRecord>, String> {
        unimplemented!("Network fetch not implemented - use demo data or CSV")
    }
}
