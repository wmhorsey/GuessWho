// ═══════════════════════════════════════════════════════════════════════════
// DATA LOADERS — Read local chemistry data files
// ═══════════════════════════════════════════════════════════════════════════
//
// Supported formats:
// - CSV (PubChem exports, custom)
// - SDF (Structure Data File - standard chemistry format)
// - SMILES (simple line format)
// - TSV (tab-separated)

use crate::pipeline::{ChemicalRecord, ChemistrySource, DataRecord, DataSource};
use crate::resonance::Resonance;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

// ═══════════════════════════════════════════════════════════════════════════
// CSV LOADER — Most common format from PubChem downloads
// ═══════════════════════════════════════════════════════════════════════════

/// Load chemistry data from a CSV file
/// 
/// Expected columns (flexible - will auto-detect):
/// - CID or id: compound identifier
/// - name or title: common name
/// - MolecularFormula or formula: e.g., C6H12O6
/// - CanonicalSMILES or smiles: SMILES notation
/// - MolecularWeight or mass: molecular weight
/// 
/// Example PubChem download URL (manual download):
/// https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/2244,1983,3672/CSV
pub fn load_csv(path: impl AsRef<Path>) -> Result<ChemistrySource, String> {
    let path = path.as_ref();
    let file = File::open(path)
        .map_err(|e| format!("Failed to open '{}': {}", path.display(), e))?;
    let reader = BufReader::new(file);
    
    let name = path.file_stem()
        .and_then(|s| s.to_str())
        .unwrap_or("CSV Data")
        .to_string();
    
    let mut source = ChemistrySource::new(&name);
    let mut lines = reader.lines();
    
    // Parse header to find column indices
    let header_line = lines.next()
        .ok_or("Empty CSV file")?
        .map_err(|e| format!("Failed to read header: {}", e))?;
    
    let headers: Vec<String> = parse_csv_line(&header_line)
        .iter()
        .map(|s| s.to_lowercase())
        .collect();
    
    let col_id = find_column(&headers, &["cid", "id", "compound_id", "pubchem_cid"]);
    let col_name = find_column(&headers, &["name", "title", "iupacname", "compound_name"]);
    let col_formula = find_column(&headers, &["molecularformula", "formula", "molecular_formula"]);
    let col_smiles = find_column(&headers, &["canonicalsmiles", "smiles", "isomericsmiles"]);
    let col_mass = find_column(&headers, &["molecularweight", "mass", "molecular_weight", "mw"]);
    
    let mut count = 0;
    let mut line_num = 1;
    
    for line_result in lines {
        line_num += 1;
        let line = match line_result {
            Ok(l) => l,
            Err(_) => continue,
        };
        
        if line.trim().is_empty() {
            continue;
        }
        
        let fields = parse_csv_line(&line);
        
        let record = ChemicalRecord {
            id: col_id.and_then(|i| fields.get(i)).cloned().unwrap_or_else(|| line_num.to_string()),
            name: col_name.and_then(|i| fields.get(i)).cloned().unwrap_or_default(),
            formula: col_formula.and_then(|i| fields.get(i)).cloned().unwrap_or_default(),
            smiles: col_smiles.and_then(|i| fields.get(i)).cloned().filter(|s| !s.is_empty()),
            inchi: None,
            mass: col_mass.and_then(|i| fields.get(i)).and_then(|s| s.parse().ok()),
            properties: Vec::new(),
        };
        
        // Only add if we have meaningful data
        if !record.id.is_empty() && (!record.name.is_empty() || !record.formula.is_empty()) {
            source.add(record);
            count += 1;
        }
    }
    
    // Only print if loaded a significant amount
    if count >= 10 {
        println!("   Loaded {} compounds from {}", count, path.display());
    }
    Ok(source)
}

/// Parse a CSV line handling quoted fields
fn parse_csv_line(line: &str) -> Vec<String> {
    let mut fields = Vec::new();
    let mut current = String::new();
    let mut in_quotes = false;
    
    for c in line.chars() {
        match c {
            '"' => in_quotes = !in_quotes,
            ',' if !in_quotes => {
                fields.push(current.trim().to_string());
                current = String::new();
            }
            _ => current.push(c),
        }
    }
    fields.push(current.trim().to_string());
    fields
}

/// Find column index by trying multiple possible names
fn find_column(headers: &[String], names: &[&str]) -> Option<usize> {
    for name in names {
        if let Some(i) = headers.iter().position(|h| h.contains(name)) {
            return Some(i);
        }
    }
    None
}

// ═══════════════════════════════════════════════════════════════════════════
// SDF LOADER — Structure Data File format (standard chemistry)
// ═══════════════════════════════════════════════════════════════════════════

/// Load chemistry data from an SDF file
/// 
/// SDF is the standard format for chemical structures.
/// Each record contains:
/// - Molecule name (first line)
/// - Atom/bond block (MOL format)
/// - Properties as key-value pairs after "M  END"
/// 
/// Download from PubChem:
/// https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/2244/SDF
pub fn load_sdf(path: impl AsRef<Path>) -> Result<ChemistrySource, String> {
    let path = path.as_ref();
    let file = File::open(path)
        .map_err(|e| format!("Failed to open '{}': {}", path.display(), e))?;
    let reader = BufReader::new(file);
    
    let name = path.file_stem()
        .and_then(|s| s.to_str())
        .unwrap_or("SDF Data")
        .to_string();
    
    let mut source = ChemistrySource::new(&name);
    
    let mut current_name = String::new();
    let mut properties: Vec<(String, String)> = Vec::new();
    let mut current_property_name = String::new();
    let mut in_mol_block = false;
    let mut count = 0;
    
    for line_result in reader.lines() {
        let line = match line_result {
            Ok(l) => l,
            Err(_) => continue,
        };
        
        // Record separator
        if line.starts_with("$$$$") {
            // Save current record
            if !current_name.is_empty() || !properties.is_empty() {
                let mut record = ChemicalRecord {
                    id: count.to_string(),
                    name: current_name.clone(),
                    formula: String::new(),
                    smiles: None,
                    inchi: None,
                    mass: None,
                    properties: Vec::new(),
                };
                
                // Extract standard properties
                for (key, value) in &properties {
                    let key_lower = key.to_lowercase();
                    if key_lower.contains("formula") {
                        record.formula = value.clone();
                    } else if key_lower.contains("smiles") {
                        record.smiles = Some(value.clone());
                    } else if key_lower.contains("inchi") && !key_lower.contains("key") {
                        record.inchi = Some(value.clone());
                    } else if key_lower.contains("weight") || key_lower.contains("mass") {
                        record.mass = value.parse().ok();
                    } else if key_lower.contains("cid") || key_lower.contains("id") {
                        record.id = value.clone();
                    }
                }
                
                source.add(record);
                count += 1;
            }
            
            // Reset for next record
            current_name = String::new();
            properties = Vec::new();
            current_property_name = String::new();
            in_mol_block = false;
            continue;
        }
        
        // First line of a new record is the molecule name
        if !in_mol_block && current_name.is_empty() && !line.trim().is_empty() {
            current_name = line.trim().to_string();
            in_mol_block = true;
            continue;
        }
        
        // Property header
        if line.starts_with("> <") {
            // Extract property name from "> <PROPERTY_NAME>"
            if let Some(end) = line.rfind('>') {
                current_property_name = line[3..end].to_string();
            }
            continue;
        }
        
        // Property value (line after property header)
        if !current_property_name.is_empty() && !line.trim().is_empty() {
            properties.push((current_property_name.clone(), line.trim().to_string()));
            current_property_name = String::new();
        }
    }
    
    println!("   Loaded {} compounds from {}", count, path.display());
    Ok(source)
}

// ═══════════════════════════════════════════════════════════════════════════
// SMILES LOADER — Simple line format (one molecule per line)
// ═══════════════════════════════════════════════════════════════════════════

/// Load from a SMILES file
/// 
/// Format: SMILES<tab or space>NAME (one per line)
/// Example:
/// CCO    Ethanol
/// CC(C)=O    Acetone
pub fn load_smiles(path: impl AsRef<Path>) -> Result<ChemistrySource, String> {
    let path = path.as_ref();
    let file = File::open(path)
        .map_err(|e| format!("Failed to open '{}': {}", path.display(), e))?;
    let reader = BufReader::new(file);
    
    let name = path.file_stem()
        .and_then(|s| s.to_str())
        .unwrap_or("SMILES Data")
        .to_string();
    
    let mut source = ChemistrySource::new(&name);
    let mut count = 0;
    
    for line_result in reader.lines() {
        let line = match line_result {
            Ok(l) => l,
            Err(_) => continue,
        };
        
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        
        // Split on tab or multiple spaces
        let parts: Vec<&str> = line.split(|c| c == '\t' || c == ' ')
            .filter(|s| !s.is_empty())
            .collect();
        
        if parts.is_empty() {
            continue;
        }
        
        let smiles = parts[0].to_string();
        let mol_name = if parts.len() > 1 {
            parts[1..].join(" ")
        } else {
            format!("Compound_{}", count + 1)
        };
        
        let record = ChemicalRecord {
            id: (count + 1).to_string(),
            name: mol_name,
            formula: String::new(), // Could compute from SMILES
            smiles: Some(smiles),
            inchi: None,
            mass: None,
            properties: Vec::new(),
        };
        
        source.add(record);
        count += 1;
    }
    
    println!("   Loaded {} compounds from {}", count, path.display());
    Ok(source)
}

// ═══════════════════════════════════════════════════════════════════════════
// AUTO LOADER — Detect format from extension
// ═══════════════════════════════════════════════════════════════════════════

/// Automatically load based on file extension
pub fn load_auto(path: impl AsRef<Path>) -> Result<ChemistrySource, String> {
    let path = path.as_ref();
    let ext = path.extension()
        .and_then(|s| s.to_str())
        .map(|s| s.to_lowercase())
        .unwrap_or_default();
    
    match ext.as_str() {
        "csv" => load_csv(path),
        "tsv" => load_csv(path), // TSV is handled by CSV parser
        "sdf" | "mol" => load_sdf(path),
        "smi" | "smiles" => load_smiles(path),
        _ => Err(format!("Unknown file format: .{}", ext)),
    }
}

// ═══════════════════════════════════════════════════════════════════════════
// DIRECTORY LOADER — Load all chemistry files from a directory (recursive)
// ═══════════════════════════════════════════════════════════════════════════

/// Load all chemistry files from a directory (recursive)
pub fn load_directory(dir: impl AsRef<Path>) -> Result<Vec<ChemistrySource>, String> {
    let dir = dir.as_ref();
    
    if !dir.is_dir() {
        return Err(format!("'{}' is not a directory", dir.display()));
    }
    
    let mut sources = Vec::new();
    load_directory_recursive(dir, &mut sources)?;
    Ok(sources)
}

fn load_directory_recursive(dir: &Path, sources: &mut Vec<ChemistrySource>) -> Result<(), String> {
    let entries = std::fs::read_dir(dir)
        .map_err(|e| format!("Failed to read directory: {}", e))?;
    
    for entry in entries {
        let entry = match entry {
            Ok(e) => e,
            Err(_) => continue,
        };
        
        let path = entry.path();
        
        // Recurse into subdirectories
        if path.is_dir() {
            load_directory_recursive(&path, sources)?;
            continue;
        }
        
        if !path.is_file() {
            continue;
        }
        
        let ext = path.extension()
            .and_then(|s| s.to_str())
            .map(|s| s.to_lowercase())
            .unwrap_or_default();
        
        if matches!(ext.as_str(), "csv" | "tsv" | "sdf" | "mol" | "smi" | "smiles") {
            match load_auto(&path) {
                Ok(source) => sources.push(source),
                Err(e) => eprintln!("   Warning: {}", e),
            }
        }
    }
    
    Ok(())
}

// ═══════════════════════════════════════════════════════════════════════════
// SAMPLE DATA GENERATOR — Create test files
// ═══════════════════════════════════════════════════════════════════════════

/// Generate sample CSV file for testing
pub fn generate_sample_csv(path: impl AsRef<Path>) -> Result<(), String> {
    use std::io::Write;
    
    let path = path.as_ref();
    let mut file = File::create(path)
        .map_err(|e| format!("Failed to create file: {}", e))?;
    
    writeln!(file, "CID,Name,MolecularFormula,CanonicalSMILES,MolecularWeight").unwrap();
    
    // Common compounds
    let compounds = [
        ("5793", "Glucose", "C6H12O6", "OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O", "180.16"),
        ("5988", "Sucrose", "C12H22O11", "OC[C@H]1O[C@@](CO)(O[C@H]2O[C@H](CO)[C@@H](O)[C@H](O)[C@H]2O)[C@@H](O)[C@@H]1O", "342.30"),
        ("6134", "Fructose", "C6H12O6", "OC[C@H]1O[C@](O)(CO)[C@@H](O)[C@@H]1O", "180.16"),
        ("2244", "Aspirin", "C9H8O4", "CC(=O)OC1=CC=CC=C1C(O)=O", "180.16"),
        ("1983", "Acetaminophen", "C8H9NO2", "CC(=O)NC1=CC=C(O)C=C1", "151.16"),
        ("3672", "Ibuprofen", "C13H18O2", "CC(C)CC1=CC=C(C=C1)C(C)C(O)=O", "206.28"),
        ("2519", "Caffeine", "C8H10N4O2", "Cn1cnc2c1c(=O)n(c(=O)n2C)C", "194.19"),
        ("5960", "Alanine", "C3H7NO2", "CC(N)C(O)=O", "89.09"),
        ("6106", "Glycine", "C2H5NO2", "NCC(O)=O", "75.07"),
        ("702", "Ethanol", "C2H6O", "CCO", "46.07"),
        ("887", "Methanol", "CH4O", "CO", "32.04"),
        ("180", "Acetone", "C3H6O", "CC(C)=O", "58.08"),
        ("8030", "Benzene", "C6H6", "c1ccccc1", "78.11"),
        ("1140", "Toluene", "C7H8", "Cc1ccccc1", "92.14"),
        ("6137", "Leucine", "C6H13NO2", "CC(C)CC(N)C(O)=O", "131.17"),
        ("6306", "Tryptophan", "C11H12N2O2", "NC(CC1=CNC2=CC=CC=C12)C(O)=O", "204.23"),
        ("54670067", "Vitamin C", "C6H8O6", "OC[C@H](O)[C@H]1OC(=O)C(O)=C1O", "176.12"),
        ("962", "Water", "H2O", "O", "18.02"),
        ("6324", "Acetic Acid", "C2H4O2", "CC(=O)O", "60.05"),
        ("7500", "Naphthalene", "C10H8", "c1ccc2ccccc2c1", "128.17"),
    ];
    
    for (cid, name, formula, smiles, mass) in compounds {
        writeln!(file, "{},\"{}\",{},{},{}", cid, name, formula, smiles, mass).unwrap();
    }
    
    println!("✓ Generated sample CSV: {}", path.display());
    Ok(())
}

/// Generate sample SMILES file for testing
pub fn generate_sample_smiles(path: impl AsRef<Path>) -> Result<(), String> {
    use std::io::Write;
    
    let path = path.as_ref();
    let mut file = File::create(path)
        .map_err(|e| format!("Failed to create file: {}", e))?;
    
    writeln!(file, "# Common organic compounds").unwrap();
    writeln!(file, "CCO\tEthanol").unwrap();
    writeln!(file, "CC(C)=O\tAcetone").unwrap();
    writeln!(file, "c1ccccc1\tBenzene").unwrap();
    writeln!(file, "Cc1ccccc1\tToluene").unwrap();
    writeln!(file, "CC(=O)O\tAcetic acid").unwrap();
    writeln!(file, "CO\tMethanol").unwrap();
    writeln!(file, "CC(=O)OC1=CC=CC=C1C(O)=O\tAspirin").unwrap();
    writeln!(file, "CC(=O)NC1=CC=C(O)C=C1\tAcetaminophen").unwrap();
    writeln!(file, "NCC(O)=O\tGlycine").unwrap();
    writeln!(file, "CC(N)C(O)=O\tAlanine").unwrap();
    
    println!("✓ Generated sample SMILES: {}", path.display());
    Ok(())
}
