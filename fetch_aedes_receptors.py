import os
import requests
import pandas as pd
from Bio import PDB
from Bio.SeqUtils import ProtParam
import pubchempy as pcp
from rdkit import Chem
from rdkit.Chem import Descriptors, Draw

# ==========================================
# CONFIGURATION & DIRECTORY SETUP
# ==========================================
DATA_DIR = "./aedes_data"
PDB_DIR = os.path.join(DATA_DIR, "structures")
LIGAND_DIR = os.path.join(DATA_DIR, "ligands")

os.makedirs(PDB_DIR, exist_ok=True)
os.makedirs(LIGAND_DIR, exist_ok=True)

# Correct Canonical UniProt IDs for Aedes aegypti Receptors
TARGET_RECEPTORS = {
    "AaegOr4": "A0A1S7UE45",   # Human-preference Odorant Receptor
    "AaegOrco": "Q177B5"      # Universal Odorant Co-Receptor
}

TARGET_LIGANDS = [
    "Sulcatone",       # 6-methyl-5-hepten-2-one (Key human skin odorant)
    "1-Octen-3-ol",    # Mosquito attractant compound
    "Lactic acid",     # Human sweat volatile component
    "DEET"             # Standard synthetic repellent control
]

# ==========================================
# 1. ROBUST ALPHAFOLD STRUCTURE RETRIEVAL
# ==========================================
def fetch_alphafold_structure(uniprot_id: str, target_name: str) -> str:
    """Retrieves 3D protein structures from AlphaFold DB using its REST API."""
    api_url = f"https://alphafold.ebi.ac.uk/api/prediction/{uniprot_id}"
    print(f"[*] Querying AlphaFold API for {target_name} ({uniprot_id})...")
    
    try:
        response = requests.get(api_url)
        if response.status_code == 200:
            data = response.json()
            if data and len(data) > 0:
                # Get the PDB download URL from AlphaFold API metadata
                pdb_url = data[0].get("pdbUrl")
                if pdb_url:
                    pdb_res = requests.get(pdb_url)
                    output_path = os.path.join(PDB_DIR, f"{target_name}_{uniprot_id}.pdb")
                    with open(output_path, "wb") as f:
                        f.write(pdb_res.content)
                    print(f"[+] Saved PDB structure to {output_path}")
                    return output_path
        
        # Fallback direct URL format if API metadata isn't returned
        fallback_url = f"https://alphafold.ebi.ac.uk/files/AF-{uniprot_id}-F1-model_v4.pdb"
        fb_res = requests.get(fallback_url)
        if fb_res.status_code == 200:
            output_path = os.path.join(PDB_DIR, f"{target_name}_{uniprot_id}.pdb")
            with open(output_path, "wb") as f:
                f.write(fb_res.content)
            print(f"[+] Saved PDB (via fallback) to {output_path}")
            return output_path
            
    except Exception as e:
        print(f"[-] Error fetching structure for {uniprot_id}: {e}")
        
    print(f"[-] Could not retrieve PDB for {uniprot_id}.")
    return None

def analyze_protein_structure(pdb_path: str):
    """Parses sequence and structural metrics from a PDB file using Biopython."""
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure("receptor", pdb_path)
    
    ppb = PDB.PPBuilder()
    peptides = ppb.build_peptides(structure)
    
    if not peptides:
        return {}
    
    sequence = str(peptides[0].get_sequence())
    analyser = ProtParam.ProteinAnalysis(sequence)
    
    return {
        "Residue_Count": len(sequence),
        "Molecular_Weight_Da": round(analyser.molecular_weight(), 2),
        "Isoelectric_Point": round(analyser.isoelectric_point(), 2),
        "Aromaticity": round(analyser.aromaticity(), 4)
    }

# ==========================================
# 2. FETCH LIGANDS & COMPUTE RDKIT METRICS
# ==========================================
def fetch_and_analyze_ligands(ligand_names: list):
    """Retrieves compounds from PubChem and calculates 2D descriptors using RDKit."""
    ligand_data = []
    rdkit_mols = []
    
    print("\n[*] Processing Volatile Ligands via PubChem & RDKit...")
    for name in ligand_names:
        compounds = pcp.get_compounds(name, "name")
        if not compounds:
            continue
            
        cmpd = compounds[0]
        smiles = getattr(cmpd, "smiles", getattr(cmpd, "isomeric_smiles", None))
        mol = Chem.MolFromSmiles(smiles)
        
        if mol is None:
            continue
            
        rdkit_mols.append((name, mol))
        
        metrics = {
            "Compound": name,
            "PubChem_CID": cmpd.cid,
            "SMILES": smiles,
            "MW": round(Descriptors.MolWt(mol), 2),
            "LogP": round(Descriptors.MolLogP(mol), 2),
            "H-Donor_Count": Descriptors.NumHDonors(mol),
            "H-Acceptor_Count": Descriptors.NumHAcceptors(mol),
            "TPSA": round(Descriptors.TPSA(mol), 2)
        }
        ligand_data.append(metrics)
        print(f"[+] Processed: {name} | SMILES: {smiles}")

    if rdkit_mols:
        img_path = os.path.join(LIGAND_DIR, "ligand_structures.png")
        img = Draw.MolsToGridImage(
            [m for _, m in rdkit_mols],
            legends=[n for n, _ in rdkit_mols],
            molsPerRow=2,
            subImgSize=(300, 300)
        )
        img.save(img_path)
        print(f"[+] Saved ligand 2D grid image to {img_path}")

    return pd.DataFrame(ligand_data)

# ==========================================
# MAIN EXECUTION PIPELINE
# ==========================================
if __name__ == "__main__":
    print("="*60)
    print("AEDES AEGYPTI RECEPTOR & LIGAND DATA RETRIEVAL PIPELINE")
    print("="*60)
    
    protein_summary = []
    for target_name, uniprot_id in TARGET_RECEPTORS.items():
        pdb_path = fetch_alphafold_structure(uniprot_id, target_name)
        if pdb_path:
            metrics = analyze_protein_structure(pdb_path)
            metrics["Target"] = target_name
            metrics["UniProt_ID"] = uniprot_id
            protein_summary.append(metrics)
            
    df_receptors = pd.DataFrame(protein_summary)
    df_ligands = fetch_and_analyze_ligands(TARGET_LIGANDS)
    
    print("\n" + "="*60)
    print("RECEPTOR PROTEIN SUMMARY (AlphaFold / Biopython)")
    print("="*60)
    print(df_receptors.to_string(index=False) if not df_receptors.empty else "No receptor data available.")
    
    print("\n" + "="*60)
    print("LIGAND CHEMICAL DESCRIPTORS (PubChem / RDKit)")
    print("="*60)
    print(df_ligands.to_string(index=False))