#!/usr/bin/env python3

import argparse
from pathlib import Path
from typing import Dict, List, Tuple, Optional, Union
import yaml  # type: ignore
import numpy as np
import gemmi  # type: ignore

def load_cif_structure(cif_path: Path) -> gemmi.Structure:
    """Load a CIF file into a gemmi Structure object."""
    structure = gemmi.read_structure(str(cif_path))
    structure.remove_alternative_conformations()
    structure.remove_hydrogens()
    structure.remove_waters()
    return structure

def extract_sequence_from_structure(structure: gemmi.Structure, chain_id: str) -> str:
    """Extract amino acid sequence from a specific chain in the structure."""
    three_to_one = {
        'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E', 'PHE': 'F',
        'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LYS': 'K', 'LEU': 'L',
        'MET': 'M', 'ASN': 'N', 'PRO': 'P', 'GLN': 'Q', 'ARG': 'R',
        'SER': 'S', 'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y'
    }
    
    sequence = ""
    for model in structure:
        for chain in model:
            if chain.name == chain_id:
                for residue in chain:
                    res_name = residue.name
                    if res_name in three_to_one:
                        sequence += three_to_one[res_name]
    return sequence

def count_chain_tokens(structure: gemmi.Structure, chain_id: str) -> int:
    """
    Count the number of tokens in a chain.
    For proteins, this is the number of residues.
    For ligands, this is the number of atoms.
    """
    for model in structure:
        for chain in model:
            if chain.name == chain_id:
                # Check if it's a protein chain
                sequence = extract_sequence_from_structure(structure, chain_id)
                if sequence:
                    return len(sequence)
                # If not a protein, count atoms
                atom_count = 0
                for residue in chain:
                    atom_count += len(residue)
                return atom_count
    return 0

def get_chain_info(structure: gemmi.Structure, chain_id: str, design_mask: Optional[List[bool]] = None) -> Dict[str, Union[str, Dict[str, str]]]:
    """
    Get information about a chain, determining if it's a protein or ligand.
    Returns a dict with either:
    - {'protein': {'sequence', 'designed', 'msa'}} for proteins
    - {'ligand': {'ccd'}} for ligands
    """
    # First try to get protein sequence
    sequence = extract_sequence_from_structure(structure, chain_id)
    
    # If we got a protein sequence, return protein info
    if sequence:
        if design_mask is not None:
            designed = ''.join('D' if x else '.' for x in design_mask)
        else:
            designed = '.' * len(sequence)
        return {
            'protein': {
                'id': chain_id,
                'sequence': sequence,
                'designed': designed
            }
        }
    
    # If no protein sequence found, look for ligand/CCD
    for model in structure:
        for chain in model:
            if chain.name == chain_id:
                # For ligands, we expect a single residue with a CCD code
                if len(chain) == 1:
                    return {
                        'ligand': {
                            'id': chain_id,
                            'ccd': chain[0].name
                        }
                    }
                # For multi-residue non-protein chains, use the first residue's name
                elif len(chain) > 1:
                    return {
                        'ligand': {
                            'id': chain_id,
                            'ccd': chain[0].name
                        }
                    }
    
    raise ValueError(f"Could not determine type of chain {chain_id}")

def get_all_chain_ids(structure: gemmi.Structure) -> List[str]:
    """Get all chain IDs from the structure in alphanumeric order."""
    # We only look at the first model as that's typically what we work with
    if len(structure) > 0:
        return sorted([chain.name for chain in structure[0]])
    return []

def split_design_mask(structure: gemmi.Structure, design_mask: np.ndarray) -> Dict[str, List[bool]]:
    """
    Split the concatenated design mask into per-chain masks.
    Chains are assumed to be concatenated in alphanumeric order.
    """
    chain_ids = get_all_chain_ids(structure)
    chain_masks: Dict[str, List[bool]] = {}
    
    # Get token counts for each chain
    start_idx = 0
    for chain_id in chain_ids:
        n_tokens = count_chain_tokens(structure, chain_id)
        if start_idx + n_tokens <= len(design_mask):
            chain_masks[chain_id] = design_mask[start_idx:start_idx + n_tokens].tolist()
            start_idx += n_tokens
        else:
            raise ValueError(f"Design mask length {len(design_mask)} is too short for all chains")
    
    if start_idx != len(design_mask):
        raise ValueError(f"Design mask length {len(design_mask)} doesn't match total number of tokens {start_idx}")
    
    return chain_masks

def process_design_files(design_dir: Path) -> List[Dict]:
    """Process all CIF and NPZ files in the design directory."""
    designs: List[Dict] = []
    
    # Group CIF and NPZ files by their base names
    file_groups: Dict[str, Dict[str, Path]] = {}
    for file_path in design_dir.glob("batch0_sample*"):
        # Remove extension before calculating base name
        name_without_ext = file_path.stem
        base_name = "_".join(name_without_ext.split("_")[:3])  # batch0_sample0_rank0
        if base_name not in file_groups:
            file_groups[base_name] = {}
        if file_path.suffix == '.cif':
            file_groups[base_name]['cif'] = file_path
        elif file_path.suffix == '.npz':
            if file_path.name.endswith('_metadata.npz'):
                file_groups[base_name]['npz'] = file_path

    # Process each group of files
    for base_name, files in file_groups.items():
        if 'cif' not in files or 'npz' not in files:
            continue
            
        try:
            # Load structure from CIF
            structure = load_cif_structure(files['cif'])
            
            # Load metadata from NPZ
            metadata = np.load(files['npz'])
            design_mask = metadata['design_mask']
            
            # Get all chain IDs in alphanumeric order (same as design mask concatenation)
            chain_ids = get_all_chain_ids(structure)
            if not chain_ids:
                print(f"Warning: No chains found in {base_name}")
                continue
            
            try:
                # Split design mask into per-chain masks
                chain_masks = split_design_mask(structure, design_mask)
            except ValueError as e:
                print(f"Warning: {str(e)} in {base_name}")
                continue
            
            sequences = []
            # Process each chain in alphanumeric order
            for chain_id in chain_ids:
                chain_info = get_chain_info(structure, chain_id, chain_masks[chain_id])
                sequences.append(chain_info)
            
            # Create design entry
            design = {
                'name': base_name,
                'sequences': sequences
            }
            designs.append(design)
            
        except Exception as e:
            print(f"Error processing {base_name}: {str(e)}")
            continue
    
    return designs

def main() -> None:
    parser = argparse.ArgumentParser(description='Generate YAML from design directory')
    parser.add_argument('design_dir', type=str, help='Directory containing CIF and NPZ files')
    parser.add_argument('output_yaml', type=str, help='Output YAML file path')
    args = parser.parse_args()
    
    design_dir = Path(args.design_dir)
    output_yaml = Path(args.output_yaml)
    
    # Create output directory if it doesn't exist
    output_yaml.parent.mkdir(parents=True, exist_ok=True)
    
    # Process designs
    designs = process_design_files(design_dir)
    
    # Create final YAML structure
    yaml_data = {
        'version': 1,
        'experiment_name': design_dir.name,
        'designs': designs
    }
    
    # Write YAML file
    with open(output_yaml, 'w') as f:
        yaml.safe_dump(yaml_data, f, sort_keys=False)
    
    print(f"Successfully processed {len(designs)} designs")
    print(f"YAML file written to: {output_yaml}")

if __name__ == '__main__':
    main()
