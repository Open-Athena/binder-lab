import pandas as pd
import numpy as np
import json
from pathlib import Path
from typing import Optional, Dict, Any, Union
import warnings

import gemmi

def load_structure_predictions(csv_path: Union[str, Path], 
                        parse_structures: bool = True,
                        parse_confidence: bool = True, 
                        parse_npz: bool = True,
                        base_dir: Optional[Union[str, Path]] = None,
                        add_summaries: bool = True) -> pd.DataFrame:
    """
    Load predictions CSV file and parse all associated data files.
    
    Parameters:
    -----------
    csv_path : str or Path
        Path to the predictions CSV file
    parse_structures : bool
        Whether to parse CIF structure files using gemmi
    parse_confidence : bool
        Whether to parse confidence JSON files  
    parse_npz : bool
        Whether to parse NPZ files (plddt, pae, pde)
    base_dir : str or Path, optional
        Base directory for resolving relative paths. If None, uses CSV file directory
        
    Returns:
    --------
    pd.DataFrame
        DataFrame with parsed data added as new columns
    """
    csv_path = Path(csv_path)
    if base_dir is None:
        base_dir = csv_path.parent
    else:
        base_dir = Path(base_dir)
    
    # Load the CSV
    df = pd.read_csv(csv_path)
    print(f"Loaded {len(df)} predictions from {csv_path}")
    
    if len(df) == 0:
        return df
    
    # Parse structures
    if parse_structures:
        print("Parsing CIF structures...")
        df['structure'] = df.apply(lambda row: _parse_cif_structure(base_dir / row['cif_path']), axis=1)
    
    # Parse confidence files
    if parse_confidence:
        print("Parsing confidence files...")
        df['confidence'] = df.apply(lambda row: _parse_confidence_file(base_dir / row['confidence_path']), axis=1)
    
    # Parse NPZ files
    if parse_npz:
        print("Parsing NPZ files...")
        for npz_type in ['plddt', 'pae', 'pde']:
            col_name = f'{npz_type}_path'
            if col_name in df.columns:
                print(f"  Parsing {npz_type} files...")
                df[npz_type] = df.apply(lambda row: _parse_npz_file(base_dir / row[col_name]) if not pd.isna(row[col_name]) else None, axis=1)
    
    if add_summaries:
        df = add_summary_columns(df)

    print(f"Parsing complete. DataFrame now has {len(df.columns)} columns")
    return df


def _parse_cif_structure(cif_path: Path) -> Optional[gemmi.Structure]:
    """Parse a CIF structure file using gemmi."""        
    try:
        if not cif_path.exists():
            print(f"Warning: CIF file not found: {cif_path}")
            return None
            
        structure = gemmi.read_structure(str(cif_path))
        return structure
    except Exception as e:
        print(f"Error parsing CIF file {cif_path}: {e}")
        return None


def _parse_confidence_file(json_path: Path) -> Optional[Dict[str, Any]]:
    """Parse a confidence JSON file."""
    try:
        if not json_path.exists():
            print(f"Warning: Confidence file not found: {json_path}")
            return None
            
        with open(json_path, 'r') as f:
            confidence_data = json.load(f)
        return confidence_data
    except Exception as e:
        print(f"Error parsing confidence file {json_path}: {e}")
        return None


def _parse_npz_file(npz_path: Path) -> Optional[Dict[str, np.ndarray]]:
    """Parse an NPZ file containing numpy arrays."""
    try:
        if not npz_path.exists():
            print(f"Warning: NPZ file not found: {npz_path}")
            return None
            
        npz_data = np.load(npz_path)
        # Convert to dict for easier access
        return {key: npz_data[key] for key in npz_data.files}
    except Exception as e:
        print(f"Error parsing NPZ file {npz_path}: {e}")
        return None


def get_structure_info(structure: gemmi.Structure) -> Dict[str, Any]:
    """
    Extract useful information from a gemmi Structure object.
    
    Parameters:
    -----------
    structure : gemmi.Structure
        Parsed structure from gemmi
        
    Returns:
    --------
    dict
        Dictionary with structure information
    """
    if structure is None:
        return {}
    
    info = {
        'name': structure.name,
        'num_models': len(structure),
        'num_chains': 0,
        'num_residues': 0,
        'num_atoms': 0,
        'chain_ids': [],
    }
    
    # Count chains, residues, atoms
    for model in structure:
        info['num_chains'] += len(model)
        for chain in model:
            if chain.name not in info['chain_ids']:
                info['chain_ids'].append(chain.name)
            info['num_residues'] += len(chain)
            for residue in chain:
                info['num_atoms'] += len(residue)
    
    return info


def get_confidence_summary(confidence_data: Dict[str, Any]) -> Dict[str, Any]:
    """
    Extract summary statistics from confidence data.
    
    Parameters:
    -----------
    confidence_data : dict
        Parsed confidence JSON data
        
    Returns:
    --------
    dict
        Summary statistics
    """
    if confidence_data is None:
        return {}
    
    summary = {}
    
    # Extract common confidence metrics if they exist
    for key in ['mean_plddt', 'ptm', 'iptm', 'ranking_confidence']:
        if key in confidence_data:
            summary[key] = confidence_data[key]
    
    return summary


def get_plddt_summary(plddt_data: Dict[str, np.ndarray]) -> Dict[str, float]:
    """
    Extract summary statistics from pLDDT data.
    
    Parameters:
    -----------
    plddt_data : dict
        Parsed pLDDT NPZ data
        
    Returns:
    --------
    dict
        pLDDT summary statistics
    """
    if plddt_data is None:
        return {}
    
    summary = {}
    
    # Look for pLDDT array (common key names)
    plddt_array = None
    for key in ['plddt', 'lddt', 'confidence']:
        if key in plddt_data:
            plddt_array = plddt_data[key]
            break
    
    if plddt_array is not None:
        summary.update({
            'plddt_mean': float(np.mean(plddt_array)),
            'plddt_median': float(np.median(plddt_array)),
            'plddt_min': float(np.min(plddt_array)),
            'plddt_max': float(np.max(plddt_array)),
            'plddt_std': float(np.std(plddt_array)),
            'plddt_high_conf': float(np.mean(plddt_array > 90)),  # Fraction > 90
            'plddt_low_conf': float(np.mean(plddt_array < 50))    # Fraction < 50
        })
    
    return summary


def add_summary_columns(df: pd.DataFrame) -> pd.DataFrame:
    """
    Add summary columns extracted from parsed data.
    
    Parameters:
    -----------
    df : pd.DataFrame
        DataFrame with parsed structure and confidence data
        
    Returns:
    --------
    pd.DataFrame
        DataFrame with additional summary columns
    """
    df = df.copy()
    
    # Add structure info
    if 'structure' in df.columns:
        print("Adding structure summary columns...")
        structure_info = df['structure'].apply(get_structure_info)
        structure_df = pd.json_normalize(structure_info)
        
        # Add with prefix to avoid column name conflicts
        for col in structure_df.columns:
            df[f'struct_{col}'] = structure_df[col]
    
    # Add confidence summary
    if 'confidence' in df.columns:
        print("Adding confidence summary columns...")
        confidence_summary = df['confidence'].apply(get_confidence_summary)
        confidence_df = pd.json_normalize(confidence_summary)
        
        for col in confidence_df.columns:
            df[f'conf_{col}'] = confidence_df[col]
    
    # Add pLDDT summary
    if 'plddt' in df.columns:
        print("Adding pLDDT summary columns...")
        plddt_summary = df['plddt'].apply(get_plddt_summary)
        plddt_df = pd.json_normalize(plddt_summary)
        
        for col in plddt_df.columns:
            df[col] = plddt_df[col]
    
    return df


