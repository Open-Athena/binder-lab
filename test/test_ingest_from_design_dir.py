#!/usr/bin/env python3

import pytest
import tempfile
import yaml
from pathlib import Path
import sys
import shutil
import numpy as np

# Add the scripts directory to the path so we can import the module
sys.path.insert(0, str(Path(__file__).parent.parent / "scripts"))

from ingest_from_design_dir import (
    load_cif_structure,
    extract_sequence_from_structure,
    count_chain_tokens,
    get_chain_info,
    get_all_chain_ids,
    split_design_mask,
    process_design_files,
    main
)


class TestIngestFromDesignDir:
    """Test suite for ingest_from_design_dir.py"""
    
    @pytest.fixture
    def test_data_dir(self):
        """Path to the test data directory"""
        return Path(__file__).parent / "data" / "design_dir_cif_npz" / "oqo-1"
    
    @pytest.fixture
    def temp_output_dir(self):
        """Create a temporary directory for test outputs"""
        with tempfile.TemporaryDirectory() as temp_dir:
            yield Path(temp_dir)
    
    def test_load_cif_structure(self, test_data_dir):
        """Test loading CIF files into gemmi structures"""
        cif_path = test_data_dir / "batch0_sample0_rank0.cif"
        structure = load_cif_structure(cif_path)
        
        assert len(structure) == 1  # Should have one model
        assert len(structure[0]) == 2  # Should have two chains (A and B)
        
        chain_names = [chain.name for chain in structure[0]]
        assert "A" in chain_names
        assert "B" in chain_names
    
    def test_extract_sequence_from_structure(self, test_data_dir):
        """Test extracting protein sequences from structures"""
        cif_path = test_data_dir / "batch0_sample0_rank0.cif"
        structure = load_cif_structure(cif_path)
        
        # Test protein chain A
        sequence_a = extract_sequence_from_structure(structure, "A")
        assert len(sequence_a) == 100  # Should have 100 residues
        assert all(aa in "ACDEFGHIKLMNPQRSTVWY" for aa in sequence_a)  # Valid amino acids
        
        # Test ligand chain B (should return empty string for non-protein)
        sequence_b = extract_sequence_from_structure(structure, "B")
        assert sequence_b == ""  # Ligand chain should not have protein sequence
    
    def test_count_chain_tokens(self, test_data_dir):
        """Test counting tokens in protein and ligand chains"""
        cif_path = test_data_dir / "batch0_sample0_rank0.cif"
        structure = load_cif_structure(cif_path)
        
        # Protein chain should count residues
        tokens_a = count_chain_tokens(structure, "A")
        assert tokens_a == 100
        
        # Ligand chain should count atoms
        tokens_b = count_chain_tokens(structure, "B")
        assert tokens_b > 0  # Should have some atoms
    
    def test_get_chain_info_protein(self, test_data_dir):
        """Test getting chain info for protein chains"""
        cif_path = test_data_dir / "batch0_sample0_rank0.cif"
        structure = load_cif_structure(cif_path)
        
        # Test without design mask
        info = get_chain_info(structure, "A")
        assert "protein" in info
        assert info["protein"]["id"] == "A"
        assert len(info["protein"]["sequence"]) == 100
        assert info["protein"]["designed"] == "0" * 100  # Default all zeros
        assert info["protein"]["msa"] == "empty"
        
        # Test with design mask
        design_mask = [True] * 50 + [False] * 50  # First 50 designed
        info = get_chain_info(structure, "A", design_mask)
        assert info["protein"]["designed"] == "1" * 50 + "0" * 50
    
    def test_get_chain_info_ligand(self, test_data_dir):
        """Test getting chain info for ligand chains"""
        cif_path = test_data_dir / "batch0_sample0_rank0.cif"
        structure = load_cif_structure(cif_path)
        
        info = get_chain_info(structure, "B")
        assert "ligand" in info
        assert info["ligand"]["id"] == "B"
        assert info["ligand"]["ccd"] == "OQO"  # Expected ligand code
    
    def test_get_all_chain_ids(self, test_data_dir):
        """Test getting all chain IDs in alphanumeric order"""
        cif_path = test_data_dir / "batch0_sample0_rank0.cif"
        structure = load_cif_structure(cif_path)
        
        chain_ids = get_all_chain_ids(structure)
        assert chain_ids == ["A", "B"]  # Should be in alphanumeric order
    
    def test_split_design_mask(self, test_data_dir):
        """Test splitting concatenated design mask into per-chain masks"""
        cif_path = test_data_dir / "batch0_sample0_rank0.cif"
        structure = load_cif_structure(cif_path)
        
        # Load the actual design mask from the NPZ file
        npz_path = test_data_dir / "batch0_sample0_rank0_metadata.npz"
        metadata = np.load(npz_path)
        design_mask = metadata['design_mask']
        
        chain_masks = split_design_mask(structure, design_mask)
        
        # Should have masks for both chains
        assert "A" in chain_masks
        assert "B" in chain_masks
        
        # Check lengths match expected token counts
        assert len(chain_masks["A"]) == 100  # Protein residues
        assert len(chain_masks["B"]) > 0     # Ligand atoms
        
        # Total length should match original mask
        total_length = len(chain_masks["A"]) + len(chain_masks["B"])
        assert total_length == len(design_mask)
    
    def test_process_design_files(self, test_data_dir):
        """Test processing all design files in a directory"""
        designs = process_design_files(test_data_dir)
        
        # Should process both sample files
        assert len(designs) == 2
        
        # Check structure of first design
        design0 = designs[0]
        assert "name" in design0
        assert "sequences" in design0
        assert design0["name"] in ["batch0_sample0_rank0", "batch0_sample1_rank0"]
        
        # Should have two sequences (protein and ligand)
        sequences = design0["sequences"]
        assert len(sequences) == 2
        
        # Find protein and ligand sequences
        protein_seq = None
        ligand_seq = None
        for seq in sequences:
            if "protein" in seq:
                protein_seq = seq["protein"]
            elif "ligand" in seq:
                ligand_seq = seq["ligand"]
        
        assert protein_seq is not None
        assert ligand_seq is not None
        
        # Validate protein sequence structure
        assert protein_seq["id"] == "A"
        assert len(protein_seq["sequence"]) == 100
        assert len(protein_seq["designed"]) == 100
        assert protein_seq["msa"] == "empty"
        
        # Validate ligand sequence structure
        assert ligand_seq["id"] == "B"
        assert ligand_seq["ccd"] == "OQO"
    
    def test_main_function(self, test_data_dir, temp_output_dir):
        """Test the main function end-to-end"""
        output_yaml = temp_output_dir / "test_output.yaml"
        
        # Mock sys.argv to simulate command line arguments
        import sys
        original_argv = sys.argv
        try:
            sys.argv = ["ingest_from_design_dir.py", str(test_data_dir), str(output_yaml)]
            main()
        finally:
            sys.argv = original_argv
        
        # Check that output file was created
        assert output_yaml.exists()
        
        # Load and validate the YAML content
        with open(output_yaml, 'r') as f:
            yaml_data = yaml.safe_load(f)
        
        # Validate top-level structure
        assert yaml_data["version"] == 1
        assert yaml_data["experiment_name"] == "oqo-1"
        assert "designs" in yaml_data
        assert len(yaml_data["designs"]) == 2
        
        # Validate first design
        design = yaml_data["designs"][0]
        assert "name" in design
        assert "sequences" in design
        assert len(design["sequences"]) == 2
    
    def test_invalid_design_mask_length(self, test_data_dir):
        """Test error handling for invalid design mask length"""
        cif_path = test_data_dir / "batch0_sample0_rank0.cif"
        structure = load_cif_structure(cif_path)
        
        # Create a design mask that's too short
        invalid_mask = np.array([True] * 50)  # Too short
        
        with pytest.raises(ValueError, match="Design mask length.*is too short"):
            split_design_mask(structure, invalid_mask)
        
        # Create a design mask that's too long
        invalid_mask = np.array([True] * 200)  # Too long
        
        with pytest.raises(ValueError, match="Design mask length.*doesn't match"):
            split_design_mask(structure, invalid_mask)
    
    def test_nonexistent_chain(self, test_data_dir):
        """Test error handling for nonexistent chains"""
        cif_path = test_data_dir / "batch0_sample0_rank0.cif"
        structure = load_cif_structure(cif_path)
        
        # Test extracting sequence from nonexistent chain
        sequence = extract_sequence_from_structure(structure, "Z")
        assert sequence == ""
        
        # Test counting tokens from nonexistent chain
        tokens = count_chain_tokens(structure, "Z")
        assert tokens == 0
        
        # Test getting chain info from nonexistent chain
        with pytest.raises(ValueError, match="Could not determine type of chain Z"):
            get_chain_info(structure, "Z")


def test_cli_integration(tmp_path):
    """Test command line interface integration"""
    # Copy test data to temp directory
    test_data_source = Path(__file__).parent / "data" / "design_dir_cif_npz" / "oqo-1"
    test_data_dir = tmp_path / "test_data"
    shutil.copytree(test_data_source, test_data_dir)
    
    output_yaml = tmp_path / "output.yaml"
    
    # Test running the script via command line
    import subprocess
    import sys
    
    script_path = Path(__file__).parent.parent / "scripts" / "ingest_from_design_dir.py"
    result = subprocess.run([
        sys.executable, str(script_path), 
        str(test_data_dir), str(output_yaml)
    ], capture_output=True, text=True)
    
    assert result.returncode == 0
    assert "Successfully processed 2 designs" in result.stdout
    assert output_yaml.exists()


if __name__ == "__main__":
    pytest.main([__file__])
