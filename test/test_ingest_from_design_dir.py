#!/usr/bin/env python3

import pytest
import tempfile
import yaml
from pathlib import Path
import sys
import shutil
import subprocess
import gemmi


class TestIngestFromDesignDirIntegration:
    """Integration test suite for ingest_from_design_dir.py script"""
    
    @pytest.fixture
    def test_data_dir(self):
        """Path to the test data directory"""
        return Path(__file__).parent / "data" / "design_dir_cif_npz" / "oqo-1"
    
    @pytest.fixture
    def temp_output_dir(self):
        """Create a temporary directory for test outputs"""
        with tempfile.TemporaryDirectory() as temp_dir:
            yield Path(temp_dir)
    
    def _extract_protein_sequence_from_cif(self, cif_path, chain_id="A"):
        """Extract protein sequence from CIF file"""
        structure = gemmi.read_structure(str(cif_path))
        structure.remove_alternative_conformations()
        structure.remove_hydrogens()
        structure.remove_waters()
        
        sequence = ""
        for model in structure:
            for chain in model:
                if chain.name == chain_id:
                    for residue in chain:
                        try:
                            # Use gemmi's built-in residue info to get one-letter code
                            res_info = gemmi.find_tabulated_residue(residue.name)
                            if res_info.is_amino_acid():
                                sequence += res_info.one_letter_code
                        except RuntimeError:
                            raise
        return sequence
    
    def _get_expected_sequences(self, test_data_dir):
        """Load CIF files and extract expected sequences"""
        expected_sequences = {}
        
        # Extract sequences from both CIF files
        for cif_file in ["batch0_sample0_rank0.cif", "batch0_sample1_rank0.cif"]:
            cif_path = test_data_dir / cif_file
            design_name = cif_file.replace(".cif", "")
            protein_sequence = self._extract_protein_sequence_from_cif(cif_path, "A")
            expected_sequences[design_name] = protein_sequence
            
        return expected_sequences
    
    def test_script_runs_successfully(self, test_data_dir, temp_output_dir):
        """Test that the script runs successfully and produces output"""
        output_yaml = temp_output_dir / "test_output.yaml"
        script_path = Path(__file__).parent.parent / "scripts" / "ingest_from_design_dir.py"
        
        # Run the script
        result = subprocess.run([
            sys.executable, str(script_path),
            str(test_data_dir), str(output_yaml)
        ], capture_output=True, text=True)
        
        # Check that script ran successfully
        assert result.returncode == 0, f"Script failed with error: {result.stderr}"
        assert "Successfully processed 2 designs" in result.stdout
        assert output_yaml.exists(), "Output YAML file was not created"
    
    def test_yaml_output_structure_and_content(self, test_data_dir, temp_output_dir):
        """Test that the YAML output has correct structure and expected content"""
        output_yaml = temp_output_dir / "test_output.yaml"
        script_path = Path(__file__).parent.parent / "scripts" / "ingest_from_design_dir.py"
        
        # Get expected sequences from CIF files
        expected_sequences = self._get_expected_sequences(test_data_dir)
        
        # Run the script
        result = subprocess.run([
            sys.executable, str(script_path),
            str(test_data_dir), str(output_yaml)
        ], capture_output=True, text=True)
        
        assert result.returncode == 0
        
        # Load and validate the YAML content
        with open(output_yaml, 'r') as f:
            yaml_data = yaml.safe_load(f)
        
        # Validate top-level structure
        assert yaml_data["version"] == 1
        assert yaml_data["experiment_name"] == "oqo-1"
        assert "designs" in yaml_data
        assert len(yaml_data["designs"]) == 2
        
        # Expected data properties
        expected_ligand_ccd = "OQO"
        expected_protein_length = 100
        
        # Validate both designs
        design_names = {design["name"] for design in yaml_data["designs"]}
        assert "batch0_sample0_rank0" in design_names
        assert "batch0_sample1_rank0" in design_names
        
        for design in yaml_data["designs"]:
            assert "name" in design
            assert "sequences" in design
            assert len(design["sequences"]) == 2  # Should have protein and ligand
            
            design_name = design["name"]
            expected_protein_sequence = expected_sequences[design_name]
            
            # Find protein and ligand sequences
            protein_seq = None
            ligand_seq = None
            for seq in design["sequences"]:
                if "protein" in seq:
                    protein_seq = seq["protein"]
                elif "ligand" in seq:
                    ligand_seq = seq["ligand"]
            
            assert protein_seq is not None, "Missing protein sequence"
            assert ligand_seq is not None, "Missing ligand sequence"
            
            # Validate protein sequence structure and content
            assert protein_seq["id"] == "A"
            assert protein_seq["sequence"] == expected_protein_sequence, f"Protein sequence mismatch for {design_name}"
            assert len(protein_seq["designed"]) == expected_protein_length
            # Design mask should be all dots and D's
            assert all(c in "D." for c in protein_seq["designed"])
            
            # Validate ligand sequence structure
            assert ligand_seq["id"] == "B"
            assert ligand_seq["ccd"] == expected_ligand_ccd
    
    def test_script_with_nonexistent_directory(self, temp_output_dir):
        """Test that script handles nonexistent input directory gracefully"""
        nonexistent_dir = temp_output_dir / "nonexistent"
        output_yaml = temp_output_dir / "test_output.yaml"
        script_path = Path(__file__).parent.parent / "scripts" / "ingest_from_design_dir.py"
        
        # Run the script with nonexistent directory
        result = subprocess.run([
            sys.executable, str(script_path),
            str(nonexistent_dir), str(output_yaml)
        ], capture_output=True, text=True)
        
        # Should succeed but process 0 designs
        assert result.returncode == 0
        assert "Successfully processed 0 designs" in result.stdout
        assert output_yaml.exists()  # YAML file is still created


def test_cli_integration_with_copy(tmp_path):
    """Test command line interface integration with copied test data"""
    # Copy test data to temp directory
    test_data_source = Path(__file__).parent / "data" / "design_dir_cif_npz" / "oqo-1"
    test_data_dir = tmp_path / "test_data"
    shutil.copytree(test_data_source, test_data_dir)
    
    output_yaml = tmp_path / "output.yaml"
    script_path = Path(__file__).parent.parent / "scripts" / "ingest_from_design_dir.py"
    
    # Test running the script via command line
    result = subprocess.run([
        sys.executable, str(script_path), 
        str(test_data_dir), str(output_yaml)
    ], capture_output=True, text=True)
    
    assert result.returncode == 0
    assert "Successfully processed 2 designs" in result.stdout
    assert output_yaml.exists()
    
    # Quick validation of output structure
    with open(output_yaml, 'r') as f:
        yaml_data = yaml.safe_load(f)
    
    assert yaml_data["version"] == 1
    assert yaml_data["experiment_name"] == "test_data"
    assert len(yaml_data["designs"]) == 2


if __name__ == "__main__":
    pytest.main([__file__])
