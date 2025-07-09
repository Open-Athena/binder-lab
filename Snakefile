import yaml
from pathlib import Path
import os

# Save original cwd
original_cwd = os.getcwd()

# Configuration
configfile: f"{config.get('workdir', '.')}/config.yaml"

working_dir = Path(os.path.abspath(config.get('workdir', '.')))

# Get boltz predictors from config
def get_boltz_predictors():
    boltz_predictors = []
    for predictor in config.get('structure_predictors', []):
        if predictor.get('tool') == 'boltz':
            boltz_predictors.append(predictor)
    return boltz_predictors

# Helper function to get predictor config
def get_predictor_config(predictor_name):
    for pred in config.get('structure_predictors', []):
        if pred.get('name') == predictor_name and pred.get('tool') == 'boltz':
            return pred
    raise ValueError(f"Could not find boltz predictor config for {predictor_name}")

# Helper function to get predictor names
def get_boltz_predictor_names():
    return [predictor['name'] for predictor in get_boltz_predictors()]

# Main workflow entry point
rule all:
    input:
        working_dir / "structure_predictions.csv"

# Rule to prepare individual YAML files for each boltz predictor
rule prepare_boltz_inputs:
    input:
        designs=working_dir / "designs.yaml"
    output:
        [directory(working_dir / f"{predictor['name']}_inputs") for predictor in get_boltz_predictors()]
    run:
        import yaml
        from pathlib import Path
        
        # Load designs from YAML
        with open(input.designs, 'r') as f:
            designs_data = yaml.safe_load(f)
        
        # Get boltz predictors
        boltz_predictors = get_boltz_predictors()
        
        # Create separate input directories for each boltz predictor
        for i, predictor in enumerate(boltz_predictors):
            output_dir = Path(str(output[i]))
            output_dir.mkdir(exist_ok=True)
            
            # Process each design
            for design in designs_data.get('designs', []):
                design_name = design['name']
                sequences = design['sequences']
                
                # Convert to boltz format
                boltz_yaml = {
                    'version': 1,
                    'sequences': []
                }
                
                # Process each sequence in the design
                for seq in sequences:
                    # Handle protein sequences - check if 'protein' key exists 
                    if 'protein' in seq:
                        protein_data = seq['protein']
                        protein_seq = {
                            'protein': {
                                'id': protein_data['id'],
                                'sequence': protein_data['sequence']
                            }
                        }
                        
                        # If any residues are designed, use empty MSA
                        if 'designed' in protein_data and 'D' in protein_data['designed']:
                            protein_seq['protein']['msa'] = 'empty'
                        
                        boltz_yaml['sequences'].append(protein_seq)
                    
                    # Handle ligand sequences - check if 'ligand' key exists
                    elif 'ligand' in seq:
                        ligand_data = seq['ligand']
                        ligand_seq = {
                            'ligand': {
                                'id': ligand_data['id'],
                                'ccd': ligand_data['ccd']
                            }
                        }
                        boltz_yaml['sequences'].append(ligand_seq)
                
                # Write individual YAML file for this design
                output_file = output_dir / f"{design_name}.yaml"
                with open(output_file, 'w') as f:
                    yaml.safe_dump(boltz_yaml, f, sort_keys=False)
            
            print(f"Created {len(designs_data.get('designs', []))} YAML files in {output_dir} for {predictor['name']}")

# Rule to run boltz predict on the prepared YAML files
rule boltz_predict:
    input:
        boltz_inputs=working_dir / "{predictor}_inputs"
    output:
        directory(working_dir / "{predictor}_predictions")
    container:
        f"{config.get('boltz_image', 'resources/boltz/boltz.sif')}"
    resources:
        gpu=1
    params:
        version=lambda wildcards: get_predictor_config(wildcards.predictor).get('version', 2),
        recycles=lambda wildcards: get_predictor_config(wildcards.predictor).get('recycles', 3),
        diffusion_samples=lambda wildcards: get_predictor_config(wildcards.predictor).get('diffusion_samples', 5),
        resources_dir="/resources/boltz",
    shell:
        """
        # Remove snakemake timestamp file that interferes with boltz
        rm -f {input.boltz_inputs}/.snakemake_timestamp
        
        # Build base command
        CMD="boltz predict {input.boltz_inputs} --out_dir {output} --use_msa_server"
        
        # Add version-specific flags
        if [ "{params.version}" = "1" ]; then
            CMD="$CMD --model boltz1"
        fi
        
        # Add recycling steps and number of models
        CMD="$CMD --recycling_steps {params.recycles} --diffusion_samples {params.diffusion_samples}"

        # Check that cache directory exists
        if [ ! -d "{params.resources_dir}" ]; then
            echo "Error: Cache directory {params.resources_dir} does not exist!"
            exit 1
        fi

        # Set up cache directory (model weights etc)
        CMD="$CMD --cache {params.resources_dir}"

        # Set up devices
        CMD="$CMD --devices {resources.gpu}"

        echo "Running: $CMD"
        $CMD
        """

# Rule to collect all boltz prediction file paths into a CSV
rule collect_prediction_data:
    input:
        prediction_dir=working_dir / "{predictor}_predictions",
        designs_file=working_dir / "designs.yaml"
    output:
        working_dir / "{predictor}_predictions_data.csv"
    run:
        import pandas as pd
        from pathlib import Path
        import yaml
        import json
        
        predictor_name = wildcards.predictor
        prediction_dir = Path(input.prediction_dir)
        
        print(f"Collecting prediction file paths for {predictor_name}...")
        
        # Get predictor config to determine expected files
        pred_config = get_predictor_config(predictor_name)
        diffusion_samples = pred_config.get('diffusion_samples', 5)
        
        # Load designs from YAML
        with open(input.designs_file, 'r') as f:
            designs_data = yaml.safe_load(f)
        
        print(f"Found {len(designs_data.get('designs', []))} designs")
        print(f"Expecting {diffusion_samples} models per design")
        
        # Create list to store all data
        data_rows = []
        
        # Process each design and model combination
        for design in designs_data.get('designs', []):
            design_name = design['name']
            design_dict_str = json.dumps(design, sort_keys=True)  # Convert design dict to string
            
            for model_idx in range(diffusion_samples):
                # Compute expected file paths based on designs.yaml and config
                base_pattern = f"{design_name}_model_{model_idx}"
                
                # Expected file paths (relative to prediction_dir) - common to all versions
                expected_files = {
                    'cif': f"boltz_results_{predictor_name}_inputs/predictions/{design_name}/{design_name}_model_{model_idx}.cif",
                    'confidence': f"boltz_results_{predictor_name}_inputs/predictions/{design_name}/confidence_{design_name}_model_{model_idx}.json",
                    'plddt': f"boltz_results_{predictor_name}_inputs/predictions/{design_name}/plddt_{design_name}_model_{model_idx}.npz",
                }
                
                # Additional files only for boltz2 and higher
                if pred_config.get('version', 2) >= 2:
                    expected_files.update({
                        'pae': f"boltz_results_{predictor_name}_inputs/predictions/{design_name}/pae_{design_name}_model_{model_idx}.npz", 
                        'pde': f"boltz_results_{predictor_name}_inputs/predictions/{design_name}/pde_{design_name}_model_{model_idx}.npz"
                    })
                
                # Create row data with metadata
                row_data = {
                    'predictor': predictor_name,
                    'design_name': design_name,
                    'design_dict': design_dict_str,
                    'model_idx': model_idx,
                    'base_pattern': base_pattern,
                }
                
                # Check file existence and add paths - all files must exist
                for file_type, rel_path in expected_files.items():
                    abs_path = (prediction_dir / rel_path).absolute()
                    
                    if not abs_path.exists():
                        raise FileNotFoundError(f"Expected {file_type} file not found: {abs_path}")
                    
                    # Store relative path (portable)
                    row_data[f'{file_type}_path'] = str(abs_path.relative_to(working_dir))
                
                # Add empty columns for missing file types (for consistency across versions)
                all_possible_types = ['cif', 'confidence', 'plddt', 'pae', 'pde']
                for file_type in all_possible_types:
                    if f'{file_type}_path' not in row_data:
                        row_data[f'{file_type}_path'] = None
                
                data_rows.append(row_data)
                
        
        # Create DataFrame
        df = pd.DataFrame(data_rows)

        assert len(df) > 0, "No data was collected!"
        
        print(f"Created DataFrame with {len(df)} rows and {len(df.columns)} columns")
        print(f"Columns: {list(df.columns)}")
        print(f"Designs: {sorted(df['design_name'].unique()) if 'design_name' in df.columns else 'None'}")
        print(f"All expected files verified to exist")
        
        # Save as CSV
        df.to_csv(output[0], index=False)
        print(f"Saved prediction file paths to {output[0]}")

# Rule to aggregate all predictor CSV files into one master CSV
rule aggregate_predictions:
    input:
        [working_dir / f"{predictor['name']}_predictions_data.csv" for predictor in get_boltz_predictors()]
    output:
        working_dir / "structure_predictions.csv"
    run:
        import pandas as pd
        from pathlib import Path
        
        print("Aggregating all prediction data into master CSV...")
        
        # Read all individual CSV files
        all_dfs = []
        for csv_file in input:
            df = pd.read_csv(csv_file)
            all_dfs.append(df)
            print(f"  Loaded {len(df)} rows from {csv_file}")
        
        assert len(all_dfs) > 0, "No CSV files found!"

        # Concatenate all DataFrames
        combined_df = pd.concat(all_dfs, ignore_index=True)
        
        print(f"Combined DataFrame:")
        print(f"  Total rows: {len(combined_df)}")
        print(f"  Predictors: {sorted(combined_df['predictor'].unique()) if 'predictor' in combined_df.columns else 'None'}")
        print(f"  Designs: {sorted(combined_df['design_name'].unique()) if 'design_name' in combined_df.columns else 'None'}")
        print(f"  All predictions completed successfully")
        
        # Save aggregated CSV
        combined_df.to_csv(output[0], index=False)
        print(f"Saved aggregated prediction data to {output[0]}")
        print(f"Load with: df = pd.read_csv('{output[0]}')")
