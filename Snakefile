import yaml
from pathlib import Path

# Configuration
configfile: f"{config.get('workdir', '.')}/config.yaml"

# Define workdir from config
workdir: config.get('workdir', '.')

# Get boltz predictors from config
def get_boltz_predictors():
    boltz_predictors = []
    for predictor in config.get('structure_predictors', []):
        if predictor.get('tool') == 'boltz':
            boltz_predictors.append(predictor)
    return boltz_predictors

# Main workflow entry point
rule all:
    input:
        [f"{predictor['name']}_predictions" for predictor in get_boltz_predictors()]

# Rule to prepare individual YAML files for each boltz predictor
rule prepare_boltz_inputs:
    input:
        designs="designs.yaml"
    output:
        [directory(f"{predictor['name']}_inputs") for predictor in get_boltz_predictors()]
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

# Helper function to get predictor config
def get_predictor_config(predictor_name):
    for pred in config.get('structure_predictors', []):
        if pred.get('name') == predictor_name and pred.get('tool') == 'boltz':
            return pred
    raise ValueError(f"Could not find boltz predictor config for {predictor_name}")

# Rule to run boltz predict on the prepared YAML files
rule boltz_predict:
    input:
        boltz_inputs="{predictor}_inputs"
    output:
        directory("{predictor}_predictions")
    container:
        f"{config.get('boltz_image', 'binder-lab-boltz.sif')}"
    params:
        version=lambda wildcards: get_predictor_config(wildcards.predictor).get('version', 2),
        recycles=lambda wildcards: get_predictor_config(wildcards.predictor).get('recycles', 3)
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
        
        # Add recycling steps
        CMD="$CMD --recycling_steps {params.recycles}"
        
        echo "Running: $CMD"
        $CMD
        """
