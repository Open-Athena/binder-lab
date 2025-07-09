# binder-lab
Suite for evaluation of designed protein binders

## Setup

### Copy or link AF3 weights
You will need to acquire AF3 weights from Google. Then:
```
cp /path/to/your/af3.bin.zst docker/af3
```

### Build docker containers
```
docker build -t binder-lab-af3:latest docker/af3
docker build -t binder-lab-boltz:latest docker/boltz
docker build -t binder-lab-metrics:latest docker/metrics
```

### Install the package
```
pip install -e .
```

## Running
First you need to make a yaml file with your designs.

This step will depend on the method you used to generate designs. For example, we have a
script that will generate a yaml from a directory containing cif files for designs along
with npz files with a mask indicating what residues are designed.

```
DESIGN_DIR=test/data/design_dir_cif_npz/oqo-1
# or:
# DESIGN_DIR=test/data/design_dir_cif_npz/glp1
OUT_DIR=results/$(basename $DESIGN_DIR)
mkdir -p $OUT_DIR
python scripts/ingest_from_design_dir.py $DESIGN_DIR $OUT_DIR/designs.yaml
```

```
OUT_DIR=results/$(basename $DESIGN_DIR)
mkdir -p $OUT_DIR
python scripts/ingest_from_design_dir.py $DESIGN_DIR $OUT_DIR/designs.yaml
```

Take a look at $OUT_DIR/designs.yaml and make sure everything looks right. Then copy a config
file that indicates what predictors to run and what metrics to calculate.

```
cp examples/config1.yaml $OUT_DIR/config.yaml
```

Now invoke snakemake:
```
snakemake --cores all --container-backend docker --config workdir=$OUT_DIR
```

## Run unit tests
```
pytest -sv test/
```