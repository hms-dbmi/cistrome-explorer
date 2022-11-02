### Create conda environment

```sh
conda env create -f environment.yml
```

### Activate conda environment

```sh
conda activate mira
```

### Aggregate Multivec

```sh
clodius aggregate multivec \
    --chromsizes-filename mm10.chrom.sizes \
    --starting-resolution 1000 \
    ./output/e18_mouse_brain_10x_dataset_500_random_rows.hdf5
```