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

### Copy File to Server

```sh
rsync -avz --progress filename server:filename
```

### Add to HiGlass Server

```sh
sudo docker exec higlass-container python higlass-server/manage.py ingest_tileset \
      --filename /tmp/e18_mouse_brain_10x_dataset_500_smoothed_random_rows.multires.mv5 \
      --filetype multivec \
      --datatype multivec \
      --uid e18_mouse_brain_10x_dataset_500_smoothed_random_rows \
```