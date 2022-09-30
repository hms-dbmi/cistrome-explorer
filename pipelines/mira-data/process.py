import numpy as np
import math
import h5py
import scanpy as sc

NUM_ROWS = 100
RESOLUTION = 1000

def anndata_to_multivec():
    # Download this data by running `mira.datasets.MouseBrainDataset()`
    data = sc.read_h5ad('mira-datasets/e18_10X_brain_dataset/e18_mouse_brain_10x_dataset.ad')
    atac_data = data[:, data.var.feature_types == 'Peaks']
    df = atac_data.to_df()

    # https://github.com/igvteam/igv/blob/master/genomes/sizes/mm10.chrom.sizes
    chromSizes = [
        ('chr1', 195471971),
        ('chr2', 182113224),
        ('chr3', 160039680),
        ('chr4', 156508116),
        ('chr5', 151834684),
        ('chr6', 149736546),
        ('chr7', 145441459),
        ('chr8', 129401213),
        ('chr9', 124595110),
        ('chr10', 130694993),
        ('chr11', 122082543),
        ('chr12', 120129022),
        ('chr13', 120421639),
        ('chr14', 124902244),
        ('chr15', 104043685),
        ('chr16', 98207768),
        ('chr17', 94987271),
        ('chr18', 90702639),
        ('chr19', 61431566),
        ('chrX', 171031299),
        ('chrY', 91744698),
        ('chrM', 16299)
    ]

    # filter data for initial example
    dff = df.filter(regex='chr').head(NUM_ROWS)
    # num_rows = len(dff.index)

    density_dict = { c: np.zeros((math.ceil(s / 1000) * RESOLUTION, NUM_ROWS)) for (c, s) in chromSizes }

    # "chr1:3060610-3061485"
    prev_chr_str = None
    for column in dff.columns:
        [c, interval] = column.split(':')
        [start, end] = interval.split('-')
        start = int(start)
        end = int(end)

        density_dict[c][start:end+1] = dff[column]
        
        if c != prev_chr_str:
            prev_chr_str = c
            print(prev_chr_str)

    # density_dict['chr1'][3113326]

    with h5py.File(f'./e18_mouse_brain_10x_dataset_{NUM_ROWS}_rows.hdf5', "w") as f:
        for (c, s) in chromSizes:
            print(c, s, density_dict[c].shape)
            
            density = density_dict[c]
            density = density.reshape(-1, math.ceil(s / 1000), NUM_ROWS).sum(axis=0)

            f.create_dataset(name=c, data=density, compression='gzip')

            density_dict[c] = None
            
            f.flush()

if __name__ == "__main__":
    anndata_to_multivec()
    # f = h5py.File(f'./e18_mouse_brain_10x_dataset_{NUM_ROWS}_rows.hdf5', "r")
    # print(f['chr1'].size)
    # f.close()
