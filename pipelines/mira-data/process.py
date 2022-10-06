import numpy as np
import math
import h5py
from torch import chunk
import scanpy as sc

NUM_ROWS = 500
RESOLUTION = 1000

def basic_preprocess():
    # Download this data by running `mira.datasets.MouseBrainDataset()`
    data = sc.read_h5ad('mira-datasets/e18_10X_brain_dataset/e18_mouse_brain_10x_dataset.ad')
    rna_data = data[:, data.var.feature_types == 'Gene Expression']
    atac_data = data[:, data.var.feature_types == 'Peaks']

    # Basic preprocessing steps
    rna_data.var.index = rna_data.var.index.str.upper()
    rna_data.var_names_make_unique()
    rna_data = rna_data[:, ~rna_data.var.index.str.startswith('GM')]

    sc.pp.filter_cells(rna_data, min_counts = 400)
    sc.pp.filter_genes(rna_data, min_cells=15)

    rna_data.var['mt'] = rna_data.var_names.str.startswith('MT-')
    sc.pp.calculate_qc_metrics(rna_data, qc_vars=['mt'], percent_top=None,
                            log1p=False, inplace=True)

    rna_data = rna_data[rna_data.obs.pct_counts_mt < 15, :]
    rna_data = rna_data[rna_data.obs.n_genes_by_counts < 8000, :]
    sc.pp.filter_genes(rna_data, min_cells=15)

    rna_data.raw = rna_data # save raw counts
    sc.pp.normalize_total(rna_data, target_sum=1e4)
    sc.pp.log1p(rna_data)

    sc.pp.highly_variable_genes(rna_data, min_disp = -0.1)
    rna_data.layers['norm'] = rna_data.X # save normalized count data
    rna_data.X = rna_data.raw.X # and reload raw counts
    rna_data = rna_data[:, rna_data.var.highly_variable] 
    rna_data.var['exog_feature'] = rna_data.var.highly_variable # set column "exog_features" to all genes that met dispersion threshold
    rna_data.var.highly_variable = (rna_data.var.dispersions_norm > 0.8) & rna_data.var.exog_feature # set column "highly_variable" to genes that met first criteria and dispersion > 0.8

    overlapping_barcodes = np.intersect1d(rna_data.obs_names, atac_data.obs_names) # make sure barcodes are matched between modes
    atac_data = atac_data[[i for i in overlapping_barcodes],:]

    rna_data = None
    data = None

    return atac_data.to_df()

def anndata_to_multivec():
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

    df = basic_preprocess()

    # filter data for initial example
    dff = df.filter(regex='chr').head(NUM_ROWS)
    # num_rows = len(dff.index)
    
    chromScaled =  { c: math.ceil(s / RESOLUTION) for (c, s) in chromSizes }

    with h5py.File(f'./output/e18_mouse_brain_10x_dataset_{NUM_ROWS}_rows.hdf5', "w") as f:
        prev_c = None
        for column in dff.columns:
            # e.g., "chr1:3060610-3061485"
            [c, interval] = column.split(':')
            [start, end] = interval.split('-')
            start = int(start)
            end = int(end)

            if prev_c == None:
                ss = chromScaled[c]
                density = np.zeros((ss, NUM_ROWS))
                prev_c = c
            
            if c != prev_c:
                print(f'Storing {prev_c}')
                # density = density.reshape(-1, math.ceil(chromScaled[prev_c] / RESOLUTION), NUM_ROWS).sum(axis=0)
                f.create_dataset(name=prev_c, data=density, dtype='f', compression='gzip', chunks=True)
                
                # memory
                density = None
                f.flush()

                # start new
                ss = chromScaled[c]
                density = np.zeros((ss, NUM_ROWS))

                prev_c = c
            
            density[math.floor(start / RESOLUTION) : math.ceil(end / RESOLUTION)] = dff[column]

if __name__ == "__main__":
    anndata_to_multivec()
