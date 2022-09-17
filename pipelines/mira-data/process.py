import numpy as np
import h5py
import scanpy as sc

def main():
    data = sc.read_h5ad('mira-datasets/e18_10X_brain_dataset/e18_mouse_brain_10x_dataset.ad')
    atac_data = data[:, data.var.feature_types == 'Peaks']
    df = atac_data.to_df()

    dff = df.filter(regex='chr1:').head(1)
    columns = list(dff.columns)
    num_rows = len(dff.index)

    # CHR_FILTER = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', 'X', 'Y']
    # https://github.com/igvteam/igv/blob/master/genomes/sizes/mm10.chrom.sizes
    density_dict = { c : np.zeros((195471971, num_rows)) for c in ['chr1']}

    density_dict['chr1']

    for column in ['chr1:3060610-3061485']:
        
        chr = column.split(':')[0]
        range = column.split(':')[1].split('-')
        start = int(range[0])
        end = int(range[1])
        value = 0

        print(chr, start, end)
        print(dff[column])
    #     density_dict[chr][start:end] = dff[column]#.transpose()

    print("Writing to hdf5 file")
    mv5_data = h5py.File('./chr1.hdf5', "w")

    for c in ['chr1']:
        print(c)
        mv5_data.create_dataset((c), (195471971, num_rows), "i", data=density_dict[c], compression='gzip')

    mv5_data.close()

    print('closed')

if __name__ == "__main__":
    main()