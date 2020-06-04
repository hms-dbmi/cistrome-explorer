import bbi
import h5py
import negspy.coordinates as nc
import numpy as np
import math

def bigwigs_to_multivec(
    input_bigwig_files,
    input_metadata_files,
    output_file,
    starting_resolution
):

    f = h5py.File(output_file, 'w')

    num_samples = len(input_bigwig_files)

    # Zip the input to create (bw, metadata) tuples
    zipped_input = zip(input_bigwig_files, input_metadata_files)

    # Create level zero groups
    info_group = f.create_group("info")
    resolutions_group = f.create_group("resolutions")
    chroms_group = f.create_group("chroms") # TODO: rename to "chromosomes"

    # Prepare to fill in chroms dataset
    chromosomes = nc.get_chromorder('hg38')
    num_chromosomes = len(chromosomes)
    chroms_length_arr = np.array([ nc.get_chrominfo('hg38').chrom_lengths[x] for x in chromosomes ], dtype="i8")
    chroms_name_arr = np.array(chromosomes, dtype="S23")

    # Fill in chroms dataset entries "length" and "name"
    chroms_group.create_dataset("length", data=chroms_length_arr)
    chroms_group.create_dataset("name", data=chroms_name_arr)

    
    # Prepare to fill in resolutions dataset
    resolutions = [ 1000*(2**x) for x in range(15)]
    
    # Create dataset for each resolution
    for resolution in resolutions:
        # Create each resolution group
        resolution_group = resolutions_group.create_group(str(resolution))
        resolution_values_group = resolution_group.create_group("values") # TODO: remove the unnecessary "values" layer
        
        for chr_name, chr_len in zip(chromosomes, chroms_length_arr):
            chr_shape = (math.ceil(chr_len / resolution), num_samples)
            chr_dataset = resolution_values_group.create_dataset(chr_name, chr_shape, dtype="f4")

            # Fill in dataset for each resolution
            for bw_index, bw_file in enumerate(input_bigwig_files):

                if bbi.is_bigwig(bw_file):
                    chromsizes = bbi.chromsizes(bw_file)

                    if chr_name in chromsizes.keys():
                        arr = bbi.fetch(bw_file, chr_name, 0, chr_len, chr_shape[0], summary="sum")
                        chr_dataset[:,bw_index] = arr
                    else:
                        pass
                        #print(f"{bw_file} does not have {chr_name}")
                else:
                    print(f"{bw_file} not is_bigwig")



if __name__ == "__main__":

    bigwigs_to_multivec(
        list(snakemake.input["bigwigs"]),
        list(snakemake.input["metadata"]),
        snakemake.output[0],
        snakemake.params["starting_resolution"]
    )