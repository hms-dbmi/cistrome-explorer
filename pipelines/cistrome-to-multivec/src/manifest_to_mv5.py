import h5py
import bbi
import negspy.coordinates as nc
import numpy as np
import math
import argparse
import json
from tqdm import tqdm
import resource

from utils import metadata_json_to_row_info

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
    chroms_group = f.create_group("chroms")

    # Set info attributes
    info_group.attrs['tile-size'] = 256

    # Prepare to fill in chroms dataset
    chromosomes = nc.get_chromorder('hg38')
    chromosomes = chromosomes[:25] # TODO: should more than chr1-chrM be used?
    num_chromosomes = len(chromosomes)
    chroms_length_arr = np.array([ nc.get_chrominfo('hg38').chrom_lengths[x] for x in chromosomes ], dtype="i8")
    chroms_name_arr = np.array(chromosomes, dtype="S23")

    chromosomes_set = set(chromosomes)
    chrom_name_to_length = dict(zip(chromosomes, chroms_length_arr))

    # Fill in chroms dataset entries "length" and "name"
    chroms_group.create_dataset("length", data=chroms_length_arr)
    chroms_group.create_dataset("name", data=chroms_name_arr)

    
    # Prepare to fill in resolutions dataset
    resolutions = [ 1000*(2**x) for x in range(15)]
    lowest_resolution = resolutions[-1]
    
    # Create each resolution group.
    for resolution in resolutions:
        resolution_group = resolutions_group.create_group(str(resolution))
        # TODO: remove the unnecessary "values" layer
        resolution_values_group = resolution_group.create_group("values")
        
        # Create each chromosome dataset.
        for chr_name, chr_len in zip(chromosomes, chroms_length_arr):
            chr_shape = (math.ceil(chr_len / resolution), num_samples)
            resolution_values_group.create_dataset(chr_name, chr_shape, dtype="f4", fillvalue=np.nan, compression='gzip')
    
    # Fill in data for each bigwig file.
    for bw_index, bw_file in tqdm(list(enumerate(input_bigwig_files)), desc='bigwigs'):
        if bbi.is_bigwig(bw_file):
            chromsizes = bbi.chromsizes(bw_file)
            matching_chromosomes = set(chromsizes.keys()).intersection(chromosomes_set)

            # Fill in data for each resolution of a bigwig file.
            for resolution in resolutions:
                # Fill in data for each chromosome of a resolution of a bigwig file.
                for chr_name in matching_chromosomes:
                    chr_len = chrom_name_to_length[chr_name]
                    chr_shape = (math.ceil(chr_len / resolution), num_samples)
                    arr = bbi.fetch(bw_file, chr_name, 0, chr_len, chr_shape[0], summary="sum")
                    resolutions_group[str(resolution)]["values"][chr_name][:,bw_index] = arr
        else:
            print(f"{bw_file} not is_bigwig")

        f.flush()

    f.close()

    max_mem = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    print(max_mem)
    
    # Append metadata to the top resolution row_infos attribute.
    row_infos = []
    for metadata_index, metadata_file in enumerate(input_metadata_files):
        with open(metadata_file) as mf:
            try:
                metadata_json = json.load(mf)
            except Exception as e:
                print(f"Error loading metadata file: {metadata_file}")
                print(e)
                metadata_json = None
        row_info = metadata_json_to_row_info(metadata_json)
        row_infos.append(row_info)
    
    row_infos_encoded = str(json.dumps(row_infos))
    
    f = h5py.File(output_file, 'r+')

    info_group = f["info"]
    info_group["row_infos"] = row_infos_encoded
    
    f.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Create a multivec file.')
    parser.add_argument('-i', '--input', type=str, required=True, help='The input manifest JSON file.')
    parser.add_argument('-o', '--output', type=str, required=True, help='The output multivec file.')
    parser.add_argument('-s', '--starting-resolution', type=int, default=1000, help='The starting resolution.')
    args = parser.parse_args()

    with open(args.input) as f:
        manifest_json = json.load(f)
        input_bigwig_files = manifest_json['input_bigwig_files']
        input_metadata_files = manifest_json['input_metadata_files']

    bigwigs_to_multivec(
        input_bigwig_files,
        input_metadata_files,
        args.output,
        args.starting_resolution
    )