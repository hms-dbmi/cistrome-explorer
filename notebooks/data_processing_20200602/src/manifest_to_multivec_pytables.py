import tables
import bbi
import negspy.coordinates as nc
import numpy as np
import math
import argparse
import json

from utils import metadata_json_to_row_info

atom_i8 = tables.Atom.from_dtype(np.dtype('i8'))
atom_f4 = tables.Atom.from_dtype(np.dtype('f4'), dflt=np.nan)
atom_S23 = tables.Atom.from_dtype(np.dtype('S23'))


def p2f(x):
    if type(x) == str:
        if x.endswith('%'):
            return float(x.strip('%'))
        return float(x)
    return x

def bigwigs_to_multivec(
    input_bigwig_files,
    input_metadata_files,
    output_file,
    starting_resolution
):

    f = tables.open_file(output_file, mode='w')

    num_samples = len(input_bigwig_files)

    # Zip the input to create (bw, metadata) tuples
    zipped_input = zip(input_bigwig_files, input_metadata_files)

    # Create level zero groups
    info_group = f.create_group("/", "info")
    resolutions_group = f.create_group("/", "resolutions")
    chroms_group = f.create_group("/", "chroms") # TODO: rename to "chromosomes"

    # Set info attributes
    # Note: switched this to tile_size from tile-size
    f.set_node_attr(info_group, 'tile_size', 256)

    # Prepare to fill in chroms dataset
    chromosomes = nc.get_chromorder('hg38')
    num_chromosomes = len(chromosomes)
    chroms_length_arr = np.array([ nc.get_chrominfo('hg38').chrom_lengths[x] for x in chromosomes ], dtype="i8")
    chroms_name_arr = np.array(chromosomes, dtype="S23")

    # Fill in chroms dataset entries "length" and "name"
    f.create_array(chroms_group, "length", obj=chroms_length_arr, atom=atom_i8)
    f.create_array(chroms_group, "name", obj=chroms_name_arr, atom=atom_S23)

    
    # Prepare to fill in resolutions dataset
    resolutions = [ 1000*(2**x) for x in range(15)]
    lowest_resolution = resolutions[-1]
    
    # Create dataset for each resolution
    for resolution in resolutions:
        # Create each resolution group
        resolution_group = f.create_group(resolutions_group, str(resolution))
        resolution_values_group = f.create_group(resolution_group, "values") # TODO: remove the unnecessary "values" layer
        
        for chr_name, chr_len in zip(chromosomes, chroms_length_arr):
            chr_shape = (math.ceil(chr_len / resolution), num_samples)
            chr_dataset = f.create_array(resolution_values_group, chr_name, shape=chr_shape, atom=atom_f4)

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
            f.flush()
    
    # Append metadata to the top resolution row_infos attribute.
    row_infos = []
    for metadata_index, metadata_file in enumerate(input_metadata_files):
        with open(metadata_file) as mf:
            metadata_json = json.load(mf)
        row_info = metadata_json_to_row_info(metadata_json)
        row_infos.append(row_info)
    
    row_infos_encoded = [ str(json.dumps(r)).encode() for r in row_infos ]
    # TODO: how to deal with list attributes in pytables? by default they get converted to pickles...
    f.set_node_attr(resolution_group, 'row_infos', row_infos_encoded)


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