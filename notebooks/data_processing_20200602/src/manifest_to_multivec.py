import h5py
import bbi
import negspy.coordinates as nc
import numpy as np
import math
import argparse
import json


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

    f = h5py.File(output_file, 'w')

    num_samples = len(input_bigwig_files)

    # Zip the input to create (bw, metadata) tuples
    zipped_input = zip(input_bigwig_files, input_metadata_files)

    # Create level zero groups
    info_group = f.create_group("info")
    resolutions_group = f.create_group("resolutions")
    chroms_group = f.create_group("chroms") # TODO: rename to "chromosomes"

    # Set info attributes
    info_group.attrs['tile-size'] = 256

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
    lowest_resolution = resolutions[-1]
    
    # Create dataset for each resolution
    for resolution in resolutions:
        # Create each resolution group
        resolution_group = resolutions_group.create_group(str(resolution))
        resolution_values_group = resolution_group.create_group("values") # TODO: remove the unnecessary "values" layer
        
        for chr_name, chr_len in zip(chromosomes, chroms_length_arr):
            chr_shape = (math.ceil(chr_len / resolution), num_samples)
            chr_dataset = resolution_values_group.create_dataset(chr_name, chr_shape, dtype="f4", fillvalue=np.nan, compression='gzip')

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
    
    # Append metadata to the top resolution row_infos attribute.
    row_infos = []
    for metadata_index, metadata_file in enumerate(input_metadata_files):
        with open(metadata_file) as f:
            metadata_json = json.load(f)
        
        try:
            treats_0 = metadata_json["treats"][0]
        except IndexError:
            treats_0 = dict()
        
        # Flatten the metadata object
        row_info = {
            "id": metadata_json.get("id", None),
            "status": metadata_json.get("status", None),
            "treats__0__cell_line__name": treats_0.get("cell_line__name", None),
            "treats__0__cell_type__name": treats_0.get("cell_type__name", None),
            "treats__0__cell_pop__name": treats_0.get("cell_pop__name", None),
            "treats__0__disease_state__name": treats_0.get("disease_state__name", None),
            "treats__0__factor__name": treats_0.get("factor__name", None),
            "treats__0__is_correcting": treats_0.get("is_correcting", None),
            "treats__0__link": treats_0.get("link", None),
            "treats__0__name": treats_0.get("name", None),
            "treats__0__paper__journal__name": treats_0.get("paper__journal__name", None),
            "treats__0__paper__lab": treats_0.get("paper__lab", None),
            "treats__0__paper__pmid": treats_0.get("paper__pmid", None),
            "treats__0__paper__reference": treats_0.get("paper__reference", None),
            "treats__0__species__name": treats_0.get("species__name", None),
            "treats__0__strain__name": treats_0.get("strain__name", None),
            "treats__0__tissue_type__name": treats_0.get("tissue_type__name", None),
            "treats__0__unique_id": treats_0.get("unique_id", None),
            
            "qc__judge__map": metadata_json["qc"]["judge"]["map"],
            "qc__judge__peaks": metadata_json["qc"]["judge"]["peaks"],
            "qc__judge__fastqc": metadata_json["qc"]["judge"]["fastqc"],
            "qc__judge__frip": metadata_json["qc"]["judge"]["frip"],
            "qc__judge__pbc": metadata_json["qc"]["judge"]["pbc"],
            "qc__judge__motif": metadata_json["qc"]["judge"]["motif_judge"],
            "qc__judge__dhs": metadata_json["qc"]["judge"]["dhs"],

            "qc__table__treat_number": metadata_json["qc"]["table"]["treat_number"],
            "qc__table__control_number": metadata_json["qc"]["table"]["control_number"],
            "qc__table__map__0": p2f(metadata_json["qc"]["table"]["map"][0]),
            "qc__table__fastqc__0": metadata_json["qc"]["table"]["fastqc"][0],
            "qc__table__frip__0": p2f(metadata_json["qc"]["table"]["frip"][0]),
            "qc__table__map_number__0": metadata_json["qc"]["table"]["map_number"][0],
            "qc__table__pbc__0": p2f(metadata_json["qc"]["table"]["pbc"][0]),
            "qc__table__motif": metadata_json["qc"]["table"]["motif"],
            "qc__table__dhs": p2f(metadata_json["qc"]["table"]["dhs"]),
            "qc__table__raw_number__0": metadata_json["qc"]["table"]["raw_number"][0],

            "qc__table__meta_orig__intron": metadata_json["qc"]["table"]["meta_orig"]["intron"],
            "qc__table__meta_orig__inter": metadata_json["qc"]["table"]["meta_orig"]["inter"],
            "qc__table__meta_orig__exon": metadata_json["qc"]["table"]["meta_orig"]["exon"],
            "qc__table__meta_orig__promoter": metadata_json["qc"]["table"]["meta_orig"]["promoter"],
        }

        row_infos.append(row_info)
    row_infos_encoded = [ str(json.dumps(r)).encode() for r in row_infos ]
    resolutions_group[str(lowest_resolution)].attrs["row_infos"] = row_infos_encoded

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