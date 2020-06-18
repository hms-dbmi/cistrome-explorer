import json

def bigwigs_to_manifest(
    input_bigwig_files,
    input_metadata_files,
    output_file
):

    manifest_json = {
        "input_bigwig_files": input_bigwig_files,
        "input_metadata_files": input_metadata_files,
    }

    with open(output_file, "w") as f:
        json.dump(manifest_json, f)


if __name__ == "__main__":

    bigwigs_to_manifest(
        list(snakemake.input["bigwigs"]),
        list(snakemake.input["metadata"]),
        snakemake.output[0]
    )