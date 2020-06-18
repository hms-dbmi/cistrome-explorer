

def p2f(x):
    if type(x) == str:
        if x.endswith('%'):
            return float(x.strip('%'))
        return float(x)
    return x
    
def metadata_json_to_row_info(metadata_json):
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

    return row_info