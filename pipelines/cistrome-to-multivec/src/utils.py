

def p2f(x):
    if type(x) == str:
        if x.endswith('%'):
            return float(x.strip('%'))
        return float(x)
    return x

def dict_get(d, keys, default):
    curr_d = d
    try:
        for k in keys:
            curr_d = curr_d[k]
        return curr_d
    except KeyError:
        return default
    except IndexError:
        return default
    
def metadata_json_to_row_info(metadata_json):

    # Flatten the metadata object
    row_info = {
        "id": dict_get(metadata_json, ["id"], None),
        "status": dict_get(metadata_json, ["status"], None),
        "treats__0__cell_line__name": dict_get(metadata_json, ["treats", 0, "cell_line__name"], None),
        "treats__0__cell_type__name": dict_get(metadata_json, ["treats", 0, "cell_type__name"], None),
        "treats__0__cell_pop__name": dict_get(metadata_json, ["treats", 0, "cell_pop__name"], None),
        "treats__0__disease_state__name": dict_get(metadata_json, ["treats", 0, "disease_state__name"], None),
        "treats__0__factor__name": dict_get(metadata_json, ["treats", 0, "factor__name"], None),
        "treats__0__is_correcting": dict_get(metadata_json, ["treats", 0, "is_correcting"], None),
        "treats__0__link": dict_get(metadata_json, ["treats", 0, "link"], None),
        "treats__0__name": dict_get(metadata_json, ["treats", 0, "name"], None),
        "treats__0__paper__journal__name": dict_get(metadata_json, ["treats", 0, "paper__journal__name"], None),
        "treats__0__paper__lab": dict_get(metadata_json, ["treats", 0, "paper__lab"], None),
        "treats__0__paper__pmid": dict_get(metadata_json, ["treats", 0, "paper__pmid"], None),
        "treats__0__paper__reference": dict_get(metadata_json, ["treats", 0, "paper__reference"], None),
        "treats__0__species__name": dict_get(metadata_json, ["treats", 0, "species__name"], None),
        "treats__0__strain__name": dict_get(metadata_json, ["treats", 0, "strain__name"], None),
        "treats__0__tissue_type__name": dict_get(metadata_json, ["treats", 0, "tissue_type__name"], None),
        "treats__0__unique_id": dict_get(metadata_json, ["treats", 0, "unique_id"], None),
        
        "qc__judge__map": dict_get(metadata_json, ["qc", "judge", "map"], None),
        "qc__judge__peaks": dict_get(metadata_json, ["qc", "judge", "peaks"], None),
        "qc__judge__fastqc": dict_get(metadata_json, ["qc", "judge", "fastqc"], None),
        "qc__judge__frip": dict_get(metadata_json, ["qc", "judge", "frip"], None),
        "qc__judge__pbc": dict_get(metadata_json, ["qc", "judge", "pbc"], None),
        "qc__judge__motif": dict_get(metadata_json, ["qc", "judge", "motif_judge"], None),
        "qc__judge__dhs": dict_get(metadata_json, ["qc", "judge", "dhs"], None),

        "qc__table__treat_number": dict_get(metadata_json, ["qc", "table", "treat_number"], None),
        "qc__table__control_number": dict_get(metadata_json, ["qc", "table", "control_number"], None),
        "qc__table__map__0": p2f(dict_get(metadata_json, ["qc", "table", "map", 0], None)),
        "qc__table__fastqc__0": dict_get(metadata_json, ["qc", "table", "fastqc", 0], None),
        "qc__table__frip__0": p2f(dict_get(metadata_json, ["qc", "table", "frip", 0], None)),
        "qc__table__map_number__0": dict_get(metadata_json, ["qc", "table", "map_number", 0], None),
        "qc__table__pbc__0": p2f(dict_get(metadata_json, ["qc", "table", "pbc", 0], None)),
        "qc__table__motif": dict_get(metadata_json, ["qc", "table", "motif"], None),
        "qc__table__dhs": p2f(dict_get(metadata_json, ["qc", "table", "dhs"], None)),
        "qc__table__raw_number__0": dict_get(metadata_json, ["qc", "table", "raw_number", 0], None),

        "qc__table__meta_orig__intron": dict_get(metadata_json, ["qc", "table", "meta_orig", "intron"], None),
        "qc__table__meta_orig__inter": dict_get(metadata_json, ["qc", "table", "meta_orig", "inter"], None),
        "qc__table__meta_orig__exon": dict_get(metadata_json, ["qc", "table", "meta_orig", "exon"], None),
        "qc__table__meta_orig__promoter": dict_get(metadata_json, ["qc", "table", "meta_orig", "promoter"], None),
    }

    return row_info