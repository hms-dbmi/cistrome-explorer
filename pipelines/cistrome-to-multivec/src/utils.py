

def p2f(x):
    if type(x) == str:
        if x.endswith('%'):
            return float(x.strip('%'))
        return float(x)
    return x
    
def metadata_json_to_row_info(metadata_json):

    def dict_get(keys, expected_type):
        curr_d = metadata_json
        try:
            for k in keys:
                curr_d = curr_d[k]

            if curr_d == None:
                return None
            elif type(curr_d) != expected_type:
                try:
                    curr_d = expected_type(curr_d)
                    return curr_d
                except ValueError:
                    return None
            return curr_d
        except KeyError:
            return None
        except IndexError:
            return None
        return None

    # Flatten the metadata object
    row_info = {
        "id": dict_get(["id"], str),
        "status": dict_get(["status"], str),
        "treats__0__cell_line__name": dict_get(["treats", 0, "cell_line__name"], str),
        "treats__0__cell_type__name": dict_get(["treats", 0, "cell_type__name"], str),
        "treats__0__cell_pop__name": dict_get(["treats", 0, "cell_pop__name"], str),
        "treats__0__disease_state__name": dict_get(["treats", 0, "disease_state__name"], str),
        "treats__0__factor__name": dict_get(["treats", 0, "factor__name"], str),
        "treats__0__is_correcting": dict_get(["treats", 0, "is_correcting"], bool),
        "treats__0__link": dict_get(["treats", 0, "link"], str),
        "treats__0__name": dict_get(["treats", 0, "name"], str),
        "treats__0__paper__journal__name": dict_get(["treats", 0, "paper__journal__name"], str),
        "treats__0__paper__lab": dict_get(["treats", 0, "paper__lab"], str),
        "treats__0__paper__pmid": dict_get(["treats", 0, "paper__pmid"], str),
        "treats__0__paper__reference": dict_get(["treats", 0, "paper__reference"], str),
        "treats__0__species__name": dict_get(["treats", 0, "species__name"], str),
        "treats__0__strain__name": dict_get(["treats", 0, "strain__name"], str),
        "treats__0__tissue_type__name": dict_get(["treats", 0, "tissue_type__name"], str),
        "treats__0__unique_id": dict_get(["treats", 0, "unique_id"], str),
        
        "qc__judge__map": dict_get(["qc", "judge", "map"], bool),
        "qc__judge__peaks": dict_get(["qc", "judge", "peaks"], bool),
        "qc__judge__fastqc": dict_get(["qc", "judge", "fastqc"], bool),
        "qc__judge__frip": dict_get(["qc", "judge", "frip"], bool),
        "qc__judge__pbc": dict_get(["qc", "judge", "pbc"], bool),
        "qc__judge__motif": dict_get(["qc", "judge", "motif_judge"], bool),
        "qc__judge__dhs": dict_get(["qc", "judge", "dhs"], bool),

        "qc__table__treat_number": dict_get(["qc", "table", "treat_number"], int),
        "qc__table__control_number": dict_get(["qc", "table", "control_number"], int),
        "qc__table__map__0": p2f(dict_get(["qc", "table", "map", 0], str)),
        "qc__table__fastqc__0": dict_get(["qc", "table", "fastqc", 0], int),
        "qc__table__frip__0": p2f(dict_get(["qc", "table", "frip", 0], str)),
        "qc__table__map_number__0": dict_get(["qc", "table", "map_number", 0], int),
        "qc__table__pbc__0": p2f(dict_get(["qc", "table", "pbc", 0], str)),
        "qc__table__motif": dict_get(["qc", "table", "motif"], bool),
        "qc__table__dhs": p2f(dict_get(["qc", "table", "dhs"], str)),
        "qc__table__raw_number__0": dict_get(["qc", "table", "raw_number", 0], int),

        "qc__table__meta_orig__intron": dict_get(["qc", "table", "meta_orig", "intron"], float),
        "qc__table__meta_orig__inter": dict_get(["qc", "table", "meta_orig", "inter"], float),
        "qc__table__meta_orig__exon": dict_get(["qc", "table", "meta_orig", "exon"], float),
        "qc__table__meta_orig__promoter": dict_get(["qc", "table", "meta_orig", "promoter"], float),
    }

    return row_info