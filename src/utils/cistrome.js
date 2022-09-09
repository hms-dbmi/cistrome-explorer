const SEARCH_ITEM_LIMIT = 999999;

/*
 * Types for Cistrome DB Toolkit API.
 */
export const CISTROME_API_TYPES = Object.freeze({
    INTERVAL: "INTERVAL",
    GENE: "GENE",
    FACTOR: "FACTOR",
    PEAKSET: "PEAKSET"
});

/*
 * Color mapping for different types of Cistrome DB Toolkit API.
 */
export const CISTROME_API_COLORS = Object.freeze({
    INTERVAL: "#2C77B1",
    GENE: "#D6641E",
    FACTOR: "#2B9F78",
    PEAKSET: "#2B9F78"
});

/*
 * Options for the parameters of Cistrome DB APIs.
 */
export const CISTROME_DBTOOLKIT_CHROMOSOMES = [
    "chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10",
    "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", 
    "chr21", "chr22", "chrX", "chrY",
];
export const CISTROME_DBTOOLKIT_SPECIES = ["hg38", "mm10"];
export const CISTROME_DBTOOLKIT_GENE_DISTANCE = ["1kb", "10kb", "100kb"];
export const CISTROME_DBTOOLKIT_PEAK_NUMBERS = [
    "Top 1k peaks according to peak enrichment",
    "Top 10k peaks according to peak enrichment",
    "All peaks in each sample"
];

export const CISTROME_DBTOOLKIT_MAX_INTERVAL_SIZE = 2000000;

/**
 * Get a table with readable column names for an API result table.
 * @param {string} apiType Type of Cistrome DB Toolkit API (`CISTROME_API_TYPES`).
 * @returns {[array, array]} An array of column names and an array of JSON object for row information.
 */
export function getReadableTable(apiType, originalRows) {
    const readableColumns = {
        [CISTROME_API_TYPES.INTERVAL]: {
            GSM: "GEO/ENCODE ID",
            DCid: "CistromeDB ID",
            factor: "Factor",
            cellLine: "Cell Line",
            CellType: "Cell Type",
            species: "Species",
            OverlapRatio: "Overlap Ratio",
            OverlapPeakNumber: "Overlap Peak Number",
        },
        [CISTROME_API_TYPES.GENE]: {
            cellLine: "Cell Line",
            RP: "Regulatory Potential",
            DCid: "CistromeDB ID",
            Tissue: "Tissue",
            CellType: "CellType",
            factor: "Factor",
            GSM: "GEO/ENCODE ID",
            species: "Species"
        },
        [CISTROME_API_TYPES.FACTOR]: {
            title: "Title",
            external_id: "GEO/ENCODE ID",
            id: "CistromeDB ID",
            factor: "Factor",
            species: "Species",
            experiment_type: "Experiment Type",
            sample_type: "Sample Type"
        },
        [CISTROME_API_TYPES.PEAKSET]: {
            // TODO: Support Peak Set API
        }
    };
    if(!Object.values(CISTROME_API_TYPES).includes(apiType) || !readableColumns[apiType]) {
        console.warn("Unsupported Cistrome API");
        return [];
    }
    const columns = Object.values(readableColumns[apiType]);
    const rows = originalRows.map(d => {
        // Change key names using the readable columns
        const newRow = {};
        Object.keys(readableColumns[apiType]).forEach(c => {
            if(readableColumns[apiType][c]) {
                newRow[readableColumns[apiType][c]] = d[c];
            }
        });
        return newRow;
    });
    return [columns, rows];
}

/**
 * Check if interval request parameters are valid, and will be able to generate a good API url.
 * @param {string} assembly
 * @param {string} chrStartName
 * @param {number} chrStartPos
 * @param {string} chrEndName
 * @param {number} chrEndPos
 * @returns {object} Status object, for example `{ msg: "Some error", success: false }`.
 */
export function validateIntervalParams({ assembly, chrStartName, chrStartPos, chrEndName, chrEndPos }) {
    let msg;
    let success = false;
    if(!assembly) {
        msg = "No assembly value found.";
    } else if (!CISTROME_DBTOOLKIT_SPECIES.includes(assembly)) {
        msg = "Unsupported assembly encountered.";
    } else if(!chrStartName || !chrStartPos || !chrEndName || !chrEndPos) {
        msg = "Not enough chromosome parameters.";
    } else if(chrStartName !== chrEndName) {
        msg = "Interval spans across more than one chromosome.";
    } else if(!+chrStartPos || !+chrEndPos) {
        msg = "Genomic position is not a numeric value";
    } else if(+chrStartPos < 0 || +chrEndPos < 0) {
        msg = "Genomic position cannot be a negative value";
    } else if(chrEndPos - chrStartPos > CISTROME_DBTOOLKIT_MAX_INTERVAL_SIZE) {
        msg = "Search requires interval < 2 Mb";
    } else {
        msg = "Success";
        success = true;
    }
    return { msg, success };
}

/**
 * Check if gene request parameters are valid, and will be able to generate a good API url.
 * @param {string} assembly
 * @param {string} gene
 * @param {number} distance
 * @returns {object} Status object, for example `{ msg: "Some error", success: false }`.
 */
export function validateGeneParams({ assembly, gene, distance }) {
    let msg;
    let success = false;
    if(!assembly) {
        msg = "No assembly value found.";
    } else if (!CISTROME_DBTOOLKIT_SPECIES.includes(assembly)) {
        msg = "Unsupported assembly encountered.";
    } else if(!gene || gene === "") {
        msg = "No gene name suggested.";
    } else if(!distance) {
        msg = "No distance suggested.";
    } else if(!CISTROME_DBTOOLKIT_GENE_DISTANCE.includes(distance)) {
        msg = "Unsupported distance.";
    } else {
        msg = "Success";
        success = true;
    }
    return { msg, success };
}

export function validateFactorParams({ factor }) {
    let msg;
    let success = false;
    if(!factor) {
        msg = "Factor is not suggested";
    } else {
        msg = "Success";
        success = true;
    }
    return { msg, success };
}

/**
 * Check if peak set request parameters are valid, and will be able to generate a good API url.
 * @param {string} assembly
 * @param {string} gene
 * @param {number} distance
 * @returns {object} Status object, for example `{ msg: "Some error", success: false }`.
 */
export function validatePeaksetParams({ assembly, tpeak, bedFile }) {
    let msg;
    let success = false;
    if(!assembly) {
        msg = "No assembly value found.";
    } else if (!CISTROME_DBTOOLKIT_SPECIES.includes(assembly)) {
        msg = "Unsupported assembly encountered.";
    } else if(!CISTROME_DBTOOLKIT_PEAK_NUMBERS.includes(tpeak)) {
        msg = "Unsupported Peak Number encounered.";
    } else if(!bedFile) {
        msg = "No bed file suggested.";
    } else {
        msg = "Success";
        success = true;
    }
    return { msg, success };
}

/**
 * Generate a URL for the cistrome DB toolkit API for interval search.
 * @param {string} assembly
 * @param {string} chrStartName
 * @param {number} chrStartPos
 * @param {string} chrEndName
 * @param {number} chrEndPos
 * @returns {string} The url.
 */
export function makeDBToolkitIntervalAPIURL(assembly, chrStartName, chrStartPos, chrEndName, chrEndPos) {
    return `http://dbtoolkit.cistrome.org/api_interval?species=${assembly}&factor=tf&interval=${chrStartName}:${chrStartPos}-${chrEndPos}`;
}

/**
 * Generate a URL for the cistrome DB toolkit API gene search.
 * @param {string} assembly
 * @param {string} gene
 * @param {string} distance
 * @returns {string} The url.
 */
export function makeDBToolkitGeneAPIURL(assembly, gene, distance) {
    // In API, 'b' is not included in the unit for distance, so here we are removing 'b' (e.g., 10kb => 10k)
    return `http://dbtoolkit.cistrome.org/api_gene?species=${assembly}&factor=tf&transcript=${gene}&distance=${distance.replace("b", "")}`;
}

/**
 * Make an API request for Cistrome DB Toolkit based on the API type.
 * @param {string} apiType Type of Cistrome DB Toolkit API (`CISTROME_API_TYPES`).
 * @param {Object} parameter Parameter required for API request.
 */
export function requestDBToolkitAPI(apiType, parameter) {
    return {
        [CISTROME_API_TYPES.INTERVAL]: requestByInterval,
        [CISTROME_API_TYPES.GENE]: requestByGene,
        [CISTROME_API_TYPES.FACTOR]: requestByFactor,
        [CISTROME_API_TYPES.PEAKSET]: requestByPeakset
    }[apiType](parameter);
}

/**
 * Make an API request for Cistrome DB Toolkit interval information.
 * @param {string} assembly
 * @param {string} chrStartName
 * @param {number} chrStartPos
 * @param {string} chrEndName
 * @param {number} chrEndPos
 * @returns {Promise} On success, promise resolves with the following array: `[rows, columns]`.
 */
export function requestByInterval({ assembly, chrStartName, chrStartPos, chrEndName, chrEndPos }) {
    const { 
        msg: validationMsg, 
        success: validationSuccess
    } = validateIntervalParams({ assembly, chrStartName, chrStartPos, chrEndName, chrEndPos });

    if(!validationSuccess) {
        return new Promise((resolve, reject) => {
            reject(validationMsg);
        });
    }

    const dbToolkitAPIURL = makeDBToolkitIntervalAPIURL(assembly, chrStartName, chrStartPos, chrEndName, chrEndPos);

    return fetch(dbToolkitAPIURL)
        .then((response) => {
            if (!response.ok) {
                return new Promise((resolve, reject) => {
                    reject(`Error: ${response.statusText}`);
                });
            }
            return response.json();
        })
        .then((data) => {
            const keys = Object.keys(data);

            return new Promise((resolve, reject) => {
                if(keys.length === 0) {
                    reject(`No data found for interval ${chrStartName}:${chrStartPos}-${chrEndPos}`);
                }
                // Generate data for table.
                const rows = keys.map(k => data[k]);
                const filteredRows = rows.slice(0, rows.length < SEARCH_ITEM_LIMIT ? rows.length : SEARCH_ITEM_LIMIT);
                const columns = Object.keys(data[keys[0]]);
                resolve([filteredRows, columns]);
            });
        })
        .catch(error => {
            return new Promise((resolve, reject) => {
                reject(`Error: ${error.message}`);
            });
        });
}

/**
 * Make an API request for Cistrome DB Toolkit gene information.
 * @param {string} assembly
 * @param {string} gene
 * @param {string} distance
 * @returns {Promise} On success, promise resolves with the following array: `[rows, columns]`.
 */
export function requestByGene({ assembly, gene, distance }) {
    const url = makeDBToolkitGeneAPIURL(assembly, gene, distance);
    return fetch(url)
        .then((response) => {
            if (!response.ok) {
                return new Promise((resolve, reject) => {
                    reject(`Error: ${response.statusText}`);
                });
            }
            return response.json();
        })
        .then((data) => {
            const keys = Object.keys(data);

            return new Promise((resolve, reject) => {
                if(keys.length === 0) {
                    reject(`No data found for gene ${gene}`);
                }
                // Generate data for table.
                const rows = keys.map(k => data[k]);
                const filteredRows = rows.slice(0, rows.length < SEARCH_ITEM_LIMIT ? rows.length : SEARCH_ITEM_LIMIT);
                const columns = Object.keys(data[keys[0]]);
                resolve([filteredRows, columns]);
            });
        })
        .catch(error => {
            return new Promise((resolve, reject) => {
                reject(`Error: ${error.message}`);
            });
        });
}

/**
 * Make an API request for Cistrome DB factor information.
 * @param {string} factor
 * @returns {Promise} On success, promise resolves with the following array: `[rows, columns]`.
 */
export async function requestByFactor({ factor }) {
    const factorIds = await getFactorList();
    const url = `http://develop.cistrome.org/cistrome/factors/${factorIds[factor]}?fields=samples&format=json`;
    return fetch(url)
        .then((response) => {
            if (!response.ok) {
                return new Promise((resolve, reject) => {
                    reject(`Error: ${response.statusText}`);
                });
            }
            return response.json();
        })
        .then((data) => {
            return new Promise((resolve, reject) => {
                const rows = data.factors[0]?.samples;
                if(rows.length === 0) {
                    reject(`No data found for factor ${factor}`);
                }
                // Generate data for table.
                const filteredRows = rows.filter(d => d.species !== 'mus_musculus').slice(0, rows.length < SEARCH_ITEM_LIMIT ? rows.length : SEARCH_ITEM_LIMIT);
                filteredRows.forEach(d => d['factor'] = factor);
                const columns = Object.keys(rows[0]);
                resolve([filteredRows, columns]);
            });
        })
        .catch(error => {
            return new Promise((resolve, reject) => {
                reject(`Error: ${error.message}`);
            });
        });
}

function getFactorList() {
    const factorListUrl = 'http://develop.cistrome.org/cistrome/factors?limit=10000&format=json';
    return fetch(factorListUrl)
        .then((response) => {
            if (!response.ok) {
                return new Promise((resolve, reject) => {
                    reject(`Error: ${response.statusText}`);
                });
            }
            return response.json();
        })
        .then((data) => {
            return new Promise((resolve, reject) => {
                const factors = {};
                data.factors.forEach(d => { factors[d.name] = d.id });
                resolve(factors);
            });
        }).catch(error => {
            return new Promise((resolve, reject) => {
                reject(`Error: ${error.message}`);
            });
        });
}

/**
 * Make an API request for Cistrome DB Toolkit peak set information.
 * @param {string} assembly
 * @param {string} tpeak
 * @param {string} bedFile
 * @returns {Promise} On success, promise resolves with the following array: `[rows, columns]`.
 */
export function requestByPeakset({ assembly, tpeak, bedFile }) {

    // TODO: Get csrftoken from the Toolkit website
    // let csrftoken = Cookies.get('csrftoken');
    // console.log(csrftoken);

    const formData = new FormData();
    formData.append("csrfmiddlewaretoken", "Xy1sb4N0NM7XUvaMejxinyXOSasJX6mS");
    formData.append("species", "hg38");
    formData.append("factor", "tf");
    formData.append("tpeak", "1k");
    formData.append("peak", bedFile);

    return fetch("http://dbtoolkit.cistrome.org/api_similar", {
        credentials: "omit", // https://stackoverflow.com/a/50388440
        method: "post",
        body: formData,
    })
        .then((response) => {
            // console.log('response', response);
            if (!response.ok) {
                return new Promise((resolve, reject) => {
                    reject(`Error: ${response.statusText}`);
                });
            }
            return response.json();
        })
        .then((data) => {
            // console.log('data', data);
            const keys = Object.keys(data);

            return new Promise((resolve, reject) => {
                if(keys.length === 0) {
                    reject("No data found for the given bed file");
                }
                // Generate data for table.
                const rows = keys.map(k => data[k]);
                const filteredRows = rows.slice(0, rows.length < SEARCH_ITEM_LIMIT ? rows.length : SEARCH_ITEM_LIMIT);
                const columns = Object.keys(data[keys[0]]);
                resolve([filteredRows, columns]);
            });
        })
        .catch(error => {
            // console.log('error', error);
            return new Promise((resolve, reject) => {
                reject(`Error: ${error.message}`);
            });
        });
}

export function getFactor(id) {
    return fetch(`http://develop.cistrome.org/cistrome/samples/${id}?fields=factors&format=json`)
        .then(response => {
            if (!response.ok) {
                return new Promise((resolve, reject) => {
                    reject(`Error: ${response.statusText}`);
                });
            }
            return response.json();
        }).then(data => {
            return new Promise((resolve, reject) => {
                const factors = data.samples[0]?.factors;
                // if(factors.length === 0) {
                    // reject("No factor data found for given sample");
                // }
                resolve(factors.sort((a, b) => a.value - b.value)[0]?.name ?? 'Unknown')
            });
        }).catch(error => {
            return new Promise((resolve, reject) => {
                reject(`Error: ${error.message}`);
            });
        });
}

export function getOntology(id) {
    return fetch(`http://develop.cistrome.org/cistrome/samples/${id}?fields=ontologies&format=json`)
        .then(response => {
            if (!response.ok) {
                return new Promise((resolve, reject) => {
                    reject(`Error: ${response.statusText}`);
                });
            }
            return response.json();
        }).then(data => {
            return new Promise((resolve, reject) => {
                const ontologies = data.samples[0]?.ontologies;
                // if(ontologies.length === 0) {
                    // reject("No ontology data found for given sample");
                // }
                resolve(ontologies.sort((a, b) => a.value - b.value)[0]?.term ?? 'Unknown')
            });
        }).catch(error => {
            return new Promise((resolve, reject) => {
                reject(`Error: ${error.message}`);
            });
        });
}
