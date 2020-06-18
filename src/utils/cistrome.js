export const CISTROME_DBTOOLKIT_SPECIES = ["hg38", "mm10"];
export const CISTROME_DBTOOLKIT_GENE_DISTANCE = ['1kb', '10kb', '100kb'];
export const CISTROME_DBTOOLKIT_PEAK_NUMBERS = [
    'Top 1k peaks according to peak enrichment',
    'Top 10k peaks according to peak enrichment',
    'All peaks in each sample'
];
export const CISTROME_DBTOOLKIT_MAX_INTERVAL_SIZE = 2000000;

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
    } else if(!gene || gene === '') {
        msg = "No gene name suggested.";
    } else if(!distance) {
        msg = "No distance suggested.";
    } else if(!CISTROME_DBTOOLKIT_GENE_DISTANCE.includes(distance)) {
        msg = "Unsupported distance.";
    } else {
        msg = "Success";
        success = true;
    }
    console.log({msg, success});
    return { msg, success };
}

/**
 *  Generate a URL for the cistrome DB toolkit API for interval search.
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
 *  Generate a URL for the cistrome DB toolkit API gene search.
 * @param {string} assembly
 * @param {string} gene
 * @param {string} distance
 * @returns {string} The url.
 */
export function makeDBToolkitGeneAPIURL(assembly, gene, distance) {
    // In API, 'b' is not included in the unit for distance, so here we are removing 'b' (e.g., 10kb => 10k)
    return `http://dbtoolkit.cistrome.org/api_gene?species=${assembly}&factor=tf&transcript=${gene}&distance=${distance.replace('b', '')}`;
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
export function requestIntervalTFs(assembly, chrStartName, chrStartPos, chrEndName, chrEndPos) {
    const { 
        msg: validationMsg, 
        success: validationSuccess
    } = validateIntervalParams({ assembly, chrStartName, chrStartPos, chrEndName, chrEndPos });

    if(!validationSuccess) {
        return new Promise((resolve, reject) => {
            reject(validationMsg);
        })
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
                const filteredRows = rows.slice(0, rows.length < 100 ? rows.length : 100);
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
export function requestGeneTFs({ assembly, gene, distance }) {
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
                const filteredRows = rows.slice(0, rows.length < 100 ? rows.length : 100);
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