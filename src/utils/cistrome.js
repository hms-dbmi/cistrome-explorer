
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
export function validateIntervalParams(assembly, chrStartName, chrStartPos, chrEndName, chrEndPos) {
    let msg;
    let success = false;
    if(!assembly) {
        msg = "No assembly value found.";
    } else if (!["hg38", "mm10"].includes(assembly)) {
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
 *  Generate a URL for the cistrome DB toolkit API.
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
    } = validateIntervalParams(assembly, chrStartName, chrStartPos, chrEndName, chrEndPos);

    if(!validationSuccess) {
        return new Promise((resolve, reject) => {
            reject(validationMsg);
        })
    }

    const dbToolkitAPIURL = makeDBToolkitIntervalAPIURL(assembly, chrStartName, chrStartPos, chrEndName, chrEndPos);

    return fetch(dbToolkitAPIURL, { cache: "force-cache" })
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