import merge from 'lodash/merge';
import cloneDeep from 'lodash/cloneDeep';
import omit from 'lodash/omit';
import Ajv from 'ajv';
import { insertItemToArray } from './array.js'

export const DEFAULT_OPTIONS_KEY = "default";

const baseSchema = {
    "$schema": "http://json-schema.org/draft-07/schema#",
    "definitions": {
        "trackOptions": {
            "type": "object",
            "additionalProperties": false,
            "properties": {
                "viewId": { "type": "string" },
                "trackId": { "type": "string" },
                "rowInfoAttributes": {
                    "type": "array",
                    "items": { "$ref": "#/definitions/fieldInfo" }
                },
                "rowSort": {
                    "type": "array",
                    "items": { "$ref": "#/definitions/sortInfo" }
                },
                "rowFilter": {
                    "type": "array",
                    "items": { "$ref": "#/definitions/filterInfo" }
                },
                "rowHighlight": {
                    "type": "object",
                    "oneOf":[ 
                        { "required": ["field", "type", "contains"] },
                        { "required": ["field", "type", "range"] }
                    ],
                    "properties": {
                        "field": {
                            "type": "string",
                            "description": "The name of a data field"
                        },
                        "type": {
                            "type": "string",
                            "enum": ["nominal", "quantitative"],
                            "description": "The data type of a field"
                        },
                        "contains": {
                            "type": "string",
                            "description": "The substring to search for"
                        },
                        "range": {
                            "type": "array",
                            "description": "Min and max values"
                        }
                    }
                }
            }
        },
        "fieldInfo": {
            "type": "object",
            "required": ["field", "type"],
            "properties": {
                "field": {
                    "type": ["string", "array"],
                    "description": "The data field name(s)"
                },
                "type": {
                    "type": "string",
                    "enum": ["nominal", "quantitative", "url", "tree"],
                    "description": "The data type of a field"
                },
                "position": {
                    "type": "string",
                    "enum": ["left", "right"],
                    "description": "The position to show a data attribute relative to a higlass track"
                },
                "title": {
                    "type": "string",
                    "description": "The name of a data field to alternatively use for displaying urls"
                }
            }
        },
        "sortInfo": {
            "type": "object",
            "required": ["field", "order"],
            "properties": {
                "field": {
                    "type": ["string", "array"],
                    "description": "The name of a data field"
                },
                "type": {
                    "type": "string",
                    "enum": ["nominal", "quantitative", "url", "tree"],
                    "description": "The data type of a field"
                },
                "order": {
                    "type": "string",
                    "enum": ["descending", "ascending"],
                    "description": "The order of sorting"
                }
            }
        },
        "filterInfo": {
            "type": "object",
            "oneOf":[ 
                { "required": ["field", "type", "notOneOf"] },
                { "required": ["field", "type", "subtree"] },
                { "required": ["field", "type", "range"] }
            ],
            "properties": {
                "field": {
                    "type": ["string", "array"],
                    "description": "The name of a data field"
                },
                "type": {
                    "type": "string",
                    "enum": ["nominal", "quantitative", "tree"],
                    "description": "The data type of a field"
                },
                "notOneOf": {
                    // TODO: Add tests for this.
                    "type": "array",
                    "description": "A string of the value of the field should not be included in the filtered data"
                },
                "subtree": {
                    "type": "array",
                    "description": "The substring to search for"
                },
                "range": {
                    "type": "array",
                    "description": "Min and max values"
                }
            }
        }
    }
};

const optionsArraySchema = merge(cloneDeep(baseSchema), {
    "title": "CistromeHGW wOptions (array)",
    "type": "array",
    "items": {
        "type": "object",
        "$ref": "#/definitions/trackOptions"
    },
    "definitions": {
        "trackOptions": {
            "required": ["viewId", "trackId"]
        }
    }
});

const optionsObjectSchema = merge(cloneDeep(baseSchema), {
    "title": "CistromeHGW wOptions (object)",
    "type": "object",
    "$ref": "#/definitions/trackOptions",
    "definitions": {
        "trackOptions": {
            "properties": {
                "viewId": {
                    "default": DEFAULT_OPTIONS_KEY
                },
                "trackId": {
                    "default": DEFAULT_OPTIONS_KEY
                }
            }
        }
    }
});

/**
 * Validate the CistromeHGW `options` prop.
 * @param {(object|object[]|null)} optionsRaw The raw value of the options prop.
 * @returns {boolean} True if the options prop value was valid.
 */
export function validateWrapperOptions(optionsRaw) {
    let validate;
    if(Array.isArray(optionsRaw)) {
        validate = new Ajv({ extendRefs: true }).compile(optionsArraySchema);
    } else if(typeof optionsRaw === "object") {
        validate = new Ajv({ extendRefs: true }).compile(optionsObjectSchema);
    } else if(!optionsRaw) {
        return true;
    }
    const valid = validate(optionsRaw);

    if (validate.errors) {
        console.warn(JSON.stringify(validate.errors, null, 2));
    }

    return valid;
}

/**
 * Process the CistromeHGW `options` prop by mapping track IDs to objects containing values for all possible option attributes.
 * @param {(object|object[]|null)} optionsRaw The raw value of the options prop.
 * @returns {object} A processed options object, mapping track IDs to options objects, and merging with defaults.
 */
export function processWrapperOptions(optionsRaw) {

    // Important Descriptions about Wrapper Options (i.e., optionsRaw):
    //  * Single 'view' can contain multiple 'tracks,' but not vice versa.
    //  * Unlike HiGlass View Config, individual options for each combination of {viewId, trackId} are stored in an 1D array of JSON objects,
    //    instead of using a nested format.
    //  * Both the viewId and trackId or only a trackId can be DEAULT_OPTIONS_KEY (i.e., only the viewId cannot be DEFAULT_OPTIONS_KEY).
    //  * An option of both viewId and trackId being DEFAULT_OPTIONS_KEY is a global option, which affects to any other tracks in any views.
    //  * An option of only a trackId being DEFAULT_OPTIONS_KEY is a track-global option, which affects to any tracks in a certain view.
    //  * Wrapper options may or may not contain options for all of the individual {viewId, trackId} combinations.

    // Set up the default options:
    const options = {
        [DEFAULT_OPTIONS_KEY]: {
            rowInfoAttributes: []
        }
    };

    // Validate the raw options:
    const valid = validateWrapperOptions(optionsRaw);
    if(!valid) {
        console.warn("Invalid Wrapper Options in processWrapperOptions().");
        return options;
    }

    // Process the raw options by merging into the processed options object:
    if(Array.isArray(optionsRaw)) {
        // Check for global options.
        const globalDefaults = optionsRaw.find(o => (o.viewId === DEFAULT_OPTIONS_KEY && o.trackId === DEFAULT_OPTIONS_KEY));
        if(globalDefaults) {
            options[DEFAULT_OPTIONS_KEY] = merge(
                cloneDeep(options[DEFAULT_OPTIONS_KEY]), 
                omit(globalDefaults, ['viewId', 'trackId'])
            );
        }

        // Check for view-specific, track-global options.
        optionsRaw.forEach((trackOptions) => {
            if(trackOptions.viewId !== DEFAULT_OPTIONS_KEY && trackOptions.trackId === DEFAULT_OPTIONS_KEY) {
                options[trackOptions.viewId] = {
                    [DEFAULT_OPTIONS_KEY]: merge(
                        cloneDeep(options[DEFAULT_OPTIONS_KEY]), 
                        omit(trackOptions, ['viewId', 'trackId'])
                    )
                };
            }
        });

        // Check for view-specific and track-specific options.
        optionsRaw.forEach((trackOptions) => {
            if(trackOptions.viewId !== DEFAULT_OPTIONS_KEY && trackOptions.trackId !== DEFAULT_OPTIONS_KEY) {
                if(!options[trackOptions.viewId]) {
                    options[trackOptions.viewId] = {
                        [DEFAULT_OPTIONS_KEY]: cloneDeep(options[DEFAULT_OPTIONS_KEY])
                    };
                }
                options[trackOptions.viewId][trackOptions.trackId] = merge(
                    cloneDeep(options[trackOptions.viewId][DEFAULT_OPTIONS_KEY]), 
                    omit(trackOptions, ['viewId', 'trackId'])
                );
            }
        });
    } else if(typeof optionsRaw === "object") {
        options[DEFAULT_OPTIONS_KEY] = merge(cloneDeep(options[DEFAULT_OPTIONS_KEY]), optionsRaw);
    }

    return options;
}

/**
 * Get the options for a specific track, using its viewId and trackId.
 * @param {object} options A _processed_ options object.
 * @param {string} viewId The viewId for the track of interest.
 * @param {string} trackId The trackId for the track of interest.
 * @returns {object} The options object for the track, or the default options object.
 */
export function getTrackWrapperOptions(options, viewId, trackId) {
    if(options[viewId]) {
        if(options[viewId][trackId]) {
            return options[viewId][trackId];
        } else if(options[viewId][DEFAULT_OPTIONS_KEY]) {
            return options[viewId][DEFAULT_OPTIONS_KEY];
        }
    }
    return options[DEFAULT_OPTIONS_KEY];
}

/**
 * Add options for a specific track, using its viewId and trackId.
 * @param {object} options A _processed_ options object.
 * @param {object} optionsToAdd Options for a specific track.
 * @param {string} viewId The viewId for the track of interest.
 * @param {string} trackId The trackId for the track of interest.
 */
export function addTrackWrapperOptions(options, optionsToAdd, viewId, trackId) {
    if(!options[viewId]) {
        return {
            ...options,
            [viewId]: { 
                [trackId]: optionsToAdd 
            }
        }
    } else {
        return {
            ...options,
            [viewId]: {
                ...options[viewId], 
                [trackId]: optionsToAdd
            }
        }
    }
}

/**
 * Get sub-options for a specific track, using its viewId and trackId.
 * @param {object} options A _processed_ options object to update.
 * @param {string} key The key of sub-options object.
 * @param {string} viewId The viewId for the track of interest.
 * @param {string} trackId The trackId for the track of interest.
 * @returns {object} The options object for the track, or the default options object.
 */
export function getWrapperSubOptions(options, key, viewId, trackId) {
    if(!options[viewId] || (options[viewId] && !options[viewId][trackId])) {
        // Update global defaults if there is no track specific options.
        return options[DEFAULT_OPTIONS_KEY][key];
    } else {
        // Update track specific options.
        return options[viewId][trackId][key];
    }
}

/**
 * Update options for a specific track, using its viewId and trackId.
 * @param {object} options A _processed_ options object to update.
 * @param {object} subOptions A sub-options object to replace or insert.
 * @param {string} key The key of sub-options object.
 * @param {string} viewId The viewId for the track of interest.
 * @param {string} trackId The trackId for the track of interest.
 * @param {boolean} isReplace Determine to replace or insert the subOptions.
 * @returns {object} The options object for the track, or the default options object.
 */
export function updateWrapperOptions(options, subOptions, key, viewId, trackId, { isReplace }) {
    if(!options[viewId] || (options[viewId] && !options[viewId][trackId])) {
        // Update global defaults if there is no track specific options.
        return {
            ...options,
            [DEFAULT_OPTIONS_KEY]: {
                ...options[DEFAULT_OPTIONS_KEY],
                [key]: isReplace ? subOptions : insertItemToArray(options[DEFAULT_OPTIONS_KEY][key], 0, subOptions)
                        
            }
        };
    } else {
        // Update track specific options.
        return {
            ...options,
            [viewId]: {
                ...options[viewId],
                [trackId]: {
                    ...options[viewId][trackId],
                    [key]: isReplace ? subOptions : insertItemToArray(options[viewId][trackId][key], 0, subOptions)
                }
            }
        };
    }
}