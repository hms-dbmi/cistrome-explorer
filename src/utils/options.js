import merge from 'lodash/merge';
import cloneDeep from 'lodash/cloneDeep';
import omit from 'lodash/omit';
import Ajv from 'ajv';
import { insertItemToArray, modifyItemInArray } from './array.js'

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
                "colToolsPosition": {
                    "type": "string",
                    "enum": ["hidden", "bottom", "top"]
                },
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
                    "required": ["field", "type", "contains"],
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
            "required": ["field", "type", "contains"],
            "properties": {
                "field": {
                    "type": ["string", "array"],
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
            colToolsPosition: "bottom",
            rowInfoAttributes: []
        }
    };

    // Validate the raw options:
    const valid = validateWrapperOptions(optionsRaw);
    if(!valid) {
        console.log("WARNING: Invalid Wrapper Options.");
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
 * Update a part of options, such as rowSort or rowHighlight.
 * @param {(object|object[]|null)} options The raw value of the options prop.
 * @param {object} subOption New sub-option to replace.
 * @param {string} key Key of the sub-option to replace in "options."
 * @param {boolean} isReplace Replace a suboption with new one or add to array
 * @returns {(object|object[])} The new options object or array.
 */
export function updateGlobalOptionsWithKey(options, subOption, key, { isReplace }) {
    let newOptions;

    if(Array.isArray(options)){
        let optionsRaw = options.slice();
        let globalDefaults = optionsRaw.find(o => (o.viewId === DEFAULT_OPTIONS_KEY && o.trackId === DEFAULT_OPTIONS_KEY));
        
        // If there is no globar defaults, add one.
        if(!globalDefaults) {
            globalDefaults = {
                viewId: DEFAULT_OPTIONS_KEY,
                trackId: DEFAULT_OPTIONS_KEY
            };
            optionsRaw = insertItemToArray(optionsRaw, 0, globalDefaults);
        }

        const index = optionsRaw.indexOf(globalDefaults);
        newOptions = modifyItemInArray(optionsRaw, index, {
            ...globalDefaults,
            [key]: isReplace ? subOption : insertItemToArray(globalDefaults[key], 0, subOption)
        });
    } else {
        newOptions = {
            ...options,
            [key]: isReplace ? subOption : insertItemToArray(options[key], 0, subOption)
        };
    }
    return newOptions;
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
    const _viewId = options[viewId] ? viewId : DEFAULT_OPTIONS_KEY;
    if(_viewId === DEFAULT_OPTIONS_KEY) {
        // Global defaults (1D array).
        return {
            ...options,
            [_viewId]: {
                ...options[_viewId],
                [key]: isReplace ? subOptions : insertItemToArray(options[_viewId][key], 0, subOptions)
            }
        };
    } else {
        // 2D array.
        const _trackId = options[_viewId][trackId] ? trackId : DEFAULT_OPTIONS_KEY;
        return {
            ...options,
            [_viewId]: {
                ...options[_viewId],
                [_trackId]: {
                    ...options[_viewId][_trackId],
                    [key]: isReplace ? subOptions : insertItemToArray(options[_viewId][_trackId][key], 0, subOptions)
                }
            }
        };
    }
}