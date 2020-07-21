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
                "rowAggregate": {
                    "type": "array",
                    "items": { "$ref": "#/definitions/aggregateInfo" }
                },
                "rowHighlight": {
                    "type": "object",
                    "oneOf":[ 
                        { "required": ["field", "type", "contains"] },
                        { "required": ["field", "type", "range"] },
                        { "required": ["field", "type", "subtree"] },
                        { "required": ["field", "type", "minSimilarity"] },
                    ],
                    "properties": {
                        "field": {
                            "type": "string",
                            "description": "The name of a data field"
                        },
                        "type": {
                            "type": "string",
                            "enum": ["nominal", "quantitative", "tree"],
                            "description": "The data type of a field"
                        },
                        "contains": {
                            "type": "string",
                            "description": "The substring to search for"
                        },
                        "range": {
                            "type": "array",
                            "description": "Min and max values"
                        },
                        "subtree": {
                            "type": "array",
                            "description": "The subtree to search for"
                        },
                        "minSimilarity": {
                            "type": "number",
                            "description": "A similarity threshold"
                        },
                    }
                },
                "rowZoom": { "$ref": "#/definitions/zoomInfo" }
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
                "aggFunction": {
                    "type": "string",
                    "enum": ["max", "min", "mean", "sum", "mostCommon", "leastCommon", "concat", "count", "uniqueCount"],
                    "description": "A funtion to apply when aggregating values"
                },
                "position": {
                    "type": "string",
                    "enum": ["left", "right"],
                    "description": "The position to show a data attribute relative to a higlass track"
                },
                "title": {
                    "type": "string",
                    "description": "The name of a data field to alternatively use for displaying urls"
                },
                "resolveYScale": {
                    "type": "boolean",
                    "description": "Determine if the scale of y axis should be independent to the adjacently placed tracks"
                },
                "sort": {
                    "type": "string",
                    "enum": ["descending", "ascending"],
                    "description": "The order of sorting. This only works if `resolveYScale` is set to `true`"
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
                { "required": ["field", "type", "range"] },
                { "required": ["field", "type", "subtree"] },
                { "required": ["field", "type", "subtree", "minSimilarity"] },
                { "required": ["field", "type", "minSimilarity"] },
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
                    "type": "array",
                    "description": "An array of values of a field that should not be included in the filtered data"
                },
                "range": {
                    "type": "array",
                    "description": "Min and max values"
                },
                "subtree": {
                    "type": "array",
                    "description": "The subtree to search for"
                },
                "minSimilarity": {
                    "type": "number",
                    "description": "A similarity threshold"
                },
            }
        },
        "aggregateInfo": {
            "type": "object",
            "required": ["field", "type"],
            "properties": {
                "field": {
                    "type": "string",
                    "description": "The name of a data field"
                },
                "type": {
                    "type": "string",
                    "enum": ["nominal"], // TODO: Support other types as well.
                    "description": "The data type of a field"
                },
                "oneOf": {
                    "type": "array",
                    "description": "An array of values of a field that should be aggregated"
                }
            }
        },
        "zoomInfo": {
            "type": "object",
            "required": ["level", "top", "numRows"],
            "properties": {
                "level": {
                    "type": "number",
                    "description": "The fraction of rows visible in the zoom window. Between 0 and 1."
                },
                "top": {
                    "type": "number",
                    "description": "The top offset of the zoom window. Between 0 and 1."
                },
                "numRows": {
                    "type": "number",
                    "description": "The number of rows to take into account (after filtering but before zooming)."
                }
            }
        }
    }
};

const optionsArraySchema = merge(cloneDeep(baseSchema), {
    "title": "HiGlassWithMetadata wOptions (array)",
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
    "title": "HiGlassWithMetadata wOptions (object)",
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
 * Get a key in `rowHighlight` object that indicate a certain condition (e.g., `contains`) for highlighting.
 * @param {string} type The field type.
 * @returns {string} The key of `rowHighlight` object that indicate a certain condition.
 */
export function getHighlightKeyByFieldType(type, condition) {
  switch(type) {
    case "quantitative":
        return "range";
    case "nominal":
        return "contains";
    case "tree":
        return Array.isArray(condition) ? "subtree" : "minSimilarity";
  }
}

/**
 * Validate the HiGlassWithMetadata `options` prop.
 * @param {(object|object[]|null)} options The raw value of the options prop.
 * @returns {boolean} True if the options prop value was valid.
 */
export function validateWrapperOptions(options) {
    let validate;
    if(Array.isArray(options)) {
        validate = new Ajv({ extendRefs: true }).compile(optionsArraySchema);
    } else if(typeof options === "object") {
        validate = new Ajv({ extendRefs: true }).compile(optionsObjectSchema);
    } else if(!options) {
        return true;
    }
    const valid = validate(options);

    if (validate.errors) {
        console.warn(JSON.stringify(validate.errors, null, 2));
    }

    return valid;
}

export function isProcessedWrapperOptions(options) {
    // TODO: We can make a json schema for the processed options
    // to more thoroughly check if this `options` is processed one.
    
    // `DEFAULT_OPTIONS_KEY` is always included in the processed wrapper options, but not in raw options.
    return Object.keys(options).indexOf(DEFAULT_OPTIONS_KEY) !== -1;
}
/**
 * Process the HiGlassMeta `options` prop by mapping track IDs to objects containing values for all possible option attributes.
 * @param {(object|object[]|null)} options The raw value of the options prop.
 * @returns {object} A processed options object, mapping track IDs to options objects, and merging with defaults.
 */
export function processWrapperOptions(options) {
    if(isProcessedWrapperOptions(options)) {
        // This means we do not have to process the options.
        return options;
    }

    // Set up the default options:
    const newOptions = {
        [DEFAULT_OPTIONS_KEY]: {
            rowInfoAttributes: [],
            rowFilter: [],
            rowSort: [],
            rowHighlight: {},
            rowZoom: {},
        }
    };

    // Validate the raw options:
    const valid = validateWrapperOptions(options);
    if(!valid) {
        console.warn("Invalid Wrapper Options in processWrapperOptions().");
        return newOptions;
    }

    // Process the raw options by merging into the processed options object:
    if(Array.isArray(options)) {
        // Check for global options.
        const globalDefaults = options.find(o => (o.viewId === DEFAULT_OPTIONS_KEY && o.trackId === DEFAULT_OPTIONS_KEY));
        if(globalDefaults) {
            newOptions[DEFAULT_OPTIONS_KEY] = merge(
                cloneDeep(newOptions[DEFAULT_OPTIONS_KEY]), 
                omit(globalDefaults, ['viewId', 'trackId'])
            );
        }

        // Check for view-specific, track-global options.
        options.forEach((trackOptions) => {
            if(trackOptions.viewId !== DEFAULT_OPTIONS_KEY && trackOptions.trackId === DEFAULT_OPTIONS_KEY) {
                newOptions[trackOptions.viewId] = {
                    [DEFAULT_OPTIONS_KEY]: merge(
                        cloneDeep(newOptions[DEFAULT_OPTIONS_KEY]), 
                        omit(trackOptions, ['viewId', 'trackId'])
                    )
                };
            }
        });

        // Check for view-specific and track-specific options.
        options.forEach((trackOptions) => {
            if(trackOptions.viewId !== DEFAULT_OPTIONS_KEY && trackOptions.trackId !== DEFAULT_OPTIONS_KEY) {
                if(!newOptions[trackOptions.viewId]) {
                    newOptions[trackOptions.viewId] = {
                        [DEFAULT_OPTIONS_KEY]: cloneDeep(newOptions[DEFAULT_OPTIONS_KEY])
                    };
                }
                newOptions[trackOptions.viewId][trackOptions.trackId] = merge(
                    cloneDeep(newOptions[trackOptions.viewId][DEFAULT_OPTIONS_KEY]), 
                    omit(trackOptions, ['viewId', 'trackId'])
                );
            }
        });
    } else if(typeof options === "object") {
        newOptions[DEFAULT_OPTIONS_KEY] = merge(cloneDeep(newOptions[DEFAULT_OPTIONS_KEY]), options);
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
