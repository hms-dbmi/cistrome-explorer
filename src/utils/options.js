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
                }
            }
        },
        "fieldInfo": {
            "type": "object",
            "required": ["field", "type"],
            "properties": {
                "field": {
                    "type": "string",
                    "description": "The name of a data field"
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
                    "type": "string",
                    "description": "The name of a data field"
                },
                "type": {
                    "type": "string",
                    "enum": ["nominal", "quantitative", "url"],
                    "description": "The data type of a field"
                },
                "order": {
                    "type": "string",
                    "enum": ["descending", "ascending"],
                    "description": "The order of sorting"
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
 * TODO
 */
export function getTrackWrapperOptions(options, viewId, trackId) {
    if(options[viewId]) {
        if(options[viewId][trackId]) {
            return options[viewId][trackId];
        } else {
            return options[viewId][DEFAULT_OPTIONS_KEY];
        }
    } else {
        return options[DEFAULT_OPTIONS_KEY];
    }
}

/**
 * Update rowSort information in options.
 * @param {(object|object[]|null)} options The raw value of the options prop.
 * @param {object} sortInfo The name and type of data field and sorting order.
 * @returns {(object|object[])} The new options object or array.
 */
export function updateRowSortOptions(options, sortInfo) {
    let optionsNewSort;
    let newRowSort = [{
        field: sortInfo.field,
        type: sortInfo.type,
        order: sortInfo.order
    }];

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
        optionsNewSort = modifyItemInArray(optionsRaw, index, {
            ...globalDefaults,
            rowSort: newRowSort
        });
    } else {
        optionsNewSort = {
            ...options,
            rowSort: newRowSort
        };
    }
    return optionsNewSort;
}