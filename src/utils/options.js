import merge from 'lodash/merge';
import cloneDeep from 'lodash/cloneDeep';
import omit from 'lodash/omit';
import Ajv from 'ajv';

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
                "rowInfoPosition": {
                    "type": "string",
                    "enum": ["hidden", "left", "right"]
                },
                "rowLinkPosition": {
                    "type": "string",
                    "enum": ["hidden", "left", "right"]
                },
                "rowLinkAttribute": {
                    "type": "string"
                },
                "rowLinkNameAttribute": {
                    "type": "string"
                },
                "colToolsPosition": {
                    "type": "string",
                    "enum": ["hidden", "bottom", "top"]
                },
                "infoAttributes": {
                    "type": "array",
                    "items": { "$ref": "#/definitions/fieldInfo" }
                }
            }
        },
        "fieldInfo": {
            "type": "object",
            "required": ["name", "type"],
            "properties": {
                "name": {
                    "type": "string",
                    "description": "The name of data field"
                },
                "type": {
                    "type": "string",
                    "enum": ["nominal", "quantitative", "tree"],
                    "description": "The data type of a field"
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
            rowInfoPosition: "right",
            rowLinkPosition: "left",
            rowLinkAttribute: "url",
            colToolsPosition: "bottom"
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