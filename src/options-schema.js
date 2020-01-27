import merge from 'lodash/merge';
import cloneDeep from 'lodash/cloneDeep';
import Ajv from 'ajv';

const DEFAULT_KEY = "__default__";

const baseSchema = {
    "$schema": "http://json-schema.org/draft-07/schema#",
    "definitions": {
        "trackOptions": {
            "type": "object",
            "additionalProperties": false,
            "properties": {
                "uid": {
                    "type": "string"
                },
                "rowInfoPosition": {
                    "type": "string",
                    "enum": ["hidden", "right", "left"]
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
            "required": ["uid"]
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
                "uid": {
                    "default": DEFAULT_KEY
                }
            }
        }
    }
});

/**
 * @param {any} optionsRaw The raw value of the options prop.
 * @returns {object} The processed wOptions object, mapping track IDs to options objects.
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

export function processWrapperOptions(optionsRaw) {
    // Set up the default options:
    const optionsProcessed = {
        default: {
            rowInfoPosition: "right",
        }
    };

    // Validate the raw options:
    const valid = validateWrapperOptions(optionsRaw);
    if(!valid) {
        return optionsProcessed;
    }

    // Process the raw options by merging into the processed options object:
    if(Array.isArray(optionsRaw)) {
        const optionsDefault = optionsRaw.find(o => (o.uid === DEFAULT_KEY));
        if(optionsDefault) {
            merge(optionsProcessed.default, optionsDefault);
        }
        optionsRaw.forEach((trackOptions) => {
            if(trackOptions.uid !== DEFAULT_KEY) {
                optionsProcessed[trackOptions.uid] = merge(cloneDeep(trackOptions), optionsProcessed.default);
            }
        });
    } else if(typeof optionsRaw === "object") {
        merge(optionsProcessed.default, optionsRaw);
    }

    return optionsProcessed;
}