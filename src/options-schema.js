import merge from 'lodash/merge';
import cloneDeep from 'lodash/cloneDeep';
import Ajv from 'ajv';

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

export const wOptionsArraySchema = merge(cloneDeep(baseSchema), {
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

export const wOptionsObjectSchema = merge(cloneDeep(baseSchema), {
    "title": "CistromeHGW wOptions (object)",
    "type": "object",
    "$ref": "#/definitions/trackOptions",
    "definitions": {
        "trackOptions": {
            "properties": {
                "uid": {
                    "default": "default"
                }
            }
        }
    }
});

/**
 * @param {any} wOptions The raw value of the wOptions prop.
 * @returns {object} The processed wOptions object, mapping track IDs to options objects.
 */
export function validateWrapperOptions(wOptions) {
    let validate;
    if(Array.isArray(wOptions)) {
        validate = new Ajv({ extendRefs: true }).compile(wOptionsArraySchema);
    } else if(typeof wOptions === "object") {
        validate = new Ajv({ extendRefs: true }).compile(wOptionsObjectSchema);
    } else if(!wOptions) {
        return true;
    }
    const valid = validate(wOptions);

    if (validate.errors) {
        console.warn(JSON.stringify(validate.errors, null, 2));
    }

    return valid;
}

export function processWrapperOptions(wOptions) {

    validateWrapperOptions(wOptions);

    const processedOptions = {
        default: {
            rowInfoPosition: "right",
        }
    };

   if(Array.isArray(wOptions)) {
        wOptions.forEach((trackOptions) => {
            if(trackOptions.uid === "default") {
                merge(processedOptions.default, trackOptions);
            } else {
                processedOptions[trackOptions.uid] = trackOptions;
            }
        });
    } else if(typeof wOptions === "object") {
        merge(processedOptions.default, wOptions);
    }

    return processedOptions;

}