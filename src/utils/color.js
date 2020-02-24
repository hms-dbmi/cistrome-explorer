/**
 * Change number (0 ~ 255) to hex string.
 * @prop {number} c Color intensity ranges from 0 to 255.
 * @returns {string} The hex string of a given intensity.
 */
export function componentToHex(c) {
    var hex = c.toString(16);
    return hex.length == 1 ? "0" + hex : hex;
}

/**
 * A function to convert from rgb to hex color.
 * @prop {number} r Red color intensity (0 ~ 255).
 * @prop {number} g Green color intensity (0 ~ 255).
 * @prop {number} b Blue color intensity (0 ~ 255).
 * @returns {string} Color in a hex string format.
 */
export function rgbToHex([r, g, b]) {
    return "#" + componentToHex(r) + componentToHex(g) + componentToHex(b);
}

/**
 * Get unique color with a number
 * @prop {number} i Index of a unique color.
 * @returns {string} The hex string of unique color
 */
export function generateNextUniqueColor(i) {
    if(i < 16777215) {
        // https://stackoverflow.com/questions/15804149/rgb-color-permutation/15804183#15804183
        return rgbToHex([i & 0xff, (i & 0xff00) >> 8, (i & 0xff0000) >> 16]);
    } else {
        console.log("WARNING: unique colors out of range.");
        return "#000000";
    }
}