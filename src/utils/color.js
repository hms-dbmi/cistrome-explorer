export function componentToHex(c) {
    var hex = c.toString(16);
    return hex.length == 1 ? "0" + hex : hex;
}

export function rgbToHex([r, g, b]) {
    return "#" + componentToHex(r) + componentToHex(g) + componentToHex(b);
}

// Unique color generation
// https://stackoverflow.com/questions/15804149/rgb-color-permutation/15804183#15804183
export function generateNextUniqueColor(i) {
    if(i < 16777215) {
        return rgbToHex([i & 0xff, (i & 0xff00) >> 8, (i & 0xff0000) >> 16]);
    } else {
        console.log("WARNING: unique colors out of range.");
        return "white";
    }
}