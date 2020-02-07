
/**
 * Get the retina ratio to be able to scale up a canvas context.
 * @param {object} context The canvas context.
 * @returns {number} The ratio.
 */
export function getRetinaRatio(context) {
    let devicePixelRatio = window.devicePixelRatio || 1;
    let backingStoreRatio = [
        context.webkitBackingStorePixelRatio,
        context.mozBackingStorePixelRatio,
        context.msBackingStorePixelRatio,
        context.oBackingStorePixelRatio,
        context.backingStorePixelRatio,
        1
    ].reduce(function(a, b) { return a || b });

    return devicePixelRatio / backingStoreRatio;
}