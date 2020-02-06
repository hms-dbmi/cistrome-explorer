import d3 from './d3.js';

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


/**
 * Clean up an HTML canvas element by removing its event listeners.
 * @param {object} canvas A `<canvas/>` element.
 */
export function teardownCanvas(canvas) {
    const canvasSelection = d3.select(canvas);
    canvasSelection.on("mousemove", null);
    canvasSelection.on("mouseout", null);
    canvasSelection.on("click", null);
}