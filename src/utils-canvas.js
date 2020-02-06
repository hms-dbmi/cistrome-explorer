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
 * Initialize an HTML canvas element, scaling it for retina screens.
 * @param {object} canvasRef A React ref corresponding to a `<canvas/>` element.
 * @returns {object} An object containing the canvas, d3.select(canvas), the context, width, and height.
 */
export function setupCanvas(canvasRef) {
    const canvas = canvasRef.current;
    const canvasSelection = d3.select(canvas);
    const context = canvas.getContext('2d');
    const ratio = getRetinaRatio(context);
    const scaledWidth = canvas.clientWidth * ratio;
    const scaledHeight = canvas.clientHeight * ratio;
    canvas.setAttribute("width", scaledWidth);
    canvas.setAttribute("height", scaledHeight);
    context.scale(ratio, ratio);

    return {
        canvas,
        canvasSelection,
        context,
        canvasWidth: canvas.clientWidth,
        canvasHeight: canvas.clientHeight
    };
}

/**
 * Clean up an HTML canvas element by removing its event listeners.
 * @param {object} canvasRef A React ref corresponding to a `<canvas/>` element.
 */
export function teardownCanvas(canvasRef) {
    const canvas = canvasRef.current;
    const canvasSelection = d3.select(canvas);
    canvasSelection.on("mousemove", null);
    canvasSelection.on("mouseout", null);
    canvasSelection.on("click", null);
}