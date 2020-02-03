import { select as d3_select } from 'd3-selection';

/**
 * Get the retina ratio to be able to scale up a canvas context.
 * @param {Object} context The canvas context.
 * @returns {float} The ratio.
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
 * @param {ref} canvasRef A React ref corresponding to a `<canvas/>` element.
 * @returns {Object} An object containing the canvas, d3_select(canvas), the context, width, and height.
 */
export function setupCanvas(canvasRef) {
    const canvas = canvasRef.current;
    const canvasSelection = d3_select(canvas);
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
 * @param {ref} canvasRef A React ref corresponding to a `<canvas/>` element.
 */
export function teardownCanvas(canvasRef) {
    const canvas = canvasRef.current;
    const canvasSelection = d3_select(canvas);
    canvasSelection.on("mousemove", null);
    canvasSelection.on("mouseout", null);
    canvasSelection.on("click", null);
}