import { select as d3_select } from 'd3-selection';

/**
 * Get the retina ratio to be able to scale up a canvas context.
 * @private
 * @param {context} c The canvas context.
 * @returns {float} The ratio.
 */
export function getRetinaRatio(c) {
    let devicePixelRatio = window.devicePixelRatio || 1;
    let backingStoreRatio = [
        c.webkitBackingStorePixelRatio,
        c.mozBackingStorePixelRatio,
        c.msBackingStorePixelRatio,
        c.oBackingStorePixelRatio,
        c.backingStorePixelRatio,
        1
    ].reduce(function(a, b) { return a || b });

    return devicePixelRatio / backingStoreRatio;
}

export function setupCanvas(canvasRef) {
    const canvas = canvasRef.current;
    const context = canvas.getContext('2d');
    const ratio = getRetinaRatio(context);
    const scaledWidth = canvas.clientWidth * ratio;
    const scaledHeight = canvas.clientHeight * ratio;
    canvas.setAttribute("width", scaledWidth);
    canvas.setAttribute("height", scaledHeight);
    context.scale(ratio, ratio);

    const canvasSelection = d3_select(canvas);

    return {
        canvas,
        canvasSelection,
        context,
        canvasWidth: canvas.clientWidth,
        canvasHeight: canvas.clientHeight
    };
}

export function teardownCanvas(canvasRef) {
    const canvas = canvasRef.current;
    const canvasSelection = d3_select(canvas);
    canvasSelection.on("mousemove", null);
    canvasSelection.on("mouseout", null);
    canvasSelection.on("click", null);
}