

/**
 * Common function for rendering the title text for a vertical visualization component.
 * @param {string} text The title text content.
 * @param {object} options
 * @param {Two} options.two
 * @param {boolean} options.isLeft
 * @param {boolean} options.isNominal
 * @param {number} options.width
 */
export function drawVisTitle(text, options) {
    const { two, isLeft, isNominal, width } = options;
    
    const margin = 5;
    const barAreaWidth = isNominal ? 20 : width - 20;
    const titleFontSize = 12;
    const titleLeft = (isLeft ? margin : width - margin);
    const titleRotate = isLeft ? -Math.PI/2 : Math.PI/2;

    const title = two.makeText(titleLeft, 0, 12, barAreaWidth, text);
    title.fill = "#9A9A9A";
    title.fontsize = titleFontSize;
    title.align = isLeft ? "end" : "start";
    title.baseline = "top";
    title.rotation = titleRotate;
}
