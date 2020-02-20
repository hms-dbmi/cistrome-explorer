

/**
 * Common function for rendering the title text for a vertical visualization component.
 * @param {string} text The title text content.
 * @param {object} options Options for drawing title.
 * @param {Two} options.two The object of Two.js class.
 * @param {boolean} options.isLeft In this view placed on the left side of the track?
 * @param {boolean} options.isNominal Is this view visualizes a nominal value?
 * @param {number} options.width Width of this view.
 */
export function drawVisTitle(text, options) {
    const { two, isLeft, isNominal, width, titleSuffix } = options;
    
    const margin = 5;
    const barAreaWidth = isNominal ? 20 : width - 20;
    const titleFontSize = 12;
    const titleLeft = (isLeft ? margin : width - margin);
    const titleRotate = isLeft ? -Math.PI/2 : Math.PI/2;
    const fullText = text + titleSuffix;

    const title = two.makeText(titleLeft, 0, 12, barAreaWidth, fullText);
    title.fill = "#9A9A9A";
    title.fontsize = titleFontSize;
    title.align = isLeft ? "end" : "start";
    title.baseline = "top";
    title.rotation = titleRotate;
}
