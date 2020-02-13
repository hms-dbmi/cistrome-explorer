

/**
 * Common function for rendering the title text for a vertical visualization component.
 * @prop {string} text The title text content.
 * @prop {object} options
 * @prop {Two} options.two
 * @prop {boolean} options.isLeft
 * @prop {boolean} options.isNominal
 * @prop {number} options.width
 */
export function drawVisTitle(text, { two, isLeft, isNominal, width }) {
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
