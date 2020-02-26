import { TwoText, TwoRectangle } from "./two.js";

/**
 * Common function for rendering the title text for a vertical visualization component.
 * @param {string} text The title text content.
 * @param {object} options Options for drawing title.
 * @param {Two} options.two The object of Two.js class.
 * @param {boolean} options.isLeft In this view placed on the left side of the track?
 * @param {number} options.width Width of this view.
 * @param {string} options.titleSuffix The suffix of a title, information about sorting and filtering status. Optional.
 * @param {string} options.backgroundColor If defined, creates a rect behind the text element. Optional. By default, "#FFF".
 */
export function drawVisTitle(text, options) {
    const {
        two,
        isLeft,
        width,
        height,
        titleSuffix = "",
        backgroundColor = "white"
    } = options;
    
    const margin = 4;
    const titleFontSize = 12;
    const titleLeftInitial = isLeft ? margin : (width - margin);
    const titleRotate = isLeft ? -Math.PI/2 : Math.PI/2;
    const fullText = `${text}${titleSuffix}`;

    const title = new TwoText(titleLeftInitial, 0, width, height, fullText);
    title.fill = "#9A9A9A";
    title.fontsize = titleFontSize;
    title.align = isLeft ? "end" : "start";
    title.baseline = isLeft ? "top" : "bottom";
    title.rotation = titleRotate;

    const titleDims = two.measureText(title);
    const titleLeft = isLeft ? titleLeftInitial : titleLeftInitial - titleDims.height;
    title.x = titleLeft;

    if(backgroundColor) {
        const rect = new TwoRectangle(titleLeft - 2, 0, titleDims.height + 6, titleDims.width);
        rect.fill = backgroundColor;
        rect.stroke = null;
        rect.opacity = 0.8;

        two.append(rect);
    }

    two.append(title);
}
