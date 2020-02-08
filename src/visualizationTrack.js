import { verticalBarTrack } from "./verticalBarTrack";
import { dendrogramTrack } from "./dendrogramTrack";

/**
 * Visualization for rendering an attribute based on its type.
 * @prop {object} two An instance of the Two class.
 * @prop {number} left The left position of this view.
 * @prop {number} top The top position of this view.
 * @prop {number} width The width of this view.
 * @prop {number} height The height of this view.
 * @prop {array} rowInfo Array of JSON objects, one object for each row.
 * @prop {object} attribute The name and type of data field.
 * @prop {boolean} isLeft Is this view on the left side of the track?
 */
export function visualizationTrack(props){
    const {
        two,
        left, top, width,
        attribute,
        isLeft
    } = props;

    // Data, layouts and styles
    const { name: field, type: fieldType } = attribute;
    const isNominal = fieldType === "nominal";
    const barAreaWidth = isNominal ? 20 : width - 20;
    const margin = 5;
    const titleFontSize = 12;

    // Render proper visualization
    switch(fieldType) {
        case "nominal": 
        case "quantitative":
            verticalBarTrack(props);
            break;
        case "tree":
            dendrogramTrack(props);
            break;
    }

    // Title
    const titleLeft = left + (isLeft ? margin : width - margin);
    const titleRotate = isLeft ? -Math.PI/2 : Math.PI/2;

    // Draw a title of each dimension
    const titleText = field;
    const title = two.makeText(titleLeft, top, 12, barAreaWidth, titleText);
    title.fill = "#9A9A9A";
    title.fontsize = titleFontSize;
    title.align = isLeft ? "end" : "start";
    title.baseline = "top";
    title.rotation = titleRotate;
}