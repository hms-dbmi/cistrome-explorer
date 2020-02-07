import { verticalBarTrack } from "./VerticalBarTrack";
import { dendrogramTrack } from "./DendrogramTrack";

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
    const { type } = props.attribute;

    switch(type) {
        case "nominal": 
        case "quantitative":
            verticalBarTrack(props);
            break;
        case "tree":
            dendrogramTrack(props);
            break;
    }
}