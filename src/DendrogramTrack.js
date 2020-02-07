import range from 'lodash/range';
import d3 from './utils/d3.js';
import { matrixToTree } from './utils/tree.js';

/**
 * Dendrogram for rendering an attribute of hierarchical structure.
 * @prop {object} two An instance of the Two class.
 * @prop {number} left The left position of this view.
 * @prop {number} top The top position of this view.
 * @prop {number} width The width of this view.
 * @prop {number} height The height of this view.
 * @prop {array} rowInfo Array of JSON objects, one object for each row.
 * @prop {object} attribute The name and type of data field.
 * @prop {boolean} isLeft Is this view on the left side of the track?
 */
export function dendrogramTrack(props) {
    const {
        two,
        left, top, width, height,
        rowInfo, 
        attribute,
        isLeft
    } = props;

    // Data, layouts and styles
    const { name: field, type: fieldType } = attribute;
    const hierarchyData = matrixToTree(rowInfo.map(d => d[field]));
    const root = d3.hierarchy(hierarchyData);

    const treeLayout = d3.cluster()
        .size([height, width])
        .separation(() => 1);
    treeLayout(root);

    const descendants = root.descendants();

    let pathFunction;
    pathFunction = (d) => {
        return two.makePath(
            width - (top + d.parent.y), left + d.parent.x, // M
            width - (top + d.y), left + d.parent.x, // H d.x
            width - (top + d.y), left + d.x, // M
            width - (top + d.y), left + d.parent.x // V d.parent.y
        );
    }
    
    let nodeTransformFunction;
    nodeTransformFunction = (d) => { 
        return [d.x + left, d.y + top]; 
    }

    const nodePoints = [];

    descendants.forEach((d, i) => {
        if(i > 0) {
            const path = pathFunction(d);
            path.stroke = "#555";
            path.opacity = 0.6;
            path.linewidth = 1.5;
            path.rotation = Math.PI/2;
        }
        nodePoints.push(nodeTransformFunction(d));
    });

    
}