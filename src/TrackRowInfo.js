import React from 'react';

import TrackRowInfoVisBar from './TrackRowInfoVisBar.js';
import TrackRowInfoVisLink from './TrackRowInfoVisLink.js';
import TrackRowInfoVisDendrogram from './TrackRowInfoVisDendrogram.js';

const fieldTypeToVisComponent = {
    "nominal": TrackRowInfoVisBar,
    "quantitative": TrackRowInfoVisBar,
    "url": TrackRowInfoVisLink,
    "tree": TrackRowInfoVisDendrogram
};

/**
 * Parent component for visualization of row info attribute values.
 * @prop {number} trackX The track horizontal offset.
 * @prop {number} trackY The track vertical offset.
 * @prop {number} trackWidth The track width.
 * @prop {number} trackHeight The track height.
 * @prop {array} rowInfo Array of JSON objects, one object for each row.
 * @prop {array} rowInfoAttributes Array of JSON object, one object for the names and types of each attribute.
 * @prop {string} rowInfoPosition The value of the `rowInfoPosition` option.
 * @prop {function} drawRegister The function for child components to call to register their draw functions.
 */
export default function TrackRowInfo(props) {

    const {
        trackX, trackY,
        trackWidth, trackHeight, 
        rowInfo, 
        rowInfoAttributes,
        rowInfoPosition,
        drawRegister
    } = props;

    // Dimensions
    const isLeft = rowInfoPosition === "left";
    const top = trackY;
    const unitWidth = 100;
    const width = unitWidth * rowInfoAttributes.length;
    const height = trackHeight;
    const left = isLeft ? trackX - width : trackX + trackWidth;

    // Determine position of each dimension.
    let trackProps = [], xDomain = [], xRange = [];
    rowInfoAttributes.forEach((attribute, i) => {
        const fieldInfo = isLeft ? rowInfoAttributes[rowInfoAttributes.length - i - 1] : attribute;
        let currentLeft = unitWidth * i;

        trackProps.push({
            left: currentLeft, top: 0, width: unitWidth, height: height,
            fieldInfo,
            isLeft
        });

        // Domain and range for mouse event.
        xDomain.push(currentLeft + unitWidth);
        xRange.push(fieldInfo.field);
    });

    console.log("TrackRowInfo.render");
    return (
        <div 
            className="cistrome-hgw-child"
            style={{
                top: `${top}px`,
                left: `${left}px`, 
                width: `${width}px`,
                height: `${height}px`,
            }}
        >
            {trackProps.map((d, i) => React.createElement(
                fieldTypeToVisComponent[d.fieldInfo.type],
                {
                    key: i,
                    left: d.left,
                    top: d.top,
                    width: d.width,
                    height: d.height,
                    isLeft: d.isLeft,
                    fieldInfo: d.fieldInfo,
                    rowInfo: rowInfo,
                    drawRegister: drawRegister,
                }
            ))}
        </div>
    );
}