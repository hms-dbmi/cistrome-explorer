import React, { useRef, createRef, useState, useEffect} from 'react';
import d3 from './utils/d3.js';
import { modifyItemInArray } from './utils/array.js';
import TrackRowInfoVisBar from './TrackRowInfoVisBar.js';
import TrackRowInfoVisLink from './TrackRowInfoVisLink.js';
import TrackRowInfoVisDendrogram from './TrackRowInfoVisDendrogram.js';
import './TrackResizer.scss';

const fieldTypeToVisComponent = {
    "nominal": TrackRowInfoVisBar,
    "quantitative": TrackRowInfoVisBar,
    "url": TrackRowInfoVisLink,
    "tree": TrackRowInfoVisDendrogram
};

/**
 * Parent component for visualization of row info attribute values.
 * @prop {string} viewId The viewId for the horizontal-multivec track.
 * @prop {string} trackId The trackId for the horizontal-multivec track.
 * @prop {number} trackX The track horizontal offset.
 * @prop {number} trackY The track vertical offset.
 * @prop {number} trackWidth The track width.
 * @prop {number} trackHeight The track height.
 * @prop {array} rowInfo Array of JSON objects, one object for each row.
 * @prop {array} rowInfoAttributes Array of JSON object, one object for the names and types of each attribute.
 * @prop {string} rowInfoPosition The value of the `rowInfoPosition` option.
 * @prop {function} onSortRows The function to call upon a sort interaction.
 * @prop {function} onSearchRows The function to call upon a search interaction.
 * @prop {function} drawRegister The function for child components to call to register their draw functions.
 */
export default function TrackRowInfo(props) {

    const {
        viewId, trackId,
        trackX, trackY,
        trackWidth, trackHeight, 
        rowInfo, 
        rowInfoAttributes,
        rowInfoPosition,
        onSortRows,
        onSearchRows,
        drawRegister
    } = props;

    // Dimensions
    const defaultUnitWidth = 100;
    const isLeft = rowInfoPosition === "left";
    const top = trackY;
    const height = trackHeight;
    let unitWidths = [];    // Store default unit width for each vertical tracks
    rowInfoAttributes.forEach((fieldInfo) => {
        unitWidths.push({
            field: fieldInfo.field,
            type: fieldInfo.type,
            width: defaultUnitWidth
        });
    });
    const [widths, setWidths] = useState(unitWidths);
    const width = d3.sum(widths.map(d => d.width));
    const left = isLeft ? trackX - width : trackX + trackWidth;

    const divRef = useRef();
    const resizerRef = useRef([...Array(rowInfoAttributes.length)].map(() => createRef()));
    const [resizingIndex, setResizingIndex] = useState(-1);
    
    // Determine position of each dimension.
    let trackProps = [];
    let currentLeft = 0;
    rowInfoAttributes.forEach((attribute, i) => {
        const fieldInfo = isLeft ? rowInfoAttributes[rowInfoAttributes.length - i - 1] : attribute;
        const width = widths.find(d => d.field === fieldInfo.field && d.type === fieldInfo.type).width;

        trackProps.push({
            left: currentLeft, top: 0, width, height,
            fieldInfo,
            isLeft
        });
        currentLeft += width;
    });

    console.log(trackProps);
    
    let resizers = trackProps.map((d, i) => {
        const resizerWidth = 4, resizerHeight = 10, margin = 2;
        return (
            <div
                ref={resizerRef.current[i]}
                key={i}
                className="visualization-resizer"
                style={{
                    top: `${d.top + (d.height + resizerHeight) / 2.0}px`,
                    left: `${isLeft ? d.left + margin : d.left + d.width - resizerWidth - margin}px`,
                    height: `${resizerHeight}px`,
                    width: `${resizerWidth}px`,
                    // visibility: mouseX !== null ? "visible" : "hidden"
                }}
            />
        );
    });

    useEffect(() => {
        const div = divRef.current;
        d3.select(div).on("mouseup", () => {
            console.log("MOUSE_UP"); 
            setResizingIndex(-1);
        });

        d3.select(div).on("mousemove", () => {
            if(resizingIndex !== -1) {
                const { left } = trackProps[resizingIndex];
                const [mouseX, mouseY] = d3.mouse(div);
                console.log("MOUSE_MOVE:", resizingIndex);
                let newWidth = isLeft ? trackProps[resizingIndex].width - (mouseX - left) : (mouseX - left);
                const minWidth = 50;
                if(newWidth < minWidth) {
                    newWidth = minWidth;
                }
                setWidths(modifyItemInArray(widths, resizingIndex, {
                    field: rowInfoAttributes[resizingIndex].field,
                    type: rowInfoAttributes[resizingIndex].type,
                    width: newWidth
                }));
            }
        });
        
        resizerRef.current.forEach((resizer, i) => {
            d3.select(resizer.current).on("mousedown", () => {
                const selectedIndex = isLeft ? resizerRef.current.length - i - 1 : i;
                console.log("MOUSE_DOWN:", selectedIndex);
                setResizingIndex(selectedIndex);
            });
        });
    })

    console.log("TrackRowInfo.render");
    return (
        <div
            ref={divRef}
            className="cistrome-hgw-child"
            style={{
                top: `${top}px`,
                left: `${left}px`, 
                width: `${width + 20}px`,
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
                    viewId: viewId,
                    trackId: trackId,
                    rowInfo: rowInfo,
                    onSortRows: onSortRows,
                    onSearchRows: onSearchRows,
                    drawRegister: drawRegister,
                }
            ))}
            {resizers}
        </div>
    );
}