import React, { useRef, useCallback, useEffect } from "react";
import d3 from "./utils/d3.js";
import Two from "./utils/two.js";
import { HIGHLIGHTING_COLOR } from "./utils/linking.js";

/* TODO: Add these to HiGlassMeta options. */
const BAND_TYPE = [
    'band',
    'line',
    'curved-band',
    'curved-line'
][1];
const POINT_RADIUS = 0;
const TRACK_PADDING = 0;
const BAND_PADDING = 10;

/**
 * Component for visualization of row info quantitative or nominal attribute values.
 * @prop {number} left The left position of this view.
 * @prop {number} top The top position of this view.
 * @prop {number} width The width of this view.
 * @prop {number} height The height of this view.
 * @prop {array} leftSelectedRows The array of selected indices on the left track.
 * @prop {array} rightSelectedRows The array of selected indices on the right track.
 * @prop {array} selectedRows The array of selected indices. 
 * @prop {array} highlitRows The array of highlit indices.
 * @prop {function} drawRegister The function for child components to call to register their draw functions.
 */
export default function TrackRowInfoVisBand(props) {
    const {
        left, top, width, height,
        leftSelectedRows,
        rightSelectedRows,
        selectedRows,
        highlitRows,
        drawRegister,
    } = props;

    const divRef = useRef();
    const canvasRef = useRef();

    const yScaleLeft = d3.scaleBand()
        .domain(leftSelectedRows ?? [])
        .range([0, height]);

    const yScaleRight = d3.scaleBand()
        .domain(rightSelectedRows ?? [])
        .range([0, height]);

    const bandWidth = leftSelectedRows ? height / leftSelectedRows.length : 1;

    const draw = useCallback((domElement) => {
        const two = new Two({
            width,
            height,
            domElement
        });
        
        const bandColor = (i) => highlitRows?.indexOf(i) !== -1 ? HIGHLIGHTING_COLOR : "lightgray";
        const lineWidth = (i) => highlitRows?.indexOf(i) !== -1 ? 2 : 1;
        const renderBand = {
            'band': (i) => {
                const band = two.makePath(
                    TRACK_PADDING, yScaleLeft(i),
                    TRACK_PADDING + BAND_PADDING, yScaleLeft(i),
                    width - BAND_PADDING - TRACK_PADDING, yScaleRight(i),
                    width - TRACK_PADDING, yScaleRight(i),
                    width - TRACK_PADDING, yScaleRight(i) + bandWidth,
                    width - BAND_PADDING - TRACK_PADDING, yScaleRight(i) + bandWidth,
                    TRACK_PADDING + BAND_PADDING, yScaleLeft(i) + bandWidth,
                    TRACK_PADDING, yScaleLeft(i) + bandWidth,
                )
                band.fill = bandColor(i);
                band.stroke = "transparent";
                band.linewidth = 1;
                band.opacity = 0.4;
            },
            'line': (i) => {
                if(BAND_PADDING !== 0) {
                    const lineStart = two.makeLine(
                        TRACK_PADDING, 
                        yScaleLeft(i) + bandWidth / 2.0,
                        TRACK_PADDING + BAND_PADDING, 
                        yScaleLeft(i) + bandWidth / 2.0,
                    );
                    lineStart.stroke = bandColor(i);
                    lineStart.linewidth = lineWidth(i);
                    lineStart.opacity = 0.5;

                    const lineEnd = two.makeLine(
                        width - BAND_PADDING - TRACK_PADDING, 
                        yScaleRight(i) + bandWidth / 2.0,
                        width - TRACK_PADDING,
                        yScaleRight(i) + bandWidth / 2.0,
                    );
                    lineEnd.stroke = bandColor(i);
                    lineEnd.linewidth = lineWidth(i);
                    lineEnd.opacity = 0.5;
                }

                const lineMid = two.makeLine(
                    TRACK_PADDING + BAND_PADDING,
                    yScaleLeft(i) + bandWidth / 2.0,
                    width - TRACK_PADDING - BAND_PADDING,
                    yScaleRight(i) + bandWidth / 2.0,
                );
                lineMid.stroke = bandColor(i);
                lineMid.linewidth = lineWidth(i);
                lineMid.opacity = 0.5;

                if(POINT_RADIUS) {
                    const pointLeft = two.makeCircle(0, yScaleLeft(i) + bandWidth / 2.0, POINT_RADIUS);
                    pointLeft.fill = "black";
                    pointLeft.linewidth = 1;
                    pointLeft.opacity = 0.5;

                    const pointRight = two.makeCircle(width, yScaleLeft(i) + bandWidth / 2.0, POINT_RADIUS);
                    pointRight.fill = "black";
                    pointRight.linewidth = 1;
                    pointRight.opacity = 0.5;
                }
            }
        }

        // line connections
        if(leftSelectedRows && rightSelectedRows) {
            leftSelectedRows.forEach(index => {
                renderBand[BAND_TYPE](index);
            });
        }

        // left and right axes
        const axisL = two.makeLine(TRACK_PADDING, 0, TRACK_PADDING, height);
        axisL.stroke = "black";
        axisL.linewidth = 2;
        const axisR = two.makeLine(width - TRACK_PADDING, 0, width - TRACK_PADDING, height);
        axisR.stroke = "black";
        axisR.linewidth = 2;

        two.update();
        return two.teardown;
    });
    
    useEffect(() => {
        const canvas = canvasRef.current;
        const div = divRef.current;
        const teardown = draw(canvas);

        // Clean up.
        return () => {
            teardown();
        };
    }, [top, left, width, height, leftSelectedRows]);

    drawRegister("TrackRowInfoVisBand", draw);

    return (
        <div
            ref={divRef}
            style={{
                position: 'relative',
                width: `${width}px`,
                height: `${height}px`,
            }}
        >
            <canvas
                ref={canvasRef}
                style={{
                    top: 0,
                    left: 0, 
                    width: `${width}px`,
                    height: `${height}px`,
                    position: 'relative'
                }}
            />
        </div>
    );
}