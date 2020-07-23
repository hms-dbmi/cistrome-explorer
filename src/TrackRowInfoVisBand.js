import React, { useRef, useCallback, useEffect } from "react";
import d3 from "./utils/d3.js";
import Two from "./utils/two.js";

// TODO: Add these to HiGlassMeta options.
const BAND_TYPE = [
    'band',
    'line',
    'curved-band',
    'curved-line'
][0];
const POINT_RADIUS = 0;
const TRACK_PADDING = 0;
const BAND_PADDING = 10;

/**
 * 
 * @param {array} leftSelectedRows
 */
export default function TrackRowInfoVisBand(props) {
    const {
        left, top, width, height,
        leftSelectedRows,
        rightSelectedRows,
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
                band.fill = "lightgray";
                band.stroke = "transparent";
                band.linewidth = 1;
                band.opacity = 0.3;
            },
            'line': (i) => {
                if(BAND_PADDING !== 0) {
                    const lineStart = two.makeLine(
                        TRACK_PADDING, 
                        yScaleLeft(i) + bandWidth / 2.0,
                        TRACK_PADDING + BAND_PADDING, 
                        yScaleLeft(i) + bandWidth / 2.0,
                    );
                    lineStart.stroke = "lightgray";
                    lineStart.linewidth = 1;
                    lineStart.opacity = 0.5;

                    const lineEnd = two.makeLine(
                        width - BAND_PADDING - TRACK_PADDING, 
                        yScaleRight(i) + bandWidth / 2.0,
                        width - TRACK_PADDING,
                        yScaleRight(i) + bandWidth / 2.0,
                    );
                    lineEnd.stroke = "lightgray";
                    lineEnd.linewidth = 1;
                    lineEnd.opacity = 0.5;
                }

                const lineMid = two.makeLine(
                    TRACK_PADDING + BAND_PADDING,
                    yScaleLeft(i) + bandWidth / 2.0,
                    width - TRACK_PADDING - BAND_PADDING,
                    yScaleRight(i) + bandWidth / 2.0,
                );
                lineMid.stroke = "lightgray";
                lineMid.linewidth = 1;
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
    
    drawRegister("TrackRowInfoVisBand", draw);

    useEffect(() => {
        const canvas = canvasRef.current;
        const div = divRef.current;
        const teardown = draw(canvas);

        d3.select(canvas).on("mousemove", () => {
            const [] = d3.mouse(canvas);
            // TODO:
        });

        // Handle mouse leave.
        // d3.select(canvas).on("mouseout", destroyTooltip);
        // d3.select(div).on("mouseleave", () => setHoverValue(null));

        // Clean up.
        return () => {
            teardown();
            d3.select(div).on("mouseleave", null);
        };
    }, [leftSelectedRows, rightSelectedRows, top, left, width, height]);

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