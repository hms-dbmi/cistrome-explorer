import React from 'react';
import { scaleBand } from 'd3-scale';
import { interpolateViridis } from "d3-scale-chromatic";

import './TrackRowInfo.scss';

export default function TrackRowInfo(props) {

    const {
        x0, x1, y1, height, rowInfo, 
        infoAttrPrimary, infoAttrSecondary,
        rowInfoPosition
    } = props;

    const width = 15;
    const xMargin = 80;
    const xMarginInitial = 5;

    let xPrimary, xSecondary;
    if(rowInfoPosition === "left") {
        xPrimary = x1 - xMarginInitial - width - 5;
        xSecondary = x1 - xMarginInitial - width - xMargin - width - 5;
    } else if(rowInfoPosition === "right") {
        xPrimary = x0 + x1 + xMarginInitial;
        xSecondary = x0 + x1 + xMarginInitial + width + xMargin;
    }

    const scales = {
        primary: scaleBand()
            .domain(Array.from(new Set(rowInfo.map(d => d[infoAttrPrimary]))))
            .range([0, 1]),
        secondary: scaleBand()
            .domain(Array.from(new Set(rowInfo.map(d => d[infoAttrSecondary]))))
            .range([0, 1])
    };
    
    const colorScale = interpolateViridis;

    return (
        <div>
            <div style={{
                position: 'absolute',
                top: `${y1}px`,
                left: `${xPrimary}px`, 
                width: '20px',
                height: `${height}px`,
            }}>
                {rowInfo.map((d, i) => (
                    <span 
                        key={i}
                        className="row-name"
                        style={{
                            height:`${height/rowInfo.length}px`,
                            backgroundColor: colorScale(scales.primary(d[infoAttrPrimary]))
                        }}
                        title={d[infoAttrPrimary]}
                    />
                ))}
            </div>
            <div style={{
                position: 'absolute',
                top: `${y1}px`,
                left: `${xSecondary}px`, 
                width: '20px',
                height: `${height}px`,
            }}>
                {rowInfo.map((d, i) => (
                    <span 
                        key={i}
                        className="row-name"
                        style={{
                            height:`${height/rowInfo.length}px`,
                            backgroundColor: colorScale(scales.secondary(d[infoAttrSecondary]))
                        }}
                        title={d[infoAttrSecondary]}
                    />
                ))}
            </div>
        </div>
    );
}