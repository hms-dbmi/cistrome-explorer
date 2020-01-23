import React from 'react';
import { scaleBand } from 'd3-scale';
import { interpolateViridis } from "d3-scale-chromatic";

import './CistromeGroupLabels.scss';

export default function CistromeGroupLabels(props) {

    const {
        x, y, height, rowNames
    } = props;

    const columnIndexA = 5;
    const columnIndexB = 6;

    const categoryScaleA = scaleBand()
        .domain(Array.from(new Set(rowNames.map(d => d[columnIndexA]))))
        .range([0, 1]);
    const categoryScaleB = scaleBand()
        .domain(Array.from(new Set(rowNames.map(d => d[columnIndexB]))))
        .range([0, 1]);
    
    const colorScale = interpolateViridis;

    return (
        <div>
            <div style={{
                position: 'absolute',
                top: `${y}px`,
                left: `${x + 15}px`, 
                width: '20px',
                height: `${height}px`,
            }}>
                {rowNames.map((name, i) => (
                    <span 
                        key={i}
                        className="row-name"
                        style={{
                            height:`${height/rowNames.length}px`,
                            backgroundColor: colorScale(categoryScaleA(name[columnIndexA]))
                        }}
                        title={name[columnIndexA]}
                    />
                ))}
            </div>
            <div style={{
                position: 'absolute',
                top: `${y}px`,
                left: `${x + 50}px`, 
                width: '20px',
                height: `${height}px`,
            }}>
                {rowNames.map((name, i) => (
                    <span 
                        key={i}
                        className="row-name"
                        style={{
                            height:`${height/rowNames.length}px`,
                            backgroundColor: colorScale(categoryScaleB(name[columnIndexB]))
                        }}
                        title={name[columnIndexB]}
                    />
                ))}
            </div>
        </div>
    );
}