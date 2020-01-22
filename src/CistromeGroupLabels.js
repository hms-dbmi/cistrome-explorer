import React, { useEffect, useState } from 'react';
import PubSub from 'pubsub-js';

import { scaleBand } from 'd3-scale';
import { interpolateViridis } from "d3-scale-chromatic";

import { GLOBAL_X_RANGE, TRACK_ROW_INFO, TRACK_POSITION, TRACK_DIMENSIONS } from './constants.js';

import './CistromeGroupLabels.css';

export default function CistromeGroupLabels() {

    const [x1, setX1] = useState(0);
    const [x0, setX0] = useState(0);
    const [y, setY] = useState(0);
    const [height, setHeight] = useState(30);

    const [rowNames, setRowNames] = useState([]);
    const columnIndexA = 5;
    const columnIndexB = 6;

    useEffect(() => {
        const tokenGlobalXRange = PubSub.subscribe(GLOBAL_X_RANGE, (msg, data) => {
            setX1(data[1]);
        });

        const tokenTrackPosition = PubSub.subscribe(TRACK_POSITION, (msg, data) => {
            setX0(data[0]);
            setY(data[1]);
        });

        const tokenTrackDimensions = PubSub.subscribe(TRACK_DIMENSIONS, (msg, data) => {
            setHeight(data[1]);
        });

        const tokenRowInfo = PubSub.subscribe(TRACK_ROW_INFO, (msg, data) => {
            setRowNames(data.map(d => d.split("\t")));
        });

        return (() => {
            PubSub.unsubscribe(tokenGlobalXRange);
            PubSub.unsubscribe(tokenTrackPosition);
            PubSub.unsubscribe(tokenTrackDimensions);
            PubSub.unsubscribe(tokenRowInfo);
        });
    });

    const categoryScaleA = scaleBand().domain(Array.from(new Set(rowNames.map(d => d[columnIndexA])))).range([0, 1]);
    const categoryScaleB = scaleBand().domain(Array.from(new Set(rowNames.map(d => d[columnIndexB])))).range([0, 1]);
    
    const colorScale = interpolateViridis;

    return (
        <div>
            <div style={{
                position: 'absolute',
                top: `${y}px`,
                left: `${x0 + x1 + 15}px`, 
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
                left: `${x0 + x1 + 50}px`, 
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