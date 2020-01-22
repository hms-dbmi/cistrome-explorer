import React, { useEffect, useState } from 'react';
import PubSub from 'pubsub-js';

import { scaleBand } from 'd3-scale';
import { interpolateViridis } from "d3-scale-chromatic";

import { GLOBAL_X_RANGE, GLOBAL_Y_RANGE, TRACK_ROW_INFO, TRACK_POSITION, TRACK_DIMENSIONS } from './constants.js';

import './CistromeGroupLabels.css';

function consoleSubscriber(msg, data) {
    console.log(data);
}

export default function CistromeGroupLabels() {

    const [x1, setX1] = useState(0);
    const [x0, setX0] = useState(0);
    const [y, setY] = useState(0);
    const [height, setHeight] = useState(30);

    const [rowNames, setRowNames] = useState(null);
    const columnIndex = 6;

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
            setRowNames(data.map(d => d.split("\t")).map(d => d[columnIndex]));
        });

        return (() => {
            PubSub.unsubscribe(tokenGlobalXRange);
            PubSub.unsubscribe(tokenTrackPosition);
            PubSub.unsubscribe(tokenTrackDimensions);
            PubSub.unsubscribe(tokenRowInfo);
        });
    });

    const categoryScale = scaleBand().domain(Array.from(new Set(rowNames))).range([0, 1]);
    const colorScale = interpolateViridis;


    return (
        <div style={{
            position: 'absolute',
            top: `${y}px`,
            left: `${x0 + x1}px`, 
            backgroundColor: 'silver', 
            width: '90px',
            height: `${height}px`,
        }}>
            {rowNames ? rowNames.map((name, i) => (
                <span 
                    key={i}
                    className="row-name"
                    style={{
                        height:`${height/rowNames.length}px`,
                        backgroundColor: colorScale(categoryScale(name))
                    }}
                    title={name}
                >
                </span>
            )): (null)}
        </div>
    );
}