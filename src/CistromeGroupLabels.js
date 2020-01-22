import React, { useEffect, useState } from 'react';
import PubSub from 'pubsub-js';

import { GLOBAL_X_RANGE, GLOBAL_Y_RANGE, TRACK_ROW_INFO, TRACK_POSITION, TRACK_DIMENSIONS } from './constants.js';

function mySubscriber(msg, data) {
    console.log(data);
}
export default function CistromeGroupLabels() {

    const [x, setX] = useState(0);
    const [x0, setX0] = useState(0);
    const [y, setY] = useState(0);

    const [height, setHeight] = useState(30);

    useEffect(() => {
        const tokenGlobalXRange = PubSub.subscribe(GLOBAL_X_RANGE, (msg, data) => {
            setX(data[1]);
        });

        const tokenTrackPosition = PubSub.subscribe(TRACK_POSITION, (msg, data) => {
            setX0(data[0]);
            setY(data[1]);
        });

        const tokenTrackDimensions = PubSub.subscribe(TRACK_DIMENSIONS, (msg, data) => {
            setHeight(data[1]);
        });

        return (() => {
            PubSub.unsubscribe(tokenGlobalXRange);
            PubSub.unsubscribe(tokenTrackPosition);
            PubSub.unsubscribe(tokenTrackDimensions);
        });
    });

    console.log("CistromeGroupLabels.render");

    return (
        <div style={{
            position: 'absolute',
            top: `${y}px`,
            left: `${x + x0}px`, 
            backgroundColor: 'red', 
            width: '30px',
            height: `${height}px`,
        }}></div>
    );
}