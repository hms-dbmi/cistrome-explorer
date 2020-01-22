import React, { useRef, useEffect } from 'react';
import PubSub from 'pubsub-js';

import { HiGlassComponent } from 'higlass';
import register from 'higlass-register';
import StackedBarTrack from 'higlass-multivec/es/StackedBarTrack.js';

import CistromeGroupLabels from './CistromeGroupLabels.js';

import { GLOBAL_X_RANGE, GLOBAL_Y_RANGE, TRACK_ROW_INFO, TRACK_POSITION, TRACK_DIMENSIONS } from './constants.js';

import 'higlass/dist/hglib.css';
import './CistromeHGW.css';

register({
    name: 'StackedBarTrack',
    track: StackedBarTrack,
    config: StackedBarTrack.config,
});

const hgOptions = {
    bounded: true,
    pixelPreciseMarginPadding: true,
    containerPaddingX: 0,
    containerPaddingY: 0,
    sizeMode: 'default'
};

function getHorizontalMultivecTrackId(viewObj) {
    let trackId = null;
    if(viewObj && viewObj.tracks && viewObj.tracks.center && viewObj.tracks.center.length > 0) {
        for(let trackObj of viewObj.tracks.center) {
            if(trackObj.type === "horizontal-multivec" && trackObj.uid) {
                trackId = trackObj.uid;
                break;
            }
        }
    }
    return trackId;
}

function getHorizontalMultivecViewId(viewConf) {
    let viewId = null;
    let trackId = null;
    if(viewConf && viewConf.views && viewConf.views.length > 0) {
        for(let viewObj of viewConf.views) {
            viewId = viewObj.uid;
            trackId = getHorizontalMultivecTrackId(viewObj);
            if(viewId && trackId) {
                break;
            }
        }
    }
    return [viewId, trackId];
}


/**
 * @component Cistrome HiGlass Wrapper 
 */
export default function CistromeHGW(props) {

    const { viewConfig } = props;

    const hgRef = useRef();
    
    useEffect(() => {
        hgRef.current.api.on('location', (d) => {
            PubSub.publish(GLOBAL_X_RANGE, d.xRange);
            PubSub.publish(GLOBAL_Y_RANGE, d.yRange);
        });

        hgRef.current.api.on('viewConfig', (newViewConfigString) => {
            const newViewConfig = JSON.parse(newViewConfigString);
            const [viewId, trackId] = getHorizontalMultivecViewId(newViewConfig);
            try {
                const trackObj = hgRef.current.api.getTrackObject(viewId, trackId);
                PubSub.publish(TRACK_ROW_INFO, trackObj.tilesetInfo.row_infos);
                PubSub.publish(TRACK_POSITION, trackObj.position);
                PubSub.publish(TRACK_DIMENSIONS, trackObj.dimensions);
            } catch(e) {
    
            }
        });
    });

    return (
        <div className="cistrome-hgw">
            <HiGlassComponent
                viewConfig={viewConfig}
                options={hgOptions}
                zoomFixed={false}
                ref={hgRef}
            />
            <CistromeGroupLabels />
        </div>
    );
}