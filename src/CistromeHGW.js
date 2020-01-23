import React, { useRef, useEffect, useState, useMemo } from 'react';

import { HiGlassComponent } from 'higlass';
import register from 'higlass-register';
import StackedBarTrack from 'higlass-multivec/es/StackedBarTrack.js';

import CistromeGroupLabels from './CistromeGroupLabels.js';

import 'higlass/dist/hglib.css';
import './CistromeHGW.scss';

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

    const [x1, setX1] = useState(0);
    const [x0, setX0] = useState(0);
    const [y, setY] = useState(0);
    const [height, setHeight] = useState(30);

    const [rowNames, setRowNames] = useState([]);
    
    useEffect(() => {
        hgRef.current.api.on('location', (d) => {
            setX1(d.xRange[1]);
        });
        hgRef.current.api.on('viewConfig', (newViewConfigString) => {
            try {
                const newViewConfig = JSON.parse(newViewConfigString);
                const [viewId, trackId] = getHorizontalMultivecViewId(newViewConfig);
                const trackObj = hgRef.current.api.getTrackObject(viewId, trackId);

                setX0(trackObj.position[0]);
                setY(trackObj.position[1]);
                setHeight(trackObj.dimensions[1]);
                setRowNames(trackObj.tilesetInfo.row_infos.map(d => d.split("\t")));

            } catch(e) {
                console.log(e);
            }
        });

        return () => {
            hgRef.current.api.off('location');
            hgRef.current.api.off('viewConfig');
        };
    }, [hgRef, setX0, setX1, setY, setHeight, setRowNames]);

    const hgComponent = useMemo(() => {
        console.log("HiGlassComponent.render");
        return (
        <HiGlassComponent
            viewConfig={viewConfig}
            options={hgOptions}
            zoomFixed={false}
            ref={hgRef}
        />
    )}, [viewConfig, hgOptions]);

    console.log("CistromeHGW.render");

    return (
        <div className="cistrome-hgw">
            {hgComponent}
            <CistromeGroupLabels 
                rowNames={rowNames}
                x={x0 + x1}
                y={y}
                height={height}
            />
        </div>
    );
}