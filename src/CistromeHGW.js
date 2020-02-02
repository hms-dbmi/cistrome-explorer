import React, { useRef, useEffect, useState, useMemo } from 'react';

import { HiGlassComponent } from 'higlass';
import register from 'higlass-register';
import StackedBarTrack from 'higlass-multivec/es/StackedBarTrack.js';

import TrackWrapper from './TrackWrapper.js';
import Tooltip from './Tooltip.js';

import { processWrapperOptions, DEFAULT_OPTIONS_KEY } from './utils-options.js';
import { getHorizontalMultivecTracksFromViewConfig } from './utils-viewconf.js';

import './CistromeHGW.scss';

register({
    name: 'StackedBarTrack',
    track: StackedBarTrack,
    config: StackedBarTrack.config,
});

const hgOptionsBase = {
    bounded: true,
    pixelPreciseMarginPadding: true,
    containerPaddingX: 0,
    containerPaddingY: 0,
    sizeMode: 'default'
};


/**
 * @component Cistrome HiGlass Wrapper 
 */
export default function CistromeHGW(props) {

    const { viewConfig, options: optionsRaw } = props;

    const hgRef = useRef();

    const [x0, setX0] = useState(0);
    const [y0, setY0] = useState(0);
    const [options, setOptions] = useState({});
    const [trackIds, setTrackIds] = useState([]);

    function onViewConfig(newViewConfig) {
        const newTrackIds = getHorizontalMultivecTracksFromViewConfig(newViewConfig);
        setTrackIds(newTrackIds);
    }

    function getTrackObject(viewId, trackId) {
        try {
            return hgRef.current.api.getTrackObject(viewId, trackId);
        } catch(e) {
            return null;
        }
    }

    function getTrackWrapperOptions(viewId, trackId) {
        if(options[viewId]) {
            if(options[viewId][trackId]) {
                return options[viewId][trackId];
            } else {
                return options[viewId][DEFAULT_OPTIONS_KEY];
            }
        } else {
            return options[DEFAULT_OPTIONS_KEY];
        }
    }

    useEffect(() => {
        setOptions(processWrapperOptions(optionsRaw));
    }, [optionsRaw]);
    
    useEffect(() => {
        hgRef.current.api.on('location', (d) => {
            setX0(d.xRange[1]);
            setY0(d.yRange[1]);
        });
        hgRef.current.api.on('viewConfig', (newViewConfigString) => {
            const newViewConfig = JSON.parse(newViewConfigString);
            onViewConfig(newViewConfig);
        });

        return () => {
            hgRef.current.api.off('location');
            hgRef.current.api.off('viewConfig');
        };
    }, [hgRef, setX0, setY0]);

    const hgComponent = useMemo(() => {
        const hgOptions = {
            ...hgOptionsBase,
            onViewConfLoaded: () => {
                onViewConfig(viewConfig);
            }
        };

        console.log("HiGlassComponent.render");
        return (
            <HiGlassComponent
                viewConfig={viewConfig}
                options={hgOptions}
                zoomFixed={false}
                ref={hgRef}
            />
        );
    }, [viewConfig]);
    

    console.log("CistromeHGW.render");
    return (
        <div className="cistrome-hgw">
            {hgComponent}
            {trackIds.map(([viewId, trackId], i) => (
                <TrackWrapper 
                    key={i}
                    options={getTrackWrapperOptions(viewId, trackId)}
                    track={getTrackObject(viewId, trackId)}
                    x0={x0}
                />
            ))}
            <Tooltip />
        </div>
    );
}