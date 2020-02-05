import React, { useRef, useEffect, useState, useMemo, useCallback } from 'react';

import { HiGlassComponent } from 'higlass';
import register from 'higlass-register';
import StackedBarTrack from 'higlass-multivec/es/StackedBarTrack.js';

import TrackWrapper from './TrackWrapper.js';
import Tooltip from './Tooltip.js';

import { processWrapperOptions, DEFAULT_OPTIONS_KEY } from './utils-options.js';
import { 
    getHMTrackIdsFromViewConfig, 
    getSiblingProjectionTracksFromViewConfig,
    updateViewConfigOnSelectGenomicInterval
} from './utils-viewconf.js';

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
 * Cistrome HiGlass Wrapper, a React component that wraps around HiGlass to provide visualization features for cistrome data.
 * @prop {object} viewConfig A HiGlass viewConfig object.
 * @prop {(object|object[])} options Options for the wrapper component, for positioning child components.
 * @example
 * <CistromeHGW
 *  viewConfig={higlassViewConfig}
 *  options={{
 *      rowInfoPosition: "right",
 *      rowLinkPosition: "left"
 *  }}
 * />
 */
export default function CistromeHGW(props) {

    const { viewConfig, options: optionsRaw } = props;

    const hgRef = useRef();

    const [options, setOptions] = useState({});

    // Array of horizontal-multivec track IDs.
    const [trackIds, setTrackIds] = useState([]);
    // Mapping of horizontal-multivec track IDs to arrays of viewport-projection sibling track IDs.
    const [siblingTrackIds, setSiblingTrackIds] = useState({});

    const onViewConfig = useCallback((newViewConfig) => {
        const newTrackIds = getHMTrackIdsFromViewConfig(newViewConfig);
        const newSiblingTrackIds = {};
        for(let trackId of newTrackIds) {
            newSiblingTrackIds[trackId[1]] = getSiblingProjectionTracksFromViewConfig(newViewConfig, trackId[1]);
        }
        setTrackIds(newTrackIds);
        setSiblingTrackIds(newSiblingTrackIds);
    }, []);

    const getTrackObject = useCallback((viewId, trackId) => {
        try {
            return hgRef.current.api.getTrackObject(viewId, trackId);
        } catch(e) {
            return null;
        }
    }, []);

    const getTrackWrapperOptions = useCallback((viewId, trackId) => {
        if(options[viewId]) {
            if(options[viewId][trackId]) {
                return options[viewId][trackId];
            } else {
                return options[viewId][DEFAULT_OPTIONS_KEY];
            }
        } else {
            return options[DEFAULT_OPTIONS_KEY];
        }
    }, [options]);

    useEffect(() => {
        setOptions(processWrapperOptions(optionsRaw));
    }, [optionsRaw]);
    
    useEffect(() => {
        hgRef.current.api.on('viewConfig', (newViewConfigString) => {
            const newViewConfig = JSON.parse(newViewConfigString);
            onViewConfig(newViewConfig);
        });         

        return () => {
            hgRef.current.api.off('viewConfig');
        };
    }, [hgRef]);

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
            {trackIds.map(([viewId, trackId, combinedTrackId], i) => (
                <TrackWrapper
                    key={i}
                    options={getTrackWrapperOptions(viewId, trackId)}
                    multivecTrack={getTrackObject(viewId, trackId)}
                    combinedTrack={(combinedTrackId ? getTrackObject(viewId, combinedTrackId) : null)}
                    siblingTracks={siblingTrackIds[trackId] ? siblingTrackIds[trackId].map(d => getTrackObject(viewId, d[1])) : []}
                    onSelectGenomicInterval={() => {
                        const currViewConfig = hgRef.current.api.getViewConfig();
                        const newViewConfig = updateViewConfigOnSelectGenomicInterval(currViewConfig, viewId, trackId);
                        hgRef.current.api.setViewConfig(newViewConfig).then(() => {
                            onViewConfig(newViewConfig);
                        });
                    }}
                />
            ))}
            <Tooltip />
        </div>
    );
}