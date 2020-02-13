import React, { useRef, useEffect, useState, useMemo, useCallback } from 'react';

import { HiGlassComponent } from 'higlass';
import register from 'higlass-register';
import StackedBarTrack from 'higlass-multivec/es/StackedBarTrack.js';

import PubSub from 'pubsub-js';
import { EVENT } from './constants.js';
import TrackWrapper from './TrackWrapper.js';
import Tooltip from './Tooltip.js';

import { processWrapperOptions, DEFAULT_OPTIONS_KEY, updateRowSortOptions } from './utils/options.js';
import { 
    getHMTrackIdsFromViewConfig, 
    getSiblingVPHTrackIdsFromViewConfig,
    updateViewConfigOnSelectGenomicInterval
} from './utils/viewconf.js';

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
    const drawRef = useRef({});

    const [options, setOptions] = useState({});
    const [trackIds, setTrackIds] = useState([]);
    const [siblingTrackIds, setSiblingTrackIds] = useState({});
    const [selectedRows, setSelectedRows] = useState({});
    const [highlitRows, setHighlitRows] = useState({});

    /*
     * Function to call when the view config has changed.
     * Updates the array of `horizontal-multivec` track IDs.
     * Updates the mapping of `horizontal-multivec` track IDs 
     * to sibling `viewport-projection-horizontal` track IDs.
     */
    const onViewConfig = useCallback((newViewConfig) => {
        const newTrackIds = getHMTrackIdsFromViewConfig(newViewConfig);
        const newSiblingTrackIds = {};
        const newSelectedRows = {};
        const newHighlitRows = {};
        for(let trackId of newTrackIds) {
            // Each trackId is actually an array `[viewId, trackId]`, which is why we want trackId[1].
            newSiblingTrackIds[trackId[1]] = getSiblingVPHTrackIdsFromViewConfig(newViewConfig, trackId[1]);
            newSelectedRows[trackId[1]] = null;
            newHighlitRows[trackId[1]] = null;
        }
        setTrackIds(newTrackIds);
        setSiblingTrackIds(newSiblingTrackIds);
        setSelectedRows(newSelectedRows);
        setHighlitRows(newHighlitRows);
    }, []);

    /*
     * Function to get a track object from the higlass API.
     */
    const getTrackObject = useCallback((viewId, trackId) => {
        try {
            return hgRef.current.api.getTrackObject(viewId, trackId);
        } catch(e) {
            return null;
        }
    }, []);

    /*
     * Function to get an options object for a particular track.
     */
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
    
    /*
     * Function for child components to call to "register" their draw functions.
     * These draw functions will be called when the component is exported to SVG.
     */
    const register = useCallback((key, drawFunction) => {
        drawRef.current[key] = drawFunction;
    }, [drawRef]);

    useEffect(() => {
        setOptions(processWrapperOptions(optionsRaw));
    }, [optionsRaw]);
    
    // Change options by interactions.
    useEffect(() => {
        const sortToken = PubSub.subscribe(EVENT.SORT, (msg, data) => {
            const newOptions = processWrapperOptions(updateRowSortOptions(optionsRaw, data));
            setOptions(newOptions);
        });
        return () => {
            PubSub.unsubscribe(sortToken);
        };
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
                    selectedRows={selectedRows[trackId] ? selectedRows[trackId] : null}
                    highlitRows={highlitRows[trackId] ? highlitRows[trackId] : null}
                    onSelectGenomicInterval={() => {
                        const currViewConfig = hgRef.current.api.getViewConfig();
                        const newViewConfig = updateViewConfigOnSelectGenomicInterval(currViewConfig, viewId, trackId);
                        hgRef.current.api.setViewConfig(newViewConfig).then(() => {
                            onViewConfig(newViewConfig);
                        });
                    }}
                    register={register}
                />
            ))}
            <Tooltip />
        </div>
    );
}