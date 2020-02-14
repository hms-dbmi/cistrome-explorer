import React, { useRef, useEffect, useState, useMemo, useCallback, useContext } from 'react';
import PubSub from 'pubsub-js';

import { HiGlassComponent } from 'higlass';
import higlassRegister from 'higlass-register';
import StackedBarTrack from 'higlass-multivec/es/StackedBarTrack.js';

import { EVENT } from './utils/constants.js';
import { InfoContext, ACTION } from './utils/contexts.js';
import { selectRows } from './utils/select-rows.js';
import TrackWrapper from './TrackWrapper.js';
import Tooltip from './Tooltip.js';

import { processWrapperOptions, DEFAULT_OPTIONS_KEY, updateRowSortOptions } from './utils/options.js';
import { 
    getHMTrackIdsFromViewConfig, 
    getSiblingVPHTrackIdsFromViewConfig,
    updateViewConfigOnSelectGenomicInterval,
    updateViewConfigOnSelectRowsByTrack,
    getHMSelectedRowsFromViewConfig
} from './utils/viewconf.js';

import './CistromeHGW.scss';

higlassRegister({
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
 * CistromeHGW passes its props through, and wraps this component with the context provider.
 */
export default function CistromeHGWConsumer(props) {

    const { viewConfig, options: optionsRaw } = props;

    const hgRef = useRef();
    const drawRef = useRef({});

    const [options, setOptions] = useState({});
    const [trackIds, setTrackIds] = useState([]);
    const [siblingTrackIds, setSiblingTrackIds] = useState({});

    const context = useContext(InfoContext);

    /*
     * Function to call when the view config has changed.
     * Updates the array of `horizontal-multivec` track IDs.
     * Updates the mapping of `horizontal-multivec` track IDs 
     * to sibling `viewport-projection-horizontal` track IDs.
     */
    const onViewConfig = useCallback((newViewConfig) => {
        const newTrackIds = getHMTrackIdsFromViewConfig(newViewConfig);
        const newSiblingTrackIds = {};
        for(let trackId of newTrackIds) {
            // Each trackId is actually an object `{ viewId, trackId }`, which is why we want trackId.trackId.
            newSiblingTrackIds[trackId.trackId] = getSiblingVPHTrackIdsFromViewConfig(newViewConfig, trackId.trackId);
            context.dispatch({
                type: ACTION.SELECT_ROWS,
                trackId: trackId.trackId,
                selectedRows: getHMSelectedRowsFromViewConfig(newViewConfig, trackId.viewId, trackId.trackId)
            });

            //newHighlitRows[trackId.trackId] = null;
        }
        setTrackIds(newTrackIds);
        setSiblingTrackIds(newSiblingTrackIds);

        // TODO: dispatch to context
        //setSelectedRows(newSelectedRows);
        //setHighlitRows(newHighlitRows);
    }, []);

    // Function to get a track object from the higlass API.
    const getTrackObject = useCallback((viewId, trackId) => {
        try {
            return hgRef.current.api.getTrackObject(viewId, trackId);
        } catch(e) {
            return null;
        }
    }, []);

    // Function to get an options object for a particular track.
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

    const setTrackSelectedRows = useCallback((selectedRows, viewId, trackId) => {
        const currViewConfig = hgRef.current.api.getViewConfig();
        const newViewConfig = updateViewConfigOnSelectRowsByTrack(currViewConfig, selectedRows, viewId, trackId);
        console.log(newViewConfig);
        hgRef.current.api.setViewConfig(newViewConfig).then(() => {
            onViewConfig(newViewConfig);
        });
    }, [hgRef]);
    

    // Function for child components to call to "register" their draw functions.
    const drawRegister = useCallback((key, drawFunction) => {
        drawRef.current[key] = drawFunction;
    }, [drawRef]);

    // Process the options prop (initially, upon sorting, etc).
    useEffect(() => {
        setOptions(processWrapperOptions(optionsRaw));

        const sortToken = PubSub.subscribe(EVENT.SORT, (msg, data) => {
            const newOptionsRaw = updateRowSortOptions(optionsRaw, data);
            setOptions(processWrapperOptions(newOptionsRaw));
            // TODO: generalize. need to obtain the viewId,trackId from the sort event data.
            const newSelectedRows = selectRows(context.state['cistrome-view-2'], {}, null);
            setTrackSelectedRows(newSelectedRows, 'cistrome-view-2', 'cistrome-track-2');
        });

        return () => PubSub.unsubscribe(sortToken);
    }, [optionsRaw]);

    // Listen for higlass view config changes.
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
    
    console.log("CistromeHGWConsumer.render");
    return (
        <div className="cistrome-hgw">
            {hgComponent}
            {trackIds.map(({ viewId, trackId, trackTilesetId, combinedTrackId }, i) => (
                <TrackWrapper
                    key={i}
                    options={getTrackWrapperOptions(viewId, trackId)}
                    multivecTrack={getTrackObject(viewId, trackId)}
                    multivecTrackTilesetId={trackTilesetId}
                    combinedTrack={(combinedTrackId ? getTrackObject(viewId, combinedTrackId) : null)}
                    siblingTracks={siblingTrackIds[trackId] ? siblingTrackIds[trackId].map(d => getTrackObject(viewId, d.trackId)) : []}
                    onSelectGenomicInterval={() => {
                        const currViewConfig = hgRef.current.api.getViewConfig();
                        const newViewConfig = updateViewConfigOnSelectGenomicInterval(currViewConfig, viewId, trackId);
                        hgRef.current.api.setViewConfig(newViewConfig).then(() => {
                            onViewConfig(newViewConfig);
                        });
                    }}
                    drawRegister={drawRegister}
                />
            ))}
            <Tooltip />
        </div>
    );
}