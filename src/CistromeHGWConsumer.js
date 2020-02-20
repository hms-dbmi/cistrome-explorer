import React, { useRef, useEffect, useState, useMemo, useCallback, useContext } from 'react';
import isEqual from 'lodash/isEqual';

import { HiGlassComponent } from 'higlass';
import higlassRegister from 'higlass-register';
import StackedBarTrack from 'higlass-multivec/es/StackedBarTrack.js';

import { InfoContext, ACTION } from './utils/contexts.js';
import { selectRows, highlightRowsFromSearch } from './utils/select-rows.js';
import TrackWrapper from './TrackWrapper.js';
import Tooltip from './Tooltip.js';

import { 
    processWrapperOptions, 
    getTrackWrapperOptions, 
    updateGlobalOptionsWithKey 
} from './utils/options.js';
import { 
    getHMTrackIdsFromViewConfig, 
    getSiblingVPHTrackIdsFromViewConfig,
    updateViewConfigOnSelectGenomicInterval,
    updateViewConfigOnSelectRowsByTrack,
    getHMSelectedRowsFromViewConfig
} from './utils/viewconf.js';

import './CistromeHGWConsumer.scss';

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
 * @prop {object} viewConfig A HiGlass viewConfig object.
 * @prop {(object|object[])} options Options for the wrapper component.
 */
export default function CistromeHGWConsumer(props) {

    const { viewConfig, options: initOptionsRaw } = props;

    const hgRef = useRef();
    const drawRef = useRef({});

    const [optionsRaw, setOptionsRaw] = useState(initOptionsRaw);
    const [options, setOptions] = useState(processWrapperOptions(initOptionsRaw));
    const [trackInfos, setTrackInfos] = useState([]);
    const [siblingTrackInfos, setSiblingTrackInfos] = useState({});
    
    const context = useContext(InfoContext);

    /*
     * Function to call when the view config has changed.
     * Updates the array of `horizontal-multivec` track IDs.
     * Updates the mapping of `horizontal-multivec` track IDs 
     * to sibling `viewport-projection-horizontal` track IDs.
     */
    const onViewConfig = useCallback((newViewConfig) => {
        const newTrackInfos = getHMTrackIdsFromViewConfig(newViewConfig);
        const newSiblingTrackInfos = {};
        for(let trackInfo of newTrackInfos) {
            newSiblingTrackInfos[trackInfo.trackId] = getSiblingVPHTrackIdsFromViewConfig(newViewConfig, trackInfo.trackId);
            const newSelectedRows = getHMSelectedRowsFromViewConfig(newViewConfig, trackInfo.viewId, trackInfo.trackId);
            if(
                !context.state[trackInfo.viewId] 
                || !context.state[trackInfo.viewId][trackInfo.trackId] 
                || !isEqual(newSelectedRows, context.state[trackInfo.viewId][trackInfo.trackId].selectedRows)
            ) {
                context.dispatch({
                    type: ACTION.SELECT_ROWS,
                    viewId: trackInfo.viewId,
                    trackId: trackInfo.trackId,
                    selectedRows: newSelectedRows
                });
            }
        }
        setTrackInfos(newTrackInfos);
        setSiblingTrackInfos(newSiblingTrackInfos);
    }, []);

    // Function to get a track object from the higlass API.
    const getTrackObject = useCallback((viewId, trackId) => {
        try {
            return hgRef.current.api.getTrackObject(viewId, trackId);
        } catch(e) {
            return null;
        }
    }, []);

    // Function to get view location information from the higlass API.
    const getViewLocationObject = useCallback((viewId) => {
        try {
            const location = hgRef.current.api.getLocation(viewId);
            // TODO: remove the hard-coded xRange and yRange values below,
            // once the hgc.api.getLocation() function has been updated to return them.
            return {
                ...location,
                xRange: [0, 1762],
                yRange: [0, 791]
            };
        } catch(e) {
            return null;
        }
    }, []);

    const setTrackSelectedRows = useCallback((viewId, trackId, selectedRows) => {
        const currViewConfig = hgRef.current.api.getViewConfig();
        const newViewConfig = updateViewConfigOnSelectRowsByTrack(currViewConfig, selectedRows, viewId, trackId);
        hgRef.current.api.setViewConfig(newViewConfig).then(() => {
            onViewConfig(newViewConfig);
        });
    }, [hgRef]);
    
    const setHighlitRows = useCallback((viewId, trackId, highlitRows) => {
        context.dispatch({
            type: ACTION.HIGHLIGHT_ROWS,
            viewId,
            trackId,
            highlitRows
        });
    });

    // Function for child components to call to "register" their draw functions.
    const drawRegister = useCallback((key, drawFunction) => {
        drawRef.current[key] = drawFunction;
    }, [drawRef]);

    // Callback function for sorting.
    const onSortRows = useCallback((viewId, trackId, field, type, order) => {
        const newRowSort = [ { field, type, order } ];
        const newOptionsRaw = updateGlobalOptionsWithKey(optionsRaw, newRowSort, "rowSort", { isReplace: true });
        const newOptions = processWrapperOptions(newOptionsRaw);
        
        const trackOptions = getTrackWrapperOptions(newOptions, viewId, trackId);
        const newSelectedRows = selectRows(context.state[viewId][trackId].rowInfo, trackOptions);

        setTrackSelectedRows(viewId, trackId, newSelectedRows);
        setOptionsRaw(newOptionsRaw);
        setOptions(newOptions);
    });

    // Callback function for filtering.
    const onFilter = useCallback((viewId, trackId, field, type, contains) => {
        const isResetFilter = field === undefined;
        const newRowFilter = isResetFilter ? [] : { field, type, contains };
        let newOptionsRaw = updateGlobalOptionsWithKey(optionsRaw, newRowFilter, "rowFilter", { isReplace: isResetFilter });
        newOptionsRaw = updateGlobalOptionsWithKey(newOptionsRaw, undefined, "rowHighlight", { isReplace: true });    // Reset highlight as well.
        const newOptions = processWrapperOptions(newOptionsRaw);

        const trackOptions = getTrackWrapperOptions(newOptions, viewId, trackId);
        const newSelectedRows = selectRows(context.state[viewId][trackId].rowInfo, trackOptions);

        setHighlitRows(viewId, trackId, undefined);
        setTrackSelectedRows(viewId, trackId, newSelectedRows);
        setOptionsRaw(newOptionsRaw);
        setOptions(newOptions);
    });

    // Callback function for searching.
    const onSearchRows = useCallback((viewId, trackId, field, type, contains) => {
        const newRowHighlight = { field, type, contains };
        const newOptionsRaw = updateGlobalOptionsWithKey(optionsRaw, newRowHighlight, "rowHighlight", { isReplace: true });
        const newOptions = processWrapperOptions(newOptionsRaw);

        // Highlighting options are specified only in the wrapper options.
        const newHighlitRows = highlightRowsFromSearch(context.state[viewId][trackId].rowInfo, field, type, contains);
        setHighlitRows(viewId, trackId, newHighlitRows);
        setOptionsRaw(newOptionsRaw);
        setOptions(newOptions);
    });

    // Callback function for searching.
    const onSelectGenomicInterval = useCallback((viewId, trackId, intervalCoordinates) => {
        const newOptionsRaw = updateGlobalOptionsWithKey(optionsRaw, intervalCoordinates, "colSelect", { isReplace: false });
        console.log(newOptionsRaw);
        const newOptions = processWrapperOptions(newOptionsRaw);

        // Highlighting options are specified only in the wrapper options.
        setOptionsRaw(newOptionsRaw);
        setOptions(newOptions);
    });

    // Do initial processing of the options prop.
    useEffect(() => {
        setOptionsRaw(initOptionsRaw);
        setOptions(processWrapperOptions(initOptionsRaw));
    }, [initOptionsRaw]);

    // Listen for higlass view config changes.
    useEffect(() => {
        hgRef.current.api.on('viewConfig', (newViewConfigString) => {
            const newViewConfig = JSON.parse(newViewConfigString);
            onViewConfig(newViewConfig);
        });         

        return () => hgRef.current.api.off('viewConfig');
    }, [hgRef]);

    // We only want to render HiGlass once.
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
            {trackInfos.map(({ viewId, trackId, trackTilesetId, combinedTrackId }, i) => (
                <TrackWrapper
                    key={i}
                    options={getTrackWrapperOptions(options, viewId, trackId)}
                    viewLocation={getViewLocationObject(viewId)}
                    multivecTrack={getTrackObject(viewId, trackId)}
                    multivecTrackViewId={viewId}
                    multivecTrackTrackId={trackId}
                    multivecTrackTilesetId={trackTilesetId}
                    combinedTrack={(combinedTrackId ? getTrackObject(viewId, combinedTrackId) : null)}
                    siblingTracks={siblingTrackInfos[trackId] ? siblingTrackInfos[trackId].map(d => getTrackObject(viewId, d.trackId)) : []}
                    onSelectGenomicInterval={(interval) => {
                        onSelectGenomicInterval(viewId, trackId, interval);
                    }}
                    onSortRows={(field, type, order) => {
                        onSortRows(viewId, trackId, field, type, order);
                    }}
                    onSearchRows={(field, type, contains) => {
                        onSearchRows(viewId, trackId, field, type, contains);
                    }}
                    onFilter={(field, type, contains) => {
                        onFilter(viewId, trackId, field, type, contains);
                    }}
                    drawRegister={drawRegister}
                />
            ))}
            <Tooltip />
        </div>
    );
}