import React, { useRef, useEffect, useState, useMemo, useCallback, useContext } from 'react';
import isEqual from 'lodash/isEqual';

import { HiGlassComponent } from 'higlass';
import higlassRegister from 'higlass-register';
import StackedBarTrack from 'higlass-multivec/es/StackedBarTrack.js';

import { InfoContext, ACTION } from './utils/contexts.js';
import { selectRows, highlightRowsFromSearch } from './utils/select-rows.js';
import TrackWrapper from './TrackWrapper.js';
import ViewWrapper from './ViewWrapper.js';
import Tooltip from './Tooltip.js';
import ContextMenu, { destroyContextMenu } from './ContextMenu.js';

import { 
    DEFAULT_OPTIONS_KEY,
    processWrapperOptions, 
    getTrackWrapperOptions,
    updateWrapperOptions,
    addTrackWrapperOptions
} from './utils/options.js';
import { 
    getHMTrackIdsFromViewConfig, 
    getSiblingVPHTrackIdsFromViewConfig,
    updateViewConfigOnSelectGenomicInterval,
    updateViewConfigOnSelectRowsByTrack,
    getHMSelectedRowsFromViewConfig,
    getUniqueViewOrTrackId,
    getTrackDefFromViewConfig,
    addTrackDefToViewConfig,
    getAllViewAndTrackPairs,
    removeViewportFromViewConfig
} from './utils/viewconf.js';
import { wrapSvg } from './utils/wrap-svg.js';

import './CistromeHGWConsumer.scss';
import cloneDeep from 'lodash/cloneDeep';

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
 * @prop {function} onViewConfigChange A function to call upon change of the HiGlass view config. Optional.
 */
export default function CistromeHGWConsumer(props) {

    const {
        viewConfig,
        options: optionsRaw,
        onViewConfigChange: onViewConfigChangeProp
    } = props;

    const hgRef = useRef();
    const drawRef = useRef({});

    const [options, setOptions] = useState({});
    const [muiltivecTrackIds, setMultivecTrackIds] = useState([]);
    const [viewportTrackIds, setViewportTrackIds] = useState({});
    
    const context = useContext(InfoContext);

    // Set initial sorting, filtering, and highlighting.
    const onMetadataLoad = useCallback((viewId, trackId) => {
        if(!context.state[viewId] || !context.state[viewId][trackId]) {
            return;
        }
        const trackOptions = getTrackWrapperOptions(options, viewId, trackId);
        const rowInfo = context.state[viewId][trackId].rowInfo;
        
        // Filter and sort.
        const newSelectedRows = selectRows(rowInfo, trackOptions);
        setTrackSelectedRows(viewId, trackId, newSelectedRows);
        
        // Highlight.
        if(trackOptions.rowHighlight) {
            const { field, type, contains } = trackOptions.rowHighlight;
            const newHighlitRows = highlightRowsFromSearch(rowInfo, field, type, contains);
            setHighlitRows(viewId, trackId, newHighlitRows);
        }
    }, [options]);

    /*
     * Function to call when the view config has changed.
     * Updates the array of `horizontal-multivec` track IDs.
     * Updates the mapping of `horizontal-multivec` track IDs 
     * to sibling `viewport-projection-horizontal` track IDs.
     */
    const onViewConfig = useCallback((newViewConfig) => {
        const newTrackIds = getHMTrackIdsFromViewConfig(newViewConfig);
        
        // Add viewport projection horizontal track IDs for each view.
        const newViewportTrackIds = {};
        Array.from(new Set(newTrackIds.map(d => d.viewId))).forEach((viewId) => {    
            newViewportTrackIds[viewId] = getSiblingVPHTrackIdsFromViewConfig(newViewConfig, viewId);
        });

        for(let trackIds of newTrackIds) {
            const newSelectedRows = getHMSelectedRowsFromViewConfig(newViewConfig, trackIds.viewId, trackIds.trackId);
            if(
                !context.state[trackIds.viewId] 
                || !context.state[trackIds.viewId][trackIds.trackId] 
                || !isEqual(newSelectedRows, context.state[trackIds.viewId][trackIds.trackId].selectedRows)
            ) {
                context.dispatch({
                    type: ACTION.SELECT_ROWS,
                    viewId: trackIds.viewId,
                    trackId: trackIds.trackId,
                    selectedRows: newSelectedRows
                });
            }
        }
        setMultivecTrackIds(newTrackIds);
        setViewportTrackIds(newViewportTrackIds);
    }, []);

    // Function to get a track object from the higlass API.
    const getTrackObject = useCallback((viewId, trackId) => {
        try {
            return hgRef.current.api.getTrackObject(viewId, trackId);
        } catch(e) {
            return null;
        }
    }, [hgRef]);

    // Custom function to get the boundingbox of a view in higlass,
    // i.e., { left, top, width, height }
    const getViewBoundingBox = useCallback((viewId) => {
        try {
            const viewConfig = hgRef.current.api.getViewConfig();
            const trackIds = getAllViewAndTrackPairs(viewConfig, { onlyHorizontalMultivec: false }).filter(d => d.viewId === viewId);
            
            let viewLeft = null, viewTop = null;
            let viewWidth = null, viewBottom = null;
            trackIds.forEach(({ trackId }) => {
                const track = hgRef.current.api.getTrackObject(viewId, trackId);
                const [left, top] = track.position;
                const [width, height] = track.dimensions;

                if(!viewLeft) viewLeft = left;
                if(!viewWidth) viewWidth = width;
                if(!viewTop || top < viewTop) viewTop = top;
                if(!viewBottom || top + height > viewBottom) viewBottom = top + height;
            });
            return { left: viewLeft, top: viewTop, width: viewWidth, height: viewBottom - viewTop };
        } catch(e) {
            return null;
        }
    }, [hgRef]);

    const addNewTrack = useCallback((trackDef, viewId, position) => {
        const currViewConfig = hgRef.current.api.getViewConfig();
        const newViewConfig = addTrackDefToViewConfig(currViewConfig, trackDef, viewId, position);
        hgRef.current.api.setViewConfig(newViewConfig).then(() => {
            onViewConfig(newViewConfig);
        });
    }, [hgRef]);

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
    }, [hgRef]);

    // Function for child components to call to "register" their draw functions.
    const drawRegister = useCallback((key, draw, options) => {
        drawRef.current[key] = { draw, options };
    }, [drawRef]);

    // Listen for the `createSVG` event.
    useEffect(() => {
        hgRef.current.api.on('createSVG', (svg) => {
            return wrapSvg(svg, drawRef.current);
        });
        return () => hgRef.current.api.off('createSVG');
    }, [hgRef, drawRef]);
    
    // Callback function for sorting.
    const onSortRows = useCallback((viewId, trackId, field, type, order) => {
        const newRowSort = [ { field, type, order } ];
        const newOptions = updateWrapperOptions(options, newRowSort, "rowSort", viewId, trackId, { isReplace: true });
        
        const trackOptions = getTrackWrapperOptions(newOptions, viewId, trackId);
        const newSelectedRows = selectRows(context.state[viewId][trackId].rowInfo, trackOptions);
        
        setTrackSelectedRows(viewId, trackId, newSelectedRows);
        setOptions(newOptions);
    }, [options]);

    // Callback function for searching and highlighting.
    const onSearchRows = useCallback((viewId, trackId, field, type, contains) => {
        // TODO: Better handle with multiple fields & quantitative fields.
        if(Array.isArray(field)) return;
        const newRowHighlight = { field, type, contains };
        const newOptions = updateWrapperOptions(options, newRowHighlight, "rowHighlight", viewId, trackId, { isReplace: true });

        // Highlighting options are specified only in the wrapper options.
        const newHighlitRows = highlightRowsFromSearch(context.state[viewId][trackId].rowInfo, field, type, contains);
        setHighlitRows(viewId, trackId, newHighlitRows);
        setOptions(newOptions);
    }, [options]);

    // Callback function for filtering.
    const onFilterRows = useCallback((viewId, trackId, field, type, condition) => {
        const isResetFilter = field === undefined;
        const filterKey = type === "nominal" ? "contains" : "range";
        const newRowFilter = isResetFilter ? [] : { field, type, [filterKey]: condition };
        let newOptions = updateWrapperOptions(options, newRowFilter, "rowFilter", viewId, trackId, { isReplace: isResetFilter });
        newOptions = updateWrapperOptions(newOptions, undefined, "rowHighlight", viewId, trackId, { isReplace: true });

        const trackOptions = getTrackWrapperOptions(newOptions, viewId, trackId);
        const newSelectedRows = selectRows(context.state[viewId][trackId].rowInfo, trackOptions);

        setHighlitRows(viewId, trackId, undefined);
        setTrackSelectedRows(viewId, trackId, newSelectedRows);
        setOptions(newOptions);
    }, [options]);

    // Callback function for adding a track.
    const onAddTrack = useCallback((viewId, trackId, field, type, contains, position) => {
        if(viewId === DEFAULT_OPTIONS_KEY || trackId === DEFAULT_OPTIONS_KEY) {
            console.log("A view or track ID is a default ID.");
            return;
        }
        const newRowFilter = { field, type, contains };
        const currViewConfig = hgRef.current.api.getViewConfig();
        const newTrackId = getUniqueViewOrTrackId(currViewConfig, { 
            baseId: trackId, 
            idKey: "trackId", 
            interfix: "detail-view"
        });
        
        // Add options for the new track.
        const trackOptionsCopy = cloneDeep(getTrackWrapperOptions(options, viewId, trackId));
        let newOptions = addTrackWrapperOptions(options, trackOptionsCopy, viewId, newTrackId);

        // Add filter options.
        newOptions = updateWrapperOptions(newOptions, newRowFilter, "rowFilter", viewId, newTrackId, { isReplace: false });

        // Get new selectedRows for the new viewConfig.
        const newTrackOptions = getTrackWrapperOptions(newOptions, viewId, newTrackId);
        const newSelectedRows = selectRows(context.state[viewId][trackId].rowInfo, newTrackOptions);   // Use original rowInfo.

        // Get new viewConfig with new selectedRows.
        let newTrackDef = getTrackDefFromViewConfig(currViewConfig, viewId, trackId);
        newTrackDef = {
            ...newTrackDef,
            height: 200,
            uid: newTrackId,
            options: {
                ...newTrackDef.options,
                selectRows: newSelectedRows
            }
        }
        addNewTrack(newTrackDef, viewId, position);
        setOptions(newOptions);
    }, [options]);

    // Do initial processing of the options prop.
    useEffect(() => {
        setOptions(processWrapperOptions(optionsRaw));
    }, [optionsRaw]);

    // Destroy the context menu upon any click.
    useEffect(() => {
        document.addEventListener("click", () => { destroyContextMenu() });
    }, []);

    // Listen for higlass view config changes.
    useEffect(() => {
        hgRef.current.api.on('viewConfig', (newViewConfigString) => {
            const newViewConfig = JSON.parse(newViewConfigString);
            onViewConfig(newViewConfig);
            onViewConfigChangeProp(newViewConfigString);
        });         

        return () => hgRef.current.api.off('viewConfig');
    }, [hgRef, onViewConfigChangeProp]);

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
        <div className="chw-root">
            {hgComponent}
            {muiltivecTrackIds.map(({ viewId, trackId, trackTilesetId }, i) => (
                <TrackWrapper
                    key={i}
                    options={getTrackWrapperOptions(options, viewId, trackId)}
                    multivecTrack={getTrackObject(viewId, trackId)}
                    multivecTrackViewId={viewId}
                    multivecTrackTrackId={trackId}
                    multivecTrackTilesetId={trackTilesetId}
                    onAddTrack={(field, type, contains, position) => {
                        onAddTrack(viewId, trackId, field, type, contains, position);
                    }} 
                    onSortRows={(field, type, order) => {
                        onSortRows(viewId, trackId, field, type, order);
                    }}
                    onSearchRows={(field, type, contains) => {
                        onSearchRows(viewId, trackId, field, type, contains);
                    }}
                    onFilterRows={(field, type, condition) => {
                        onFilterRows(viewId, trackId, field, type, condition);
                    }}
                    onMetadataLoad={() => {
                        onMetadataLoad(viewId, trackId);
                    }}
                    drawRegister={drawRegister}
                />
            ))}
            {Array.from(new Set(muiltivecTrackIds.map(d => d.viewId))).map((viewId, i) => (
                <ViewWrapper
                    key={i}
                    viewBoundingBox={getViewBoundingBox(viewId)}
                    viewportTracks={viewportTrackIds[viewId] ? viewportTrackIds[viewId].map(d => getTrackObject(viewId, d.trackId)) : []}
                    multivecTrack={getTrackObject(viewId, muiltivecTrackIds.filter(d => d.viewId === viewId)[0].trackId)}
                    onSelectGenomicInterval={(startProp, endProp, uid) => {
                        const currViewConfig = hgRef.current.api.getViewConfig();
                        const newViewConfig = updateViewConfigOnSelectGenomicInterval(currViewConfig, viewId, startProp, endProp, uid);
                        hgRef.current.api.setViewConfig(newViewConfig).then(() => {
                            onViewConfig(newViewConfig);
                        });
                    }}
                    onViewportRemove={(viewportId) => {
                        const currViewConfig = hgRef.current.api.getViewConfig();
                        const newViewConfig = removeViewportFromViewConfig(currViewConfig, viewId, viewportId);
                        hgRef.current.api.setViewConfig(newViewConfig).then(() => {
                            onViewConfig(newViewConfig);
                        });
                    }}
                    drawRegister={drawRegister}
                />
            ))}
            <Tooltip />
            <ContextMenu/>
        </div>
    );
}
