import React, { useRef, useEffect, useState, useMemo, useCallback, useContext } from 'react';
import isEqual from 'lodash/isEqual';
import clamp from 'lodash/clamp';

import { HiGlassComponent } from 'higlass';
import higlassRegister from 'higlass-register';
import StackedBarTrack from 'higlass-multivec/es/StackedBarTrack.js';

import { InfoContext, ACTION } from './utils/contexts.js';
import { selectRows, highlightRowsFromSearch } from './utils/select-rows.js';
import TrackWrapper from './TrackWrapper.js';
import ViewWrapper from './ViewWrapper.js';
import Tooltip from './Tooltip.js';
import ContextMenu, { destroyContextMenu } from './ContextMenu.js';
import CistromeToolkit from './CistromeToolkit.js';

import { 
    DEFAULT_OPTIONS_KEY,
    processWrapperOptions, 
    getTrackWrapperOptions,
    getWrapperSubOptions,
    updateWrapperOptions,
    addTrackWrapperOptions,
    getHighlightKeyByFieldType
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
import { removeItemFromArray, modifyItemInArray, insertItemToArray } from './utils/array.js';

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
    const [isWheelListening, setIsWheelListening] = useState(false);
    
    const context = useContext(InfoContext);

    // Set initial sorting, filtering, and highlighting.
    const onMetadataLoad = useCallback((viewId, trackId) => {
        if(!context.state[viewId] || !context.state[viewId][trackId]) {
            return;
        }
        const trackOptions = getTrackWrapperOptions(options, viewId, trackId);
        const rowInfo = context.state[viewId][trackId].rowInfo;

        // Aggregate, filter, and sort
        const newSelectedRows = selectRows(rowInfo, trackOptions);
        setTrackSelectedRows(viewId, trackId, newSelectedRows);
        
        // Highlight
        if(trackOptions.rowHighlight) {
            const { field, type } = trackOptions.rowHighlight;
            const condition = type === "nominal" ? 
                trackOptions.rowHighlight.contains :
                trackOptions.rowHighlight.range;
            const newHighlitRows = highlightRowsFromSearch(rowInfo, field, type, condition, trackOptions);
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
    const onHighlightRows = useCallback((viewId, trackId, field, type, condition) => {
        const highlightKey = getHighlightKeyByFieldType(type, condition);
        const newRowHighlight = { field, type, [highlightKey]: condition };
        const newOptions = updateWrapperOptions(options, newRowHighlight, "rowHighlight", viewId, trackId, { isReplace: true });
        const trackOptions = getTrackWrapperOptions(newOptions, viewId, trackId);

        // Notice: Highlighting options are specified only in the wrapper options.
        const newHighlitRows = highlightRowsFromSearch(context.state[viewId][trackId].rowInfo, field, type, condition, trackOptions);
        setHighlitRows(viewId, trackId, newHighlitRows);
        setOptions(newOptions);
    }, [options]);

    // Callback function for filtering.
    const onFilterRows = useCallback((viewId, trackId, field, type, condition, isRemove) => {
        let newSubOptions = getWrapperSubOptions(options, "rowFilter", viewId, trackId);
        if(!newSubOptions) newSubOptions = [];
        let fieldOption = newSubOptions.find(d => d.field === field);
        
        if(isRemove) {
            // Remove filters in a certain row track.
            if(!fieldOption) {
                // Nothing to remove.
            } else if(type === "nominal" || type == "link") {
                // Remove elements of `notOneOf` from the cetain field in the array of `condition`.
                const fieldOptionIndex = newSubOptions.indexOf(fieldOption);
                condition.forEach(one => {
                    const index = fieldOption.notOneOf.indexOf(one);
                    fieldOption.notOneOf.splice(index, 1);
                });
                newSubOptions[fieldOptionIndex] = fieldOption;
            } else if(type === "quantitative" || type === "tree") {
                // Simply remove `range` or `subtree` and `minSimilarity` filters from the cetain field.
                newSubOptions = removeItemFromArray(newSubOptions, newSubOptions.indexOf(fieldOption));
            }
        } else {
            // Add filters in a certain row track.
            if(type === "nominal" || type === "link") {
                // Add elements of `notOneOf` in not already presented.
                if(!fieldOption) {
                    fieldOption = { field, type, notOneOf: [] }
                }
                const fieldOptionIndex = newSubOptions.indexOf(fieldOption);
                condition.forEach(one => {
                    if(fieldOption.notOneOf.indexOf(one) === -1) {
                        fieldOption.notOneOf.push(one);
                    }
                });
                if(fieldOptionIndex === -1) {
                    newSubOptions = insertItemToArray(newSubOptions, 0, fieldOption);
                } else {
                    newSubOptions = modifyItemInArray(newSubOptions, fieldOptionIndex, fieldOption);
                }
            } else if(type === "quantitative") {
                // Replace with the incoming.
                if(fieldOption) {
                    const fieldOptionIndex = newSubOptions.indexOf(fieldOption);
                    newSubOptions = modifyItemInArray(newSubOptions, fieldOptionIndex, { field, type, range: condition });
                } else {
                    newSubOptions = insertItemToArray(newSubOptions, 0, { field, type, range: condition });
                }
            } else if(type === "tree") {
                // Replace with the incoming.
                const key = Array.isArray(condition) ? "subtree" : "minSimilarity";
                if(fieldOption) {
                    const fieldOptionIndex = newSubOptions.indexOf(fieldOption);
                    newSubOptions = modifyItemInArray(newSubOptions, fieldOptionIndex, {
                        ...fieldOption, [key]: condition // We do not want to remove a `subtree` or `minSimilarity` option if exists.
                    });
                } else {
                    newSubOptions = insertItemToArray(newSubOptions, 0, { field, type, [key]: condition });
                }
            }
        }
        
        let newOptions = updateWrapperOptions(options, newSubOptions, "rowFilter", viewId, trackId, { isReplace: true });
        newOptions = updateWrapperOptions(newOptions, undefined, "rowHighlight", viewId, trackId, { isReplace: true });
        newOptions = updateWrapperOptions(newOptions, undefined, "rowZoom", viewId, trackId, { isReplace: true });

        const trackOptions = getTrackWrapperOptions(newOptions, viewId, trackId);
        const newSelectedRows = selectRows(context.state[viewId][trackId].rowInfo, trackOptions);
        
        setHighlitRows(viewId, trackId, undefined);
        setTrackSelectedRows(viewId, trackId, newSelectedRows);
        setOptions(newOptions);
    }, [options]);

    // Callback function for vertical (row) zooming.
    const onZoomRows = useCallback((viewId, trackId, y, deltaY, deltaMode) => {
        const oldWrapperOptions = getTrackWrapperOptions(options, viewId, trackId);
        const trackNumRowsTotal = context.state[viewId][trackId].rowInfo.length;
        const trackNumRowsSelected = context.state[viewId][trackId].selectedRows.length;

        const oldRowZoom = oldWrapperOptions.rowZoom || { level: 1.0, top: 0.0, numRows: trackNumRowsSelected || trackNumRowsTotal };
        const numRows = oldRowZoom.numRows;
        const rowUnit = 1.0 / numRows;
        const oldLevel = oldRowZoom.level;
        const oldTop = oldRowZoom.top;
        const oldNumRows = oldLevel / rowUnit;

        const delta = deltaY;
        const factor = 1 + delta/12;
        const newLevel = clamp(oldRowZoom.level * factor, rowUnit, 1.0);
        const levelDiff = newLevel - oldLevel;
        const absoluteY = oldTop + ((y * oldNumRows) * rowUnit);
        const newTop = clamp(absoluteY - (levelDiff/2), 0.0, 1.0 - newLevel);

        const newRowZoom = { level: newLevel, top: newTop, numRows: numRows };
        const newWrapperOptions = updateWrapperOptions(options, newRowZoom, "rowZoom", viewId, trackId, { isReplace: true });
        const newTrackOptions = getTrackWrapperOptions(newWrapperOptions, viewId, trackId);
        const newSelectedRows = selectRows(context.state[viewId][trackId].rowInfo, newTrackOptions);

        setHighlitRows(viewId, trackId, undefined); // TODO: figure out how to update highlit rows y
        setTrackSelectedRows(viewId, trackId, newSelectedRows);
        setOptions(newWrapperOptions);
    }, [options]);

    // Callback function for adding a track.
    const onAddTrack = useCallback((viewId, trackId, field, type, notOneOf, position) => {
        if(viewId === DEFAULT_OPTIONS_KEY || trackId === DEFAULT_OPTIONS_KEY) {
            console.log("A view or track ID is a default ID.");
            return;
        }
        const newRowFilter = { field, type, notOneOf };
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
        const clickHandler = () => { destroyContextMenu() };
        document.addEventListener("click", clickHandler);
        return () => document.removeEventListener("click", clickHandler);
    }, []);

    // Listen for key events in order to start listening for wheel events.
    useEffect(() => {
        const keydownHandler = (keyEvent) => {
            if(keyEvent.keyCode === 89) {
                setIsWheelListening(true);
            }
        };
        const keyupHandler = (keyEvent) => {
            if(keyEvent.keyCode === 89) {
                setIsWheelListening(false);
            }
        };
        document.addEventListener("keydown", keydownHandler);
        document.addEventListener("keyup", keyupHandler);
        return () => {
            document.removeEventListener("keydown", keydownHandler);
            document.removeEventListener("keyup", keyupHandler);
        };
    }, [hgRef]);

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

    //console.log("CistromeHGWConsumer.render");
    return (
        <div className="chw-root">
            {hgComponent}
            {muiltivecTrackIds.map(({ viewId, trackId, trackTilesetId }, i) => (
                <TrackWrapper
                    key={i}
                    isWheelListening={isWheelListening}
                    options={getTrackWrapperOptions(options, viewId, trackId)}
                    multivecTrack={getTrackObject(viewId, trackId)}
                    multivecTrackViewId={viewId}
                    multivecTrackTrackId={trackId}
                    multivecTrackTilesetId={trackTilesetId}
                    onAddTrack={(field, type, notOneOf, position) => {
                        onAddTrack(viewId, trackId, field, type, notOneOf, position);
                    }}
                    onSortRows={(field, type, order) => {
                        onSortRows(viewId, trackId, field, type, order);
                    }}
                    onHighlightRows={(field, type, condition) => {
                        onHighlightRows(viewId, trackId, field, type, condition);
                    }}
                    onFilterRows={(field, type, condition, isRemove) => {
                        onFilterRows(viewId, trackId, field, type, condition, isRemove);
                    }}
                    onZoomRows={(y, deltaY, deltaMode) => {
                        onZoomRows(viewId, trackId, y, deltaY, deltaMode);
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
            <CistromeToolkit/>
        </div>
    );
}
