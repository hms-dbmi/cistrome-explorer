import React, { forwardRef, useRef, useEffect, useState, useMemo, useCallback, useContext } from 'react';
import isEqual from 'lodash/isEqual';
import clamp from 'lodash/clamp';

import { HiGlassComponent } from 'higlass';

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
    getWrapperSubOptions,
    updateWrapperOptions,
    addTrackWrapperOptions,
    getHighlightKeyByFieldType,
    getConditionFromHighlightOption
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
    removeViewportFromViewConfig,
    getTrackIdsFromViewConfig,
    removeTopTrackFromViewConfig
} from './utils/viewconf.js';
import { wrapSvg } from './utils/wrap-svg.js';

import './HiGlassMetaConsumer.scss';
import cloneDeep from 'lodash/cloneDeep';
import { removeItemFromArray, modifyItemInArray, insertItemToArray } from './utils/array.js';
import { CLOSE } from './utils/icons.js';
import { HG38_START_POSITIONS } from './utils/chromsizes.js';

const hgOptionsBase = {
    sizeMode: 'bounded', // Stretch the height of HiGlass to its container <div/>
    pixelPreciseMarginPadding: false,
    bounded: true,
    containerPaddingX: 0,
    containerPaddingY: 0,
    viewMarginTop: 0,
    viewMarginBottom: 0,
    viewMarginLeft: 0,
    viewMarginRight: 0,
    viewPaddingTop: 0,
    viewPaddingBottom: 0,
    viewPaddingLeft: 0,
    viewPaddingRight: 0,
};

// If this string tag is contained in the track id, ignore showing metadata visualization.
export const NO_METAVIS_TAG_TRACKID = '-no-metavis-';
export const REMOVE_ALLOWED_TAG_TRACKID = '-allow-remove-';

/**
 * HiGlassMeta passes its props through, and wraps this component with the context provider.
 * @prop {object} viewConfig A HiGlass viewConfig object.
 * @prop {(object|object[])} options Options for the wrapper component.
 * @prop {function} onViewChanged A function to call upon change of the view config and option. Optional.
 * @prop {function} onGenomicIntervalSearch A function to call upon searching for TFs by using the selected interval. Optional.
 */
const HiGlassMetaConsumer = forwardRef((props, ref) => {

    const {
        viewConfig: baseViewConfig,
        options: baseOptions,
        rowInfo: baseRowInfo,
        helpActivated,
        aggregateRowBy,
        onViewChanged: onViewChangedCallback,
        onGenomicIntervalSearch: onGenomicIntervalSearchCallback,
        onGeneSearch: onGeneSearchCallBack
    } = props;

    const hgRef = useRef();
    const drawRef = useRef({});

    const [options, setOptions] = useState(processWrapperOptions(baseOptions));
    const [multivecTrackIds, setMultivecTrackIds] = useState([]);
    const [viewportTrackIds, setViewportTrackIds] = useState({});
    const [stackedBarTrackIds, setStackedBarTrackIds] = useState([]);
    const [isWheelListening, setIsWheelListening] = useState(false);
    
    const context = useContext(InfoContext);

    // Initialize instances when we receive a new demo.
    useEffect(() => {
        resetTrackContext();
        setMultivecTrackIds([]);
        setViewportTrackIds({});
        setOptions(processWrapperOptions(baseOptions));
    }, [baseOptions, baseViewConfig, baseRowInfo]);
    
    useEffect(() => {
        // DEBUG
        // console.log('updated:', baseRowInfo);
    }, [baseRowInfo]);

    // Call a callback function when `options` changed.
    useEffect(() => {
        if(onViewChangedCallback) {
            onViewChangedCallback({ 
                options: JSON.parse(JSON.stringify(options))
                // ... add more here
            });
        }
    }, [options]);
    
    useEffect(() => {
        // Update Context based on the new options
        multivecTrackIds.forEach(({ viewId, trackId }) => {
            setMetadataToContext(viewId, trackId);
        });
    }, [options]);

    // This function stores `selectedRows` and `highlitRows` to the context
    // based on the updated options.
    const setMetadataToContext = useCallback((viewId, trackId) => {
        if(!context.state[viewId] || !context.state[viewId][trackId]) {
            // This means row information for this track is not yet loaded,
            // so we cannot store any additional information now.
            return;
        }
        const rowInfo = context.state[viewId][trackId].rowInfo;
        const trackOptions = getTrackWrapperOptions(options, viewId, trackId);

        // Aggregate, filter, and sort
        const newSelectedRows = selectRows(rowInfo, trackOptions); 
        if(!isEqual(newSelectedRows, context.state[viewId][trackId].selectedRows)) {
            // Update context only when there is any actual changes
            setSelectedRowsToViewConfig(viewId, trackId, newSelectedRows);
            setSelectedRows(viewId, trackId, newSelectedRows);
        }
        
        // Highlight
        let newHighlitRows = undefined; // `undefined` resets the hightlighting
        if(trackOptions.rowHighlight) {
            const { field, type } = trackOptions.rowHighlight;
            const condition = getConditionFromHighlightOption(trackOptions.rowHighlight);
            newHighlitRows = highlightRowsFromSearch(rowInfo, field, type, condition, trackOptions);
        }
        if(!isEqual(newHighlitRows, context.state[viewId][trackId].highlitRows)) {
            // Update context only when there is any actual changes
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
        const newTrackIds = getHMTrackIdsFromViewConfig(newViewConfig)
            .filter(({trackId}) => !trackId.includes(NO_METAVIS_TAG_TRACKID));
        
        // Add viewport projection horizontal track IDs for each view.
        const newViewportTrackIds = {};
        Array.from(new Set(newTrackIds.map(d => d.viewId))).forEach((viewId) => {    
            newViewportTrackIds[viewId] = getSiblingVPHTrackIdsFromViewConfig(newViewConfig, viewId);
        });

        // Tracks that are interactively added.
        const newStackedBarTrackIds = getTrackIdsFromViewConfig(newViewConfig, REMOVE_ALLOWED_TAG_TRACKID);

        // Get selected rows from the view config and update context.
        for(let trackIds of newTrackIds) {
            const newSelectedRows = getHMSelectedRowsFromViewConfig(newViewConfig, trackIds.viewId, trackIds.trackId);
            if(
                ! context.state[trackIds.viewId] 
                || !context.state[trackIds.viewId][trackIds.trackId] 
                || !isEqual(newSelectedRows, context.state[trackIds.viewId][trackIds.trackId].selectedRows)
            ) {
                // This must not force re-rendering React components since we might be in the 
                // middle of a rendering process by a `state` change.
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
        setStackedBarTrackIds(newStackedBarTrackIds);
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
        hgRef.current.api.setViewConfig(newViewConfig);
    }, [hgRef]);

    const removeTrack = useCallback((viewId, trackId) => {
        const currViewConfig = hgRef.current.api.getViewConfig();
        const newViewConfig = removeTopTrackFromViewConfig(currViewConfig, viewId, trackId);
        hgRef.current.api.setViewConfig(newViewConfig);
    }, [hgRef]);

    // HiGlassMeta APIs that can be called outside the library.
    useEffect(() => {
        ref.current = {
            api: {
                onOptions: (newOptions) => setOptions(processWrapperOptions(newOptions)),
                onRemoveAllFilters: () => removeAllFilters(),
                onRemoveAllSort: () => removeAllSort(),
                addNewTrack: (trackDef, viewId, position) => addNewTrack(trackDef, viewId, position),
                zoomTo: (keyword) => {
                    if(!keyword.includes('chr') || !keyword.includes('-')) {
                        console.warn('Genomic interval you entered is not in a correct format.');
                        return;
                    }
                    const chrStart = HG38_START_POSITIONS.find(d => d.chr === keyword.split(':')[0])?.position;
                    if(typeof chrStart === undefined) {
                        console.warn('Chromosome name is not valid', keyword.split(':')[0]);
                        return;
                    }
                    const [s, e] = keyword.split(':')[1].split('-');
                    const start = +s + chrStart;
                    const end = +e + chrStart;
                    Array.from(new Set(multivecTrackIds.map(d => d.viewId))).map((viewId) => {
                        hgRef.current.api.zoomTo(viewId, start, end, start, end, 1000);
                    });
                },
                zoomToGene: (gene) => {
                    Array.from(new Set(multivecTrackIds.map(d => d.viewId))).map((viewId) => {
                        hgRef.current.api.zoomToGene(viewId, gene, 1000);
                    });
                },
                suggestGene: (keyword, callback) => {
                    if(multivecTrackIds.length > 0) {
                        const viewId = multivecTrackIds[0].viewId;
                        hgRef.current.api.suggestGene(viewId, keyword, (suggestions) => {
                            if(suggestions && suggestions.length > 0) {
                                callback(suggestions);
                            }
                        });
                    }
                }
            }
        }
    }, [ref, addNewTrack, multivecTrackIds]);

    const setSelectedRowsToViewConfig = useCallback((viewId, trackId, selectedRows) => {
        const currViewConfig = hgRef.current.api.getViewConfig();
        const newViewConfig = updateViewConfigOnSelectRowsByTrack(currViewConfig, selectedRows, viewId, trackId);
        hgRef.current.api.setViewConfig(newViewConfig);
    }, [hgRef]);
    
    const resetTrackContext = useCallback(() => {
        multivecTrackIds.forEach(({ viewId, trackId }) => {
            context.dispatch({
                type: ACTION.RESET,
                viewId,
                trackId
            });
        });
    }, [multivecTrackIds]);

    const setSelectedRows = useCallback((viewId, trackId, selectedRows) => {
        context.dispatch({
            type: ACTION.SELECT_ROWS_RERENDER,
            viewId,
            trackId,
            selectedRows
        });
    });

    const setHighlitRows = useCallback((viewId, trackId, highlitRows) => {
        context.dispatch({
            type: ACTION.HIGHLIGHT_ROWS_RERENDER,
            viewId,
            trackId,
            highlitRows
        });
    });

    // Function for child components to call to "register" their draw functions.
    const drawRegister = useCallback((key, draw, options) => {
        drawRef.current[key] = { draw, options };
    }, [drawRef]);

    // Clear the drawRegister object when the options prop changes.
    useEffect(() => {
        drawRef.current = {};
    }, [drawRef, baseOptions]);

    // Listen for the `createSVG` event.
    useEffect(() => {
        hgRef.current.api.on('createSVG', (svg) => {
            return wrapSvg(svg, drawRef.current);
        });
        return () => hgRef.current.api.off('createSVG');
    }, [hgRef, drawRef]);
    
    const removeAllFilters = useCallback(() => {
        let newOptions = options;
        multivecTrackIds.forEach(({ viewId, trackId }) => {
            newOptions = updateWrapperOptions(newOptions, [], "rowFilter", viewId, trackId, { isReplace: true });
        });
        setOptions(newOptions);
    }, [multivecTrackIds]);

    const removeAllSort = useCallback(() => {
        let newOptions = options;
        multivecTrackIds.forEach(({ viewId, trackId }) => {
            newOptions = updateWrapperOptions(newOptions, [], "rowSort", viewId, trackId, { isReplace: true });
        });
        setOptions(newOptions);
    }, [multivecTrackIds]);

    // Callback function for sorting.
    const onSortRows = useCallback((viewId, trackId, field, type, order, isTrackIndependent) => {
        let newOptions;
        if(isTrackIndependent) {
            // We need to modify the track-level `sort` option
            const { rowInfoAttributes } = getTrackWrapperOptions(options, viewId, trackId);
            const newRowInfoAttributes = rowInfoAttributes.map(fieldInfo => {
                if(fieldInfo.field === field && fieldInfo.type === type) {
                    return { ...fieldInfo, sort: order }
                } else {
                    return fieldInfo;
                }
            });
            newOptions = updateWrapperOptions(options, newRowInfoAttributes, "rowInfoAttributes", viewId, trackId, { isReplace: true });
        }
        else {
            const newRowSort = [ { field, type, order } ];
            newOptions = updateWrapperOptions(options, newRowSort, "rowSort", viewId, trackId, { isReplace: true });
        }
        setOptions(newOptions);
    }, [options]);

    // Callback function for searching and highlighting.
    const onHighlightRows = useCallback((viewId, trackId, field, type, condition) => {
        const highlightKey = getHighlightKeyByFieldType(type, condition);
        const newRowHighlight = condition ? { field, type, [highlightKey]: condition } : {};
        const newOptions = updateWrapperOptions(options, newRowHighlight, "rowHighlight", viewId, trackId, { isReplace: true });
        setOptions(newOptions);
    }, [options]);

    // Aggregating rows
    useEffect(() => {
        if(!multivecTrackIds || multivecTrackIds.length === 0) {
            // do not have enough information
            return;
        }
        const { viewId, trackId } = multivecTrackIds[0];
        const newRowAggregate = aggregateRowBy ? [{ field: aggregateRowBy, type: 'nominal', notOneOf: [] }] : [];
        const newOptions = updateWrapperOptions(options, newRowAggregate, "rowAggregate", viewId, trackId, { isReplace: true });
        setOptions(newOptions);
    }, [aggregateRowBy]);

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
                // Simply remove `range` or `ancestors` and `minSimilarity` filters from the cetain field.
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
                const key = Array.isArray(condition) ? "ancestors" : "minSimilarity";
                if(fieldOption) {
                    const fieldOptionIndex = newSubOptions.indexOf(fieldOption);
                    newSubOptions = modifyItemInArray(newSubOptions, fieldOptionIndex, {
                        ...fieldOption, [key]: condition // We do not want to remove a `ancestors` or `minSimilarity` option if exists.
                    });
                } else {
                    newSubOptions = insertItemToArray(newSubOptions, 0, { field, type, [key]: condition });
                }
            }
        }
        
        let newOptions = updateWrapperOptions(options, newSubOptions, "rowFilter", viewId, trackId, { isReplace: true });
        newOptions = updateWrapperOptions(newOptions, {}, "rowHighlight", viewId, trackId, { isReplace: true });
        newOptions = updateWrapperOptions(newOptions, {}, "rowZoom", viewId, trackId, { isReplace: true });
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
        setSelectedRowsToViewConfig(viewId, trackId, newSelectedRows);
        setOptions(newWrapperOptions);
    }, [options]);
    
    // Callback function for adding a track.
    const onAddTrack = useCallback((viewId, trackId, field, type, notOneOf, position, selected = '') => {
        if(viewId === DEFAULT_OPTIONS_KEY || trackId === DEFAULT_OPTIONS_KEY) {
            console.log("A view or track ID is a default ID.");
            return;
        }
        const newRowFilter = { field, type, notOneOf };
        const currViewConfig = hgRef.current.api.getViewConfig();
        const newTrackId = getUniqueViewOrTrackId(currViewConfig, { 
            baseId: trackId, 
            idKey: "trackId", 
            interfix: REMOVE_ALLOWED_TAG_TRACKID
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
            height: 30,
            uid: newTrackId,
            type: 'horizontal-stacked-bar',
            options: {
                ...newTrackDef.options,
                name: selected,
                selectRows: newSelectedRows,
                trackBorderWidth: 0,
                trackBorderColor: "white",
                // barBorder: false,
                labelColor: "black",
                backgroundColor: "#F6F6F6",
                colorScale: ["lightgray"],
            }
        }
        addNewTrack(newTrackDef, viewId, position);
        setOptions(newOptions);
    }, [options]);

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
        });

        return () => hgRef.current.api.off('viewConfig');
    }, [hgRef]);

    useEffect(() => {
        hgRef.current.api.on('geneSearch', (e) => {
            onGeneSearchCallBack(e.geneSymbol);
        });

        return () => hgRef.current.api.off('geneSearch');
    }, [hgRef]);

    // We only want to render HiGlass once.
    const hgComponent = useMemo(() => {
        const hgOptions = {
            ...hgOptionsBase,
            onViewConfLoaded: () => {
                // This is called only once and is not duplicately called 
                // with `on('viewConfig')` API.
                onViewConfig(baseViewConfig);
            }
        };

        console.log("HiGlassComponent.render");
        return (
            <HiGlassComponent
                viewConfig={baseViewConfig}
                options={hgOptions}
                zoomFixed={false}
                ref={hgRef}
            />
        );
    }, [baseViewConfig, baseRowInfo]);

    //console.log("HiGlassWithMetadataConsumer.render");
    return (
        <div className="hm-root">
            {hgComponent}
            {multivecTrackIds.map(({ viewId, trackId, trackTilesetId }, i) => (
                <TrackWrapper
                    key={i}
                    isWheelListening={isWheelListening}
                    options={getTrackWrapperOptions(options, viewId, trackId)}
                    baseRowInfo={baseRowInfo}
                    multivecTrack={getTrackObject(viewId, trackId)}
                    multivecTrackViewId={viewId}
                    multivecTrackTrackId={trackId}
                    multivecTrackTilesetId={trackTilesetId}
                    onAddTrack={(field, type, notOneOf, position, selected) => {
                        onAddTrack(viewId, trackId, field, type, notOneOf, position, selected);
                    }}
                    onSortRows={(field, type, order, isTrackIndependent) => {
                        onSortRows(viewId, trackId, field, type, order, isTrackIndependent);
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
                    onMetadataInit={() => {
                        setMetadataToContext(viewId, trackId);
                    }}
                    helpActivated={helpActivated}
                    rowAggregated={aggregateRowBy !== undefined}
                    drawRegister={drawRegister}
                />
            ))}
            {stackedBarTrackIds.map(({ viewId, trackId }, i) => {
                const track = getTrackObject(viewId, trackId);
                if(!track) {
                    // did not find any track with the given viewid and trackId
                    return;  
                }

                const [left, top] = track.position;
                const [w, h] = track.dimensions;
                return (
                    <div style={{ left: left + w, top, position: 'absolute', height: h }}
                        onClick={() => removeTrack(viewId, trackId)}
                    >
                        <svg
                            className={'hm-button'}
                            style={{ color: "rgb(171, 171, 171)", background: "none", height: '30px' }}
                            viewBox={CLOSE.viewBox}
                        >
                            <title>Remove Track</title>
                            <path d={CLOSE.path} fill="currentColor"/>
                        </svg>
                    </div>
                );
            })}
            {Array.from(new Set(multivecTrackIds.map(d => d.viewId))).map((viewId, i) => (
                <ViewWrapper
                    key={i}
                    viewBoundingBox={getViewBoundingBox(viewId)}
                    viewportTracks={viewportTrackIds[viewId] ? viewportTrackIds[viewId].map(d => getTrackObject(viewId, d.trackId)) : []}
                    multivecTrack={getTrackObject(viewId, multivecTrackIds.filter(d => d.viewId === viewId)[0].trackId)}
                    onSelectGenomicInterval={(startProp, endProp, uid) => {
                        const currViewConfig = hgRef.current.api.getViewConfig();
                        const newViewConfig = updateViewConfigOnSelectGenomicInterval(currViewConfig, viewId, startProp, endProp, uid);
                        hgRef.current.api.setViewConfig(newViewConfig);
                    }}
                    onViewportRemove={(viewportId) => {
                        const currViewConfig = hgRef.current.api.getViewConfig();
                        const newViewConfig = removeViewportFromViewConfig(currViewConfig, viewId, viewportId);
                        hgRef.current.api.setViewConfig(newViewConfig);
                    }}
                    onGenomicIntervalSearch={onGenomicIntervalSearchCallback}
                    helpActivated={helpActivated}
                    drawRegister={drawRegister}
                />
            ))}
            <Tooltip/>
            <ContextMenu/>
        </div>
    );
});

export default HiGlassMetaConsumer;
