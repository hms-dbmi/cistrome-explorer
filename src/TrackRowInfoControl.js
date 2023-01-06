import React, { useRef, useState, useEffect } from "react";
import PubSub from "pubsub-js";

import { SORT_ASC, SORT_DESC, FILTER, RESET, TOGGLE_ON, PLUS, ARROW_MOVE } from "./utils/icons.js";
import TrackRowFilter from "./TrackRowFilter.js";
import { getAggregatedValue } from "./utils/aggregate.js";
import { destroyTooltip, publishHelpTooltip } from "./Tooltip.js";
import TrackAddNewTrack from "./TrackAddNewTrack.js";

const LOCAL_EVENT_FILTER_OPEN = "filter-open";
const LOCAL_EVENT_ADD_TRACK_OPEN = LOCAL_EVENT_FILTER_OPEN;

/**
 * Component with control buttons for each vertical track (for sorting, filtering, etc).
 * @prop {boolean} isLeft Is this view on the left side of the HiGlass track?
 * @prop {boolean} isVisible The visibility of the control.
 * @prop {object} fieldInfo JSON object of the name and type of an attribute.
 * @prop {boolean} sortAsceButtonHighlit Is sorting in ascending order applied in this vertical track?
 * @prop {boolean} sortDescButtonHighlit Is sorting in descending order applied in this vertical track?
 * @prop {boolean} filterButtonHighlit Is filter applied in this vertical track?
 * @prop {function} onSortRows The function to call upon a sort interaction.
 * @prop {function} onHighlightRows The function to call upon a search interaction.
 * @prop {function} onFilterRows The function to call upon a filter interaction.
 * @prop {function} toggleMinSimBar Toggle showing the minimum similarity bar in dendrogram.
 * @prop {object[]} rowInfo The array of JSON Object containing row information.
 * @prop {object} filterInfo The options for filtering rows of the field used in this track.
 * @prop {boolean} helpActivated Whether to show help instructions or not.
 */
export default function TrackRowInfoControl(props){
    const {
        isLeft,
        isVisible, 
        top,
        field, type, alt, title, aggFunction, resolveYScale,
        sortAsceButtonHighlit,
        sortDescButtonHighlit,
        filterButtonHighlit,
        onSortRows,
        onHighlightRows,
        onFilterRows,
        toggleMinSimBar,
        rowInfo,
        filterInfo,
        helpActivated
    } = props;

    const divRef = useRef();
    const [isFiltering, setIsFiltering] = useState(false);
    const [isAddingTrack, setIsAddingTrack] = useState(false);
    const [filterTop, setFilterTop] = useState(null);
    const [FilterLeft, setFilterLeft] = useState(null);

    const controlField = (type === "url" && alt ? alt : field);
    const controlType = (type === "url" ? "nominal" : type);

    // Subscribe to the filter open events of other TrackRowInfoControl components,
    // so that only one filter is open at a time.
    useEffect(() => {
        const filterOpenToken = PubSub.subscribe(LOCAL_EVENT_FILTER_OPEN, (msg, otherDivRef) => {
            if(divRef !== otherDivRef) {
                setIsFiltering(false);
                setIsAddingTrack(false);
            }
        });

        return () => PubSub.unsubscribe(filterOpenToken);
    });

    function onSortAscClick() {
        onSortRows(controlField, controlType, "ascending", resolveYScale);
    }
    function onSortDescClick() {
        onSortRows(controlField, controlType, "descending", resolveYScale);
    }
    function onAddTrackClick(event) {
        const parentRect = divRef.current.getBoundingClientRect();
        
        setIsAddingTrack(true);
        setFilterTop(event.clientY);
        setFilterLeft(event.clientX - parentRect.x);

        PubSub.publish(LOCAL_EVENT_ADD_TRACK_OPEN, divRef);
    }
    function onAddingTrackClose () {
        setIsAddingTrack(false);
        onFilterClose();
    }
    function onFilterClick(event) {
        const parentRect = divRef.current.getBoundingClientRect();
        
        setIsFiltering(true);
        onAddingTrackClose();
        setFilterTop(event.clientY);
        setFilterLeft(event.clientX - parentRect.x);

        PubSub.publish(LOCAL_EVENT_FILTER_OPEN, divRef);
    }
    function onFilterClose() {
        setIsFiltering(false);
    }
    function onFilterChange(value) {
        onHighlightRows(controlField, controlType, value);
    }
    function onToggleMinSimBar() {
        toggleMinSimBar();
    }
    function onReset() {
        if(type === "nominal" || type == "link") {
            onFilterRows(field, type, rowInfo.map(
                d => getAggregatedValue(d, field, "nominal", aggFunction).toString()
            ), true);
        } else if(type === "quantitative" || type == "tree") {
            onFilterRows(field, type, [], true);
        }
    }

    const buttons = [];

    buttons.push({
        onClick: () => {}, // TODO:
        icon: ARROW_MOVE,
        title: "Move track",
        helpTitle: "Move track",
        helpSubtitle: "Move this track to another location."
    });
    
    buttons.push({
        onClick: onAddTrackClick,
        icon: PLUS,
        title: "Add new metadata track",
        helpTitle: "Add new metadata track",
        helpSubtitle: "Add a new track that visualize a new metadata."
    });

    if(type !== "tree") {
        buttons.push({
            onClick: onSortAscClick,
            icon: SORT_ASC,
            title: "Sort rows in ascending order",
            highlit: sortAsceButtonHighlit,
            helpTitle: "Sort Samples in Ascending Order",
            helpSubtitle: "Sort samples in ascending order based on the value shown in this track."
        });
        buttons.push({
            onClick: onSortDescClick,
            icon: SORT_DESC,
            title: "Sort rows in descending order",
            highlit: sortDescButtonHighlit,
            helpTitle: "Sort Samples in Descending Order",
            helpSubtitle: "Sort samples in descending order based on the value shown in this track."
        });
    }

    if(onHighlightRows && !Array.isArray(field)) {
        buttons.push({
            onClick: onFilterClick,
            icon: FILTER,
            title: "Filter rows",
            highlit: filterButtonHighlit,
            helpTitle: "Apply Filter",
            helpSubtitle: "Apply filters in this track to remove certain samples that are not interested."
        });
    }

    if(toggleMinSimBar) {
        buttons.push({
            onClick: onToggleMinSimBar,
            icon: TOGGLE_ON,
            title: "Show minimum similarity bar",
            highlit: filterButtonHighlit
        });
    }

    if(onFilterRows) {
        buttons.push({
            onClick: onReset,
            icon: RESET,
            title: "Remove all filters",
            highlit: false,
            helpTitle: "Remove All Filter",
            helpSubtitle: "Remove all filters that are applied in the tracks."
        });
    }

    return (
        <div>
            <div ref={divRef}
                className={"hm-button-sm-container-vertical " + (helpActivated ? "help-highlight" : "")}
                style={{
                    top: `${top + 4}px`,
                    left: "4px",
                    visibility: (isVisible || helpActivated) ? "visible" : "hidden"
                }}>
                {buttons.map((button, i) => {
                    let positionClass = "hm-button-middle";
                    if(buttons.length > 1) {
                        if(i === 0) {
                            positionClass = "hm-button-top";
                        } else if(i === buttons.length - 1) {
                            positionClass = "hm-button-bottom";
                        }
                    }
                    return (
                        <svg 
                            key={button.title}
                            className={`${button.highlit ? "hm-button-sm-hl" : "hm-button-sm"} ${positionClass}`}
                            onClick={button.onClick} 
                            viewBox={button.icon.viewBox}
                            onMouseMove={(e) => publishHelpTooltip(e,
                                button.helpTitle,
                                button.helpSubtitle,
                                helpActivated
                            )}
                            onMouseLeave={() => destroyTooltip()}
                        >
                            <title>{button.title}</title>
                            <path d={button.icon.path} fill="currentColor"/>
                        </svg>
                    );
                })}
            </div>
            {isAddingTrack ? (
                <TrackAddNewTrack // TODO: add new dialog
                    top={filterTop}
                    left={FilterLeft}
                    field={controlField}
                    type={controlType}
                    onChange={onFilterChange}
                    onFilterRows={onFilterRows}
                    onClose={onAddingTrackClose}
                    rowInfo={rowInfo}
                />
            ) : null}
            {isFiltering ? (
                <TrackRowFilter
                    isLeft={isLeft}
                    top={filterTop}
                    left={FilterLeft}
                    field={controlField}
                    type={controlType}
                    aggFunction={aggFunction}
                    onChange={onFilterChange}
                    onFilterRows={onFilterRows}
                    onClose={onFilterClose}
                    rowInfo={rowInfo}
                    filterInfo={filterInfo}
                />
            ) : null}
        </div>
    );
}
