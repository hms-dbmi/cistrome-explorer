import React, { useRef, useState, useEffect } from 'react';
import PubSub from 'pubsub-js';

import { SORT_ASC, SORT_DESC, FILTER, RESET, TOGGLE_ON } from './utils/icons.js';
import TrackRowFilter from './TrackRowFilter.js';
import { getAggregatedValue } from './utils/aggregate.js';

const LOCAL_EVENT_FILTER_OPEN = "filter-open";

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
 */
export default function TrackRowInfoControl(props){
    const {
        isLeft,
        isVisible, 
        field, type, title, aggFunction,
        sortAsceButtonHighlit,
        sortDescButtonHighlit,
        filterButtonHighlit,
        onSortRows,
        onHighlightRows,
        onFilterRows,
        toggleMinSimBar,
        rowInfo,
        filterInfo
    } = props;

    const divRef = useRef();
    const [isFiltering, setIsFiltering] = useState(false);
    const [filterTop, setFilterTop] = useState(null);
    const [FilterLeft, setFilterLeft] = useState(null);

    const controlField = (type === "url" && title ? title : field);
    const controlType = (type === "url" ? "nominal" : type);

    // Subscribe to the filter open events of other TrackRowInfoControl components,
    // so that only one filter is open at a time.
    useEffect(() => {
        const filterOpenToken = PubSub.subscribe(LOCAL_EVENT_FILTER_OPEN, (msg, otherDivRef) => {
            if(divRef !== otherDivRef) {
                setIsFiltering(false);
            }
        });

        return () => PubSub.unsubscribe(filterOpenToken);
    })

    function onSortAscClick() {
        onSortRows(controlField, controlType, "ascending");
    }
    function onSortDescClick() {
        onSortRows(controlField, controlType, "descending");
    }
    function onFilterClick(event) {
        const parentRect = divRef.current.getBoundingClientRect();
        
        setIsFiltering(true);
        setFilterTop(event.clientY - parentRect.y);
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
    if(type !== "tree") {
        buttons.push({
            onClick: onSortAscClick,
            icon: SORT_ASC,
            title: "Sort rows in ascending order",
            highlit: sortAsceButtonHighlit
        });
        buttons.push({
            onClick: onSortDescClick,
            icon: SORT_DESC,
            title: "Sort rows in descending order",
            highlit: sortDescButtonHighlit
        });
    }

    if(onHighlightRows && !Array.isArray(field)) {
        buttons.push({
            onClick: onFilterClick,
            icon: FILTER,
            title: "Filter rows",
            highlit: filterButtonHighlit
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
            highlit: false
        });
    }

    return (
        <div>
            <div ref={divRef}
                className={"hm-button-sm-container-vertical"}
                style={{
                    top: "4px",
                    left: "4px",
                    visibility: isVisible ? "visible" : "hidden"
                }}>
                {buttons.map((button, i) => {
                    let positionClass = "hm-button-middle";
                    if(buttons.length > 1) {
                        if(i === 0) {
                            positionClass = "hm-button-top"
                        } else if(i === buttons.length - 1) {
                            positionClass = "hm-button-bottom"
                        }
                    }
                    return (
                        <svg 
                            key={button.title}
                            className={`${button.highlit ? "hm-button-sm-hl" : "hm-button-sm"} ${positionClass}`}
                            onClick={button.onClick} 
                            viewBox={button.icon.viewBox}
                        >
                            <title>{button.title}</title>
                            <path d={button.icon.path} fill="currentColor"/>
                        </svg>
                    );
                })}
            </div>
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
    )
}
