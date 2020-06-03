import React, { useRef, useState, useEffect } from 'react';
import PubSub from 'pubsub-js';

import { SORT_ASC, SORT_DESC, FILTER, RESET, TOGGLE_ON } from './utils/icons.js';
import TrackRowSearch from './TrackRowSearch.js';

const LOCAL_EVENT_SEARCH_OPEN = "search-open";

/**
 * Component with control buttons for each vertical track (for sorting, searching, etc).
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
 * @prop {object[]} rowInfo Array of JSON objects, one object for each sample, without filtering/sorting based on selected rows.
 * @prop {object} filterInfo The options for filtering rows of the field used in this track.
 */
export default function TrackRowInfoControl(props){
    const {
        isLeft,
        isVisible, 
        fieldInfo,
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
    const [isSearching, setIsSearching] = useState(false);
    const [searchTop, setSearchTop] = useState(null);
    const [searchLeft, setSearchLeft] = useState(null);

    const { field, type, title, aggFunction } = fieldInfo;
    const controlField = (type === "url" && title ? title : field);
    const controlType = (type === "url" ? "nominal" : type);

    // Subscribe to the search open events of other TrackRowInfoControl components,
    // so that only one search is open at a time.
    useEffect(() => {
        const searchOpenToken = PubSub.subscribe(LOCAL_EVENT_SEARCH_OPEN, (msg, otherDivRef) => {
            if(divRef !== otherDivRef) {
                setIsSearching(false);
            }
        });

        return () => PubSub.unsubscribe(searchOpenToken);
    })

    function onSortAscClick() {
        onSortRows(controlField, controlType, "ascending");
    }
    function onSortDescClick() {
        onSortRows(controlField, controlType, "descending");
    }
    function onSearchClick(event) {
        const parentRect = divRef.current.getBoundingClientRect();
        
        setIsSearching(true);
        setSearchTop(event.clientY - parentRect.y);
        setSearchLeft(event.clientX - parentRect.x);

        PubSub.publish(LOCAL_EVENT_SEARCH_OPEN, divRef);
    }
    function onSearchClose() {
        setIsSearching(false);
    }
    function onSearchChange(value) {
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
            onClick: onSearchClick,
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
                className={"chw-button-sm-container-vertical"}
                style={{
                    top: "4px",
                    left: "4px",
                    visibility: isVisible ? "visible" : "hidden"
                }}>
                {buttons.map((button, i) => {
                    let positionClass = "chw-button-middle";
                    if(buttons.length > 1) {
                        if(i === 0) {
                            positionClass = "chw-button-top"
                        } else if(i === buttons.length - 1) {
                            positionClass = "chw-button-bottom"
                        }
                    }
                    return (
                        <svg 
                            key={button.title}
                            className={`${button.highlit ? "chw-button-sm-hl" : "chw-button-sm"} ${positionClass}`}
                            onClick={button.onClick} 
                            viewBox={button.icon.viewBox}
                        >
                            <title>{button.title}</title>
                            <path d={button.icon.path} fill="currentColor"/>
                        </svg>
                    );
                })}
            </div>
            {isSearching ? (
                <TrackRowSearch
                    isLeft={isLeft}
                    top={searchTop}
                    left={searchLeft}
                    field={controlField}
                    type={controlType}
                    aggFunction={aggFunction}
                    onChange={onSearchChange}
                    onFilterRows={onFilterRows}
                    onClose={onSearchClose}
                    rowInfo={rowInfo}
                    filterInfo={filterInfo}
                />
            ) : null}
        </div>
    )
}
