import React, { useRef, useState, useEffect } from 'react';
import PubSub from 'pubsub-js';

import { EVENT } from './utils/constants.js';
import { SEARCH, SORT_ASC, SORT_DESC } from './utils/icons.js';

import TrackRowSearch from './TrackRowSearch.js';

import './TrackRowInfoControl.scss';

const LOCAL_EVENT_SEARCH_OPEN = "search-open";

/**
 * Component with control buttons for each vertical track (for sorting, searching, etc).
 * @prop {boolean} isVisible The visibility of the control.
 * @prop {object} fieldInfo JSON object of the name and type of an attribute.
 * @prop {function} onSortRows The function to call upon a sort interaction.
 * @prop {function} onSearchRows The function to call upon a search interaction.
 */
export default function TrackRowInfoControl(props){
    const {
        isVisible, 
        fieldInfo,
        onSortRows,
        onSearchRows
    } = props;

    const divRef = useRef();
    const [isSearching, setIsSearching] = useState(false);
    const [searchTop, setSearchTop] = useState(null);
    const [searchLeft, setSearchLeft] = useState(null);

    const { field, type, title } = fieldInfo;
    console.assert(type !== "tree");
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
        onSearchRows(controlField, controlType, value);
    }

    return (
        <div>
            <div ref={divRef}
                className={"chgw-control"}
                style={{
                    visibility: isVisible ? "visible" : "hidden"
                }}>
                <svg className="chgw-button-sm chgw-button-top"
                    onClick={onSortAscClick} viewBox={SORT_ASC.viewBox}>
                    <title>Sort rows in ascending order</title>
                    <path d={SORT_ASC.path} fill="currentColor"/>
                </svg>
                <svg className="chgw-button-sm chgw-button-middle"
                    onClick={onSortDescClick} viewBox={SORT_DESC.viewBox}>
                    <title>Sort rows in descending order</title>
                    <path d={SORT_DESC.path} fill="currentColor"/>
                </svg>
                <svg className="chgw-button-sm chgw-button-bottom"
                    onClick={onSearchClick} viewBox={SEARCH.viewBox}>
                    <title>Search keywords</title>
                    <path d={SEARCH.path} fill="currentColor"/>
                </svg>
            </div>
            {isSearching ? (
                <TrackRowSearch 
                    top={searchTop}
                    left={searchLeft}
                    onChange={onSearchChange}
                    onClose={onSearchClose}
                />
            ) : null}
        </div>
    )
}