import React, { useRef, useState, useEffect } from 'react';
import PubSub from 'pubsub-js';

import { SEARCH, SORT_ASC, SORT_DESC, SORT_TREE } from './utils/icons.js';

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

    const buttons = [];
    if(type !== "tree") {
        buttons.push({
            onClick: onSortAscClick,
            icon: SORT_ASC,
            title: "Sort rows in ascending order"
        });
        buttons.push({
            onClick: onSortDescClick,
            icon: SORT_DESC,
            title: "Sort rows in descending order"
        });
    } else {
        buttons.push({
            onClick: onSortAscClick,
            icon: SORT_TREE,
            title: "Sort rows by hierarchy leaf order"
        })
    }

    if(onSearchRows) {
        buttons.push({
            onClick: onSearchClick,
            icon: SEARCH,
            title: "Search keywords"
        })
    }

    return (
        <div>
            <div ref={divRef}
                className={"chgw-control"}
                style={{
                    visibility: isVisible ? "visible" : "hidden"
                }}>
                {buttons.map((button, i) => {
                    let positionClass = "chgw-button-middle";
                    if(buttons.length > 1) {
                        if(i === 0) {
                            positionClass = "chgw-button-left"
                        } else if(i === buttons.length - 1) {
                            positionClass = "chgw-button-right"
                        }
                    }
                    return (
                        <svg 
                            key={button.title}
                            className={`chgw-button-sm ${positionClass}`}
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
                    top={searchTop}
                    left={searchLeft}
                    onChange={onSearchChange}
                    onClose={onSearchClose}
                />
            ) : null}
        </div>
    )
}