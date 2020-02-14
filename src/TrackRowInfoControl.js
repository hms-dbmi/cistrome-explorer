import React from 'react';
import PubSub from 'pubsub-js';

import { EVENT } from './utils/constants.js';
import { SEARCH, SORT_ASC, SORT_DESC } from './utils/icons.js';
import './TrackRowInfoControl.scss';

/**
 * Component with control buttons for each vertical track (for sorting, searching, etc).
 * @prop {boolean} isVisible The visibility of the control.
 * @prop {object} fieldInfo JSON object of the name and type of an attribute.
 */
export default function TrackRowInfoControl(props){
    const {
        viewId, trackId,
        isVisible, 
        fieldInfo
    } = props;

    function onSortAscClick() {
        const { field, type, title } = fieldInfo;
        if(type === "tree") return;
        PubSub.publish(EVENT.SORT, {
            field: (type === "url" && title ? title : field),
            type,
            order: "ascending",
            viewId,
            trackId,
        });
    }
    function onSortDescClick() {
        const { field, type, title } = fieldInfo;
        if(type === "tree") return;
        PubSub.publish(EVENT.SORT, {
            field: (type === "url" && title ? title : field),
            type,
            order: "descending",
            viewId,
            trackId,
        });
    }

    return (
        <div 
            className={"chgw-control"}
            style={{
                visibility: isVisible ? "visible" : "hidden"
            }}>
            <svg className="chgw-button-sm chgw-button-left"
                onClick={onSortAscClick} viewBox={SORT_ASC.viewBox}>
                <title>Sort rows in ascending order</title>
                <path d={SORT_ASC.path} fill="currentColor"/>
            </svg>
            <svg className="chgw-button-sm chgw-button-middle"
                onClick={onSortDescClick} viewBox={SORT_DESC.viewBox}>
                <title>Sort rows in descending order</title>
                <path d={SORT_DESC.path} fill="currentColor"/>
            </svg>
            <svg className="chgw-button-sm chgw-button-right"
                viewBox={SEARCH.viewBox}>
                <title>Search keywords</title>
                <path d={SEARCH.path} fill="currentColor"/>
            </svg>
        </div>
    )
}