import React from 'react';
import PubSub from 'pubsub-js';

import { EVENT } from './constants.js';
import { SEARCH, SORT_ASC, SORT_DESC } from './utils/icons.js';
import './TrackControl.scss';

/**
 * Component to control each vertical track.
 * @prop {number} left The control horizontal offset.
 * @prop {number} top The control vertical offset.
 * @prop {boolean} isVisible The visibility of the control.
 * @prop {object} fieldInfo JSON object of the name and type of an attribute.
 */
export default function TrackControl(props){
    const {
        left, top, isVisible, fieldInfo
    } = props;

    function onSortAscClick() {
        const { field, type, title } = fieldInfo;
        if(type === "tree") return;
        PubSub.publish(EVENT.SORT, {
            field: (type === "url" && title ? title : field),
            type,
            order: "ascending"
        });
    }
    function onSortDescClick() {
        const { field, type, title } = fieldInfo;
        if(type === "tree") return;
        PubSub.publish(EVENT.SORT, {
            field: (type === "url" && title ? title : field),
            type,
            order: "descending"
        });
    }

    return (
        <div 
            className={"chgw-control"}
            style={{
                top: `${top}px`,
                left: `${left}px`, 
                width: `${60}px`,
                height: `${20}px`,
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