import React, { useState, useEffect } from "react";
import PubSub from 'pubsub-js';
import { EVENT } from './utils/constants.js';
import { CLOSE } from './utils/icons.js';
import './TrackRowSearch.scss';

/**
 * Text field to serach for keywords. Subscribes to 'search' event via `PubSub`.
 * @example
 * <TrackRowSearch/>
 */
export default function TrackRowSearch(){

    const [left, setLeft] = useState(null);
    const [top, setTop] = useState(null);
    const [field, setField] = useState("");
    const [type, setType] = useState("");
    const [viewId, setViewId] = useState("");
    const [trackId, setTrackId] = useState("");

    // Styles
    const width = 180;
    const height = 30;
    const padding = 5;
    
    useEffect(() => {
        const searchToken = PubSub.subscribe(EVENT.SEARCH_OPEN, (msg, data) => {
            setLeft(data.left);
            setTop(data.top);
            setField(data.field);
            setType(data.type);
            setViewId(data.viewId);
            setTrackId(data.trackId);
        });

        return () => {
            PubSub.unsubscribe(searchToken);
        };
    });

    function onKeywordChange(e) {
        const keyword = e.target.value;
        PubSub.publish(EVENT.SEARCH_CHANGE, {
            contains: keyword, field, type, viewId, trackId
        });
    }

    function onSearchClose() {
        setLeft(null);
        setTop(null);
        setField("");
        setViewId("");
        setTrackId("");
    }

    return (
        <div
            className="cistrome-hgw-searchbox"
            style={{
                display: ((left !== null && top !== null) ? 'flex' : 'none'),
                left: left - (width + padding * 2) / 2,
                top: top - (height + padding * 2) - 30,
                padding
            }}
        >
            <input
                type="text"
                name="default name"
                placeholder="keyword"
                onChange={onKeywordChange}
                style={{ width, height }}
            />
            
            <svg className="chgw-button-sm chgw-search-close-button"
                onClick={onSearchClose} viewBox={CLOSE.viewBox}>
                <title>Close search box</title>
                <path d={CLOSE.path} fill="currentColor"/>
            </svg>
        </div>
    );
}