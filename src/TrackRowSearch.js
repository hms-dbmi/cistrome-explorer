import React, { useState, useEffect, useRef } from "react";
import PubSub from 'pubsub-js';
import { EVENT } from './utils/constants.js';
import { CLOSE } from './utils/icons.js';
import './TrackRowSearch.scss';

/**
 * Text field to serach for keywords.
 * @prop {number} top The top coordinate.
 * @prop {number} left The left coordinate.
 * @prop {function} onChange The function to call when the search keyword has changed.
 * @prop {function} onClose The function to call when the search field should be closed.
 * @example
 * <TrackRowSearch/>
 */
export default function TrackRowSearch(props) {

    const {
        top, left,
        onChange,
        onClose
    } = props;

    const inputRef = useRef();

    // Styles
    const width = 180;
    const height = 30;
    const padding = 5;
    
    useEffect(() => {
        inputRef.current.focus();
    });

    function onKeywordChange(e) {
        const keyword = e.target.value;
        onChange(keyword);
    }

    function onSearchClose() {
        onChange("");
        onClose();
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
                ref={inputRef}
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