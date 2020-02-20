import React, { useEffect, useRef, useState } from "react";
import { CLOSE, FILTER } from './utils/icons.js';
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
        field, type,
        onChange,
        onFilter,
        onClose
    } = props;

    const inputRef = useRef();
    const [keyword, setKeyword] = useState("");

    // Styles
    const width = 180;
    const height = 30;
    const padding = 5;
    
    useEffect(() => {
        inputRef.current.focus();
    });

    function onKeywordChange(e) {
        const newKeyword = e.target.value;
        setKeyword(newKeyword);
        onChange(newKeyword);
    }

    function onFilterClick() {
        const contains = keyword.toString();
        onFilter(field, type, contains);
        setKeyword("");
    }

    function onSearchClose() {
        onChange("");
        onClose();
        setKeyword("");
    }

    return (
        <div
            className="cistrome-hgw-searchbox"
            style={{
                display: ((left !== null && top !== null) ? 'flex' : 'none'),
                left: left - (width + padding * 2) / 2,
                top: top - (height + padding * 2) - 60,
                padding
            }}
        >
            <input
                ref={inputRef}
                className="cistrome-hgw-searchinput"
                type="text"
                name="default name"
                placeholder="keyword"
                onChange={onKeywordChange}
                style={{ width, height }}
            />
            
            <svg className="chgw-button-sm chgw-search-button"
                onClick={onFilterClick} viewBox={FILTER.viewBox}>
                <title>Filter rows using the keyword.</title>
                <path d={FILTER.path} fill="currentColor"/>
            </svg>
            <svg className="chgw-button-sm chgw-search-button"
                onClick={onSearchClose} viewBox={CLOSE.viewBox}>
                <title>Close search box</title>
                <path d={CLOSE.path} fill="currentColor"/>
            </svg>
        </div>
    );
}