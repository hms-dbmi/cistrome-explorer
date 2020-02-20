import React, { useEffect, useRef, useState } from "react";
import { CLOSE, FILTER, UNDO } from './utils/icons.js';
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

    function onKeyDown(e) {
        switch(e.key){
            case 'Enter':
                onFilterClick(); 
                break;
            case 'Esc':
            case 'Escape':
                onSearchClose(); 
                break;
        }
    }

    function onFilterClick() {
        const contains = keyword.toString();
        onFilter(field, type, contains);
        setKeyword("");
    }

    function onUndoClick() {
        onFilter();
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
                top: top - (height + padding * 2) - 80,
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
                onKeyDown={onKeyDown}
                style={{ width, height }}
            />
            
            <svg className="chgw-button-sm chgw-search-button"
                onClick={onFilterClick} viewBox={FILTER.viewBox}>
                <title>Filter rows using the keyword.</title>
                <path d={FILTER.path} fill="currentColor"/>
            </svg>
            <svg className="chgw-button-sm chgw-search-button"
                onClick={onUndoClick} viewBox={UNDO.viewBox}>
                <title>Remove highlights and filters.</title>
                <path d={UNDO.path} fill="currentColor"/>
            </svg>
            <svg className="chgw-button-sm chgw-search-button"
                onClick={onSearchClose} viewBox={CLOSE.viewBox}>
                <title>Close search box</title>
                <path d={CLOSE.path} fill="currentColor"/>
            </svg>
        </div>
    );
}