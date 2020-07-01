import React, { useRef, useState, useEffect, useCallback, useMemo } from "react";
import { CLOSE, FILTER, RESET, SEARCH, SQUARE_CHECK, SQUARE } from './utils/icons.js';
import d3 from "./utils/d3.js";

import "./TrackRowFilter.scss";
import RangeSlider from "./RangeSlider.js";
import { getAggregatedValue } from "./utils/aggregate.js";

const MAX_NUM_SUGGESTIONS = 200;

/**
 * Component to determine rows to filter.
 * @prop {boolean} isLeft Is this view on the left side of the HiGlass track?
 * @prop {number} top The top coordinate.
 * @prop {number} left The left coordinate.
 * @prop {string} field The name of field related to the wrapper track.
 * @prop {string} type The type of field related to the wrapper track.
 * @prop {string} aggFunction The function to apply upon row aggregation.
 * @prop {function} onChange The function to call when the search keyword or range has changed.
 * @prop {function} onFilterRows The function to call when the filter should be applied.
 * @prop {function} onClose The function to call when the search field should be closed.
 * @prop {object[]} rowInfo The array of JSON Object containing row information.
 * @prop {object} filterInfo The options for filtering rows of the field used in this track.
 * @example
 * <TrackRowFilter/>
 */
export default function TrackRowFilter(props) {

    const {
        isLeft,
        top, left,
        field, type,
        aggFunction,
        onChange,
        onFilterRows,
        onClose,
        rowInfo,
        filterInfo
    } = props;

    const moverRef = useRef();
    const keywordInputRef = useRef();
    const [keyword, setKeyword] = useState("");
    const [suggestionIndex, setSuggestionIndex] = useState(undefined);
    const [offset, setOffset] = useState({x: 0, y: 0});
    const dragStartPos = useRef(null);
    const [valueExtent, setValueExtent] = useState(d3.extent(rowInfo.map(d => getAggregatedValue(d, field, type, aggFunction))));
    const cutoffRange = useRef(valueExtent);
    const keywordUpperCase = keyword.toUpperCase();
    const [notOneOf, setNotOneOf] = useState(!filterInfo || type === "quantitative" ? [] : filterInfo.notOneOf); 
    
    // Styles
    const width = 180;
    const height = 30;
    const maxHeight = 420;
    const padding = 5;
    
    useEffect(() => {
        setValueExtent(d3.extent(rowInfo.map(d => getAggregatedValue(d, field, type, aggFunction))));
    }, [field, rowInfo]);

    useEffect(() => {
        if(type === "nominal" || type === "link") {
            keywordInputRef.current.focus();
        }
    });

    useEffect(() => {
        setNotOneOf(!filterInfo || type === "quantitative" ? [] : filterInfo.notOneOf);
    }, [filterInfo])

    const suggestions = useMemo(() => {
        let result = [];
        if(!Array.isArray(field)) {
            const fieldData = rowInfo.map(d => getAggregatedValue(d, field, type, aggFunction).toString());
            const fieldDataByKeyword = fieldData.filter(d => d.toUpperCase().includes(keywordUpperCase));
            const potentialResult = Array.from(new Set(fieldDataByKeyword));
            if(potentialResult.length < MAX_NUM_SUGGESTIONS) {
                result = potentialResult;
            } else {
                result = potentialResult.slice(0, MAX_NUM_SUGGESTIONS);
            }
        }
        // Sort in alphabetical order first.
        result.sort((a, b) => {
            return a.toUpperCase() > b.toUpperCase();
        });
        // Sort so that suggestions that _start with_ the keyword appear first.
        result.sort((a, b) => {
            return a.toUpperCase().indexOf(keywordUpperCase) - b.toUpperCase().indexOf(keywordUpperCase);
        });
        return result;
    }, [field, keyword]);

    useEffect(() => {
        if(suggestionIndex === undefined && keyword.length > 0) {
            onChange(keyword);
        }
        if(suggestionIndex !== undefined && suggestionIndex >= 0 && suggestionIndex < suggestions.length) {
            onChange(suggestions[suggestionIndex]);
        }
    }, [suggestions, suggestionIndex]);

    // Set up the d3-drag handler functions (started, ended, dragged).
    const started = useCallback(() => {
        const event = d3.event;
        dragStartPos.current = {
            x: event.sourceEvent.clientX, 
            y: event.sourceEvent.clientY
        };
    });

    const ended = useCallback(() => {
        dragStartPos.current = null;
    });

    const dragged = useCallback(() => {
        const event = d3.event;
        const diffX = offset.x + event.sourceEvent.clientX - dragStartPos.current.x;
        const diffY = offset.y + event.sourceEvent.clientY - dragStartPos.current.y;
        setOffset({x: diffX, y: diffY});
    });

    // Detect drag events for the resize element.
    useEffect(() => {
        const mover = moverRef.current;

        const drag = d3.drag()
            .on("start", started)
            .on("drag", dragged)
            .on("end", ended);

        d3.select(mover).call(drag);

        return () => d3.select(mover).on(".drag", null);
    }, [moverRef, started, dragged, ended]);
    
    function onKeywordChange(e) {
        const newKeyword = e.target.value;
        setKeyword(newKeyword);
        onChange(newKeyword);
    }

    function onRangeChange(range) {
        const [left, right] = range;
        cutoffRange.current = left < right ? range : Array.from(range).reverse();
        onChange(cutoffRange.current);
    }

    function onResetClick() {
        if(type === "nominal" || type == "link") {
            onFilterRows(field, type, rowInfo.map(d => getAggregatedValue(d, field, type, aggFunction).toString()), true);
        } else if(type === "quantitative") {
            onFilterRows(field, type, [], true);
        }
        setKeyword("");
    }

    function onFilterClose() {
        onChange("");
        onClose();
        setKeyword("");
        setSuggestionIndex(undefined);
    }

    function onFilterClick() {
        if(type === "nominal" || type === "link") {
            onFilterByKeyword();
        } else if(type === "quantitative") {
            onFilterByRange(cutoffRange.current);
        }
    }

    function onUnckeckAllClick() {
        onFilterRows(field, type, rowInfo.map(d => getAggregatedValue(d, field, type, aggFunction).toString()), false);
    }

    function onFilterByKeyword() {
        let oneOfNotOneOf = keyword.toString();
        if(suggestionIndex !== undefined) {
            oneOfNotOneOf = suggestions[suggestionIndex];
        }
        onSuggestionEnter(oneOfNotOneOf);
    }

    function onFilterByRange(range) {
        const [left, right] = range;
        cutoffRange.current = left < right ? range : Array.from(range).reverse();
        onFilterRows(field, type, cutoffRange.current);
    }

    function onSuggestionEnter(oneOfNotOneOf) {
        const isRemove = notOneOf.indexOf(oneOfNotOneOf) !== -1;
        onFilterRows(field, type, [oneOfNotOneOf], isRemove);
        setKeyword("");
    }

    function suggestionIndexIncrement() {
        let newIndex;
        if(suggestionIndex === undefined) {
            newIndex = 0; // start at the first suggestion
        } else {
            if(suggestionIndex+1 > suggestions.length-1) {
                newIndex = undefined; // out of range above
            } else {
                newIndex = suggestionIndex+1;
            }
        }
        setSuggestionIndex(newIndex);
    }

    function suggestionIndexDecrement() {
        let newIndex;
        if(suggestionIndex === undefined) {
            newIndex = suggestions.length-1; // start at the last suggestion
        } else {
            if(suggestionIndex-1 < 0) {
                newIndex = undefined; // out of range below
            } else {
                newIndex = suggestionIndex-1;
            }
        }
        setSuggestionIndex(newIndex);
    }

    function onKeyDown(e) {
        switch(e.key){
            case 'ArrowUp':
                suggestionIndexDecrement();
                break;
            case 'ArrowDown':
                suggestionIndexIncrement();
                break;
            case 'Enter':
                onFilterByKeyword();
                break;
            case 'Esc':
            case 'Escape':
                onFilterClose(); 
                break;
        }
    }

    /**
     * Returns <span> elements in which text is highlighted based on a keyword
     * @prop {string} text The suggested search text.
     * @prop {string} target The keyword to highlight, uppercase.
     */
    function SuggestionWithHighlight(props) {
        const {
            text,
            target
        } = props;
        const i0 = text.toUpperCase().indexOf(target);
        const i1 = i0 + target.length;

        const s0 = text.substring(0, i0);
        const s1 = text.substring(i0, i1);
        const s2 = text.substring(i1, text.length);
        return (
            <div>
                <input
                    style={{
                        display: "inline-block",
                        verticalAlign: "top",
                        width: "30px",
                        marginLeft: "5px",
                        marginTop: "0px"
                    }}
                    type="checkbox"
                    name="data-table-checkbox"
                    value={text}
                    checked={notOneOf.indexOf(text) === -1}
                    readOnly={true} // This checkbox only provide visual feedback.
                />
                <div
                    style={{
                        width: "100px",
                        display: "inline-block"
                    }}
                >
                    {s0}<b>{s1}</b>{s2}
                </div>
            </div>
        );
    }

    return (
      <div
        className="cisvis-filter"
        style={{
          display: left !== null && top !== null ? "flex" : "none",
          left: left - (width + padding * 2) / 2 + offset.x,
          top: top - (height + padding * 2) - 80 + offset.y,
        }}
      >
        <div
          className="cisvis-filter-box"
          style={{
            padding: padding,
            paddingLeft: "0px",
          }}
        >
          <div ref={moverRef} className="chw-button-drag">
            <div />
            <div />
            <div />
          </div>
          {type === "nominal" || type === "link" ? (
            <svg
              className="chw-button-sm chw-button-static"
              style={{ color: "gray", marginLeft: "0px" }}
              viewBox={SEARCH.viewBox}
            >
              <path d={SEARCH.path} fill="currentColor" />
            </svg>
          ) : null}
          {type === "nominal" || type === "link" ? (
            <input
              ref={keywordInputRef}
              className="cisvis-filter-box-input"
              type="text"
              name="default name"
              placeholder="Search"
              onChange={onKeywordChange}
              onKeyDown={onKeyDown}
              style={{
                width,
                height,
              }}
            />
          ) : (
            <RangeSlider
              isRight={!isLeft}
              height={height}
              valueExtent={valueExtent}
              onChange={onRangeChange}
              onFilter={onFilterByRange}
              onClose={onFilterClose}
            />
          )}
          {type === "quantitative" ? (
            <svg
              className="chw-button-sm"
              onClick={onFilterClick}
              viewBox={FILTER.viewBox}
            >
              <title>Filter rows by the range of values</title>
              <path d={FILTER.path} fill="currentColor" />
            </svg>
          ) : null}
          <svg
            className="chw-button-sm"
            onClick={onResetClick}
            viewBox={
              type === "quantitative" ? RESET.viewBox : SQUARE_CHECK.viewBox
            }
          >
            <title>Remove all filters</title>
            <path
              d={type === "quantitative" ? RESET.path : SQUARE_CHECK.path}
              fill="currentColor"
            />
          </svg>
          {type === "nominal" || type === "link" ? (
            <svg
              className="chw-button-sm"
              onClick={onUnckeckAllClick}
              viewBox={SQUARE.viewBox}
            >
              <title>Unckeck all categories</title>
              <path d={SQUARE.path} fill="currentColor" />
            </svg>
          ) : null}
          <svg
            className="chw-button-sm"
            onClick={onFilterClose}
            viewBox={CLOSE.viewBox}
          >
            <title>Close filter component</title>
            <path d={CLOSE.path} fill="currentColor" />
          </svg>
        </div>
        {type === "nominal" || type === "link" ? (
          <div
            className="cisvis-filter-suggestions"
            style={{
              top: padding + height,
              left: "45px",
              width,
              maxHeight: maxHeight,
              overflow: "auto",
              visibility: suggestions.length > 0 ? "visible" : "collapse",
            }}
          >
            <ul>
              {suggestions.map((d, i) => (
                <li
                  key={d}
                  className={
                    "cisvis-filter-suggestion-text " +
                    (i === suggestionIndex ? "active-suggestion" : "")
                  }
                  onClick={() => onSuggestionEnter(d)}
                  onMouseEnter={() => setSuggestionIndex(i)}
                  onMouseLeave={() => setSuggestionIndex(undefined)}
                >
                  <SuggestionWithHighlight text={d} target={keywordUpperCase} />
                </li>
              ))}
            </ul>
          </div>
        ) : null}
      </div>
    );
}