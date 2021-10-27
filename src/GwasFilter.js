import React, { useRef, useState, useEffect, useCallback, useMemo } from "react";
import { CLOSE, SEARCH, SQUARE_CHECK, SQUARE } from "./utils/icons.js";
import d3 from "./utils/d3.js";
import { TRAITS } from "./utils/gwas.js";
import "./TrackRowFilter.scss";

const MAX_NUM_SUGGESTIONS = 20;

export default function GwasFilter(props) {
    const {
        top,
        left,
        oneOf,
        onFilter,
        onClose,
    } = props;

    const moverRef = useRef();
    const keywordInputRef = useRef();
    const [keyword, setKeyword] = useState("");
    const [suggestionIndex, setSuggestionIndex] = useState(undefined);
    const [offset, setOffset] = useState({x: 0, y: 0});
    const dragStartPos = useRef(null);
    const keywordUpperCase = keyword.toUpperCase();
    // const [notOneOf, setNotOneOf] = useState([]); 
    
    // Styles
    const width = 180;
    const height = 30;
    const maxHeight = 420;
    const padding = 5;
    
    useEffect(() => {
        keywordInputRef.current.focus();
    });

    const suggestions = useMemo(() => {
        let result = [];
            
        const fieldDataByKeyword = TRAITS.filter(d => d.toUpperCase().includes(keywordUpperCase));
        const potentialResult = Array.from(new Set(fieldDataByKeyword));
        if(potentialResult.length < MAX_NUM_SUGGESTIONS) {
            result = potentialResult;
        } else {
            result = potentialResult.slice(0, MAX_NUM_SUGGESTIONS);
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
    }, [keyword]);

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
    }

    function onResetClick() {
        onFilter(TRAITS);
        setKeyword("");
    }

    function onFilterClose() {
        onClose();
        setKeyword("");
        setSuggestionIndex(undefined);
    }

    function onUnckeckAllClick() {
        onFilter([]);
    }

    function onFilterByKeyword() {
        let oneOfNotOneOf = keyword.toString();
        if(suggestionIndex !== undefined) {
            oneOfNotOneOf = suggestions[suggestionIndex];
        }
        onSuggestionEnter(oneOfNotOneOf);
    }

    function onSuggestionEnter(oneOfNotOneOf) {
        const isRemove = oneOf.indexOf(oneOfNotOneOf) !== -1;
        onFilter(isRemove ? oneOf.filter(d => d !== oneOfNotOneOf) : [...oneOf, oneOfNotOneOf]);
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
        case "ArrowUp":
            suggestionIndexDecrement();
            break;
        case "ArrowDown":
            suggestionIndexIncrement();
            break;
        case "Enter":
            onFilterByKeyword();
            break;
        case "Esc":
        case "Escape":
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
                    checked={oneOf.indexOf(text) !== -1}
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
            className="hm-filter"
            style={{
                display: ((left !== null && top !== null) ? "flex" : "none"),
                left: left - (width + padding * 2) / 2 + offset.x,
                top: top - (height + padding * 2) - 80 + offset.y,
            }}
        >
            <div className="hm-filter-box"
                style={{
                    padding: padding,
                    paddingLeft: "0px"
                }}>
                <div ref={moverRef} className="hm-button-drag">
                    <div/><div/><div/>
                </div>
                <svg className="hm-button-sm hm-button-static" 
                    style={{ color: "gray", marginLeft: "0px" }}
                    viewBox={SEARCH.viewBox}>
                    <path d={SEARCH.path} fill="currentColor"/>
                </svg>
                <input
                    ref={keywordInputRef}
                    className="hm-filter-box-input"
                    type="text"
                    name="default name"
                    placeholder="Search"
                    onChange={onKeywordChange}
                    onKeyDown={onKeyDown}
                    style={{ 
                        width, 
                        height 
                    }}
                />
                <svg className="hm-button-sm"
                    onClick={onResetClick} viewBox={SQUARE_CHECK.viewBox}>
                    <title>Remove all filters</title>
                    <path d={SQUARE_CHECK.path} fill="currentColor"/>
                </svg>
                <svg className="hm-button-sm"
                    onClick={onUnckeckAllClick} viewBox={SQUARE.viewBox}>
                    <title>Unckeck all categories</title>
                    <path d={SQUARE.path} fill="currentColor"/>
                </svg>
                <svg className="hm-button-sm"
                    onClick={onFilterClose} viewBox={CLOSE.viewBox}>
                    <title>Close filter component</title>
                    <path d={CLOSE.path} fill="currentColor"/>
                </svg>
            </div>   
            <div 
                className="hm-filter-suggestions"
                style={{
                    top: (padding + height),
                    left: "45px",
                    width,
                    maxHeight: maxHeight,
                    overflow: "auto",
                    visibility: suggestions.length > 0 ? "visible" : "collapse"
                }}
            >
                <ul>
                    {suggestions.map((d, i) => (
                        <li
                            key={d}
                            className={"hm-filter-suggestion-text " + (i === suggestionIndex ? "active-suggestion" : "")}
                            onClick={() => onSuggestionEnter(d)}
                            onMouseEnter={() => setSuggestionIndex(i)}
                            onMouseLeave={() => setSuggestionIndex(undefined)}
                        >
                            <SuggestionWithHighlight
                                text={d}
                                target={keywordUpperCase}
                            />
                        </li>
                    ))}
                </ul>
            </div>
        </div>
    );
}
