import React, { useEffect, useState } from "react";
import PubSub from "pubsub-js";

import { EVENT } from "./utils/constants.js";

import "./Tooltip.scss";

export function destroyTooltip() {
    PubSub.publish(EVENT.TOOLTIP, {
        x: null,
        y: null,
        content: null
    });
}

/**
 * Publish tooltip of instructions
 * @param {Object} e Mouse event.
 * @param {string} title The title of an instruction.
 * @param {string} subtitle The detailed instruction to be shown in the tooltip.
 * @param {boolean} helpActivated Whether to show tooltip instruction or not.
 */
export function publishHelpTooltip(e, title, subtitle, helpActivated) {
    if(helpActivated && title && subtitle) {
        PubSub.publish(EVENT.TOOLTIP, {
            x: e.clientX,
            y: e.clientY,
            content: <TooltipContent 
                title={"💡 " + title}
                value={subtitle}
            />,
            help: true
        });
    }
}

/**
 * Generate tooltip contents.
 * @prop {string} title The title of tooltip.
 * @prop {number} value A value to show in tooltip.
 * @prop {string} color A color related to the value/title.
 */
export function TooltipContent(props) {
    const { title, value, color: background, warning, help } = props;
    return (
        <div className="hm-tooltip-content">
            <div className="hm-tooltip-title">{title}</div>
            <div className="hm-tooltip-content-group">
                {background ? 
                    <div className="hm-tooltip-color-container">
                        <div className="hm-tooltip-color" style={{ background }}/>
                    </div>
                    : null
                }
                {value ? 
                    <div className={warning ? "hm-tooltip-warning" : ""}>{value}</div>
                    : null
                }
            </div>
        </div>
    );
}

/**
 * Tooltip component. Subscribes to 'tooltip' event via `PubSub`.
 * @example
 * <Tooltip />
 */
export default function Tooltip() {

    const [left, setLeft] = useState(null);
    const [top, setTop] = useState(null);
    const [content, setContent] = useState("");
    const [help, setHelp] = useState(false);

    useEffect(() => {
        const tooltipToken = PubSub.subscribe(EVENT.TOOLTIP, (msg, data) => {
            setLeft(data.x);
            setTop(data.y);
            setContent(data.content);
            setHelp(data.help);
        });

        return () => {
            PubSub.unsubscribe(tooltipToken);
        };
    });

    return (
        <div
            className={"hm-tooltip " + (help ? "hm-tooltip-help" : "")}
            style={{
                display: ((left !== null && top !== null) ? "inline-block" : "none"),
                top: `${Math.min(document.body.scrollHeight - 100, top)}px`,
                left: `${Math.min(document.body.scrollWidth - 200, left)}px`,
            }}
        >
            {content}
        </div>
    );
}