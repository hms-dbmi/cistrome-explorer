import React, { useEffect, useState } from 'react';
import PubSub from 'pubsub-js';

import './Tooltip.scss';

/**
 * Tooltip component. Subscribes to 'tooltip' event via `PubSub`.
 * @example
 * <Tooltip />
 */
export default function Tooltip() {

    const [left, setLeft] = useState(null);
    const [top, setTop] = useState(null);
    const [content, setContent] = useState("");

    useEffect(() => {
        const tooltipToken = PubSub.subscribe('tooltip', (msg, data) => {
            setLeft(data.x);
            setTop(data.y);
            setContent(data.content);
        });

        return () => {
            PubSub.unsubscribe(tooltipToken);
        };
    });

    return (
        <div>
            <div
                className="cistrome-hgw-tooltip"
                style={{
                    display: ((left !== null && top !== null) ? 'inline-block' : 'none'),
                    top: `${top}px`,
                    left: `${left}px`,
                }}
            >
                {content}
            </div>
        </div>
    );
};