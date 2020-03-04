import React, { useEffect, useState } from 'react';
import PubSub from 'pubsub-js';

import { EVENT, CONTEXT_MENU_TYPE } from './utils/constants.js';

import './ContextMenu.scss';

export default function ContextMenu() {

    const [left, setLeft] = useState(null);
    const [top, setTop] = useState(null);
    const [menuType, setMenuType] = useState(null);
    const [menuItemData, setMenuItemData] = useState([]);

    function ContextMenuItem(props) {
        const {
            key,
            text
        } = props;
    
        return (
            <div className="chw-context-menu-item" key={key}>
                {text}
            </div>
        );
    }

    useEffect(() => {
        const contextMenuToken = PubSub.subscribe(EVENT.CONTEXT_MENU, (msg, data) => {
            setLeft(data.x);
            setTop(data.y);
            setMenuType(data.menuType);
        });

        return () => {
            PubSub.unsubscribe(contextMenuToken);
        };
    });

    useEffect(() => {
        let contextMenuItemData = [];
        switch(menuType) {
            case CONTEXT_MENU_TYPE.NOMINAL_BAR:
                contextMenuItemData.push({text: "Highlight rows"});
                contextMenuItemData.push({text: "Filter rows"});
                break;
            default:
                break;
        }
        setMenuItemData(contextMenuItemData);
    }, [left, top, menuType]);

    return (
        <div className="chw-context-menu-container"
            style={{
                left, top
            }}>
            {menuItemData.map((d, i) => {
                return ContextMenuItem({...d, key: i})
            })}
        </div>
    )
}