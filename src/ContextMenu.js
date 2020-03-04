import React, { useEffect, useState } from 'react';
import PubSub from 'pubsub-js';

import { EVENT, CONTEXT_MENU_TYPE } from './utils/constants.js';
import { FILTER } from './utils/icons.js';

import './ContextMenu.scss';

export function destroyContextMenu() {
    PubSub.publish(EVENT.CONTEXT_MENU, {
        x: null,
        y: null,
        menuType: null
    });
}

export default function ContextMenu() {

    const [left, setLeft] = useState(null);
    const [top, setTop] = useState(null);
    const [menuType, setMenuType] = useState(null);
    const [menuItemData, setMenuItemData] = useState([]);

    function ContextMenuItem(props) {
        const {
            key,
            text,
            icon
        } = props;
    
        return (
            <div className="chw-context-menu-item" key={key}
                 style={{ display: "flex", alignItems: "center" }}>
                {icon?
                    <svg className="chgw-button-sm chgw-search-button chgw-button-static"
                        viewBox={icon.viewBox}>
                        <path d={icon.path} fill="currentColor"/>
                    </svg>
                : <svg className="chgw-button-sm chgw-search-button chgw-button-static"/>}
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
        let menuData = [];
        switch(menuType) {
            case CONTEXT_MENU_TYPE.NOMINAL_BAR:
                menuData.push({text: "Highlight rows"});
                menuData.push({text: "Filter rows", icon: FILTER});
                break;
            default:
                break;
        }
        setMenuItemData(menuData);
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