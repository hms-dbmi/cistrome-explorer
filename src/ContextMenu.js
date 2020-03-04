import React, { useEffect, useState } from 'react';
import PubSub from 'pubsub-js';

import { EVENT, CONTEXT_MENU_TYPE } from './utils/constants.js';

import './ContextMenu.scss';

export function destroyContextMenu() {
    PubSub.publish(EVENT.CONTEXT_MENU, {
        x: null,
        y: null,
        menuType: null,
        items: []
    });
}

/**
 * Context menu component. Subscribes to 'context-menu' event via `PubSub`.
 * @example
 * <ContextMenu/>
 */
export default function ContextMenu() {
    const [left, setLeft] = useState(null);
    const [top, setTop] = useState(null);
    const [menuItemData, setMenuItemData] = useState([]);
    
    useEffect(() => {
        const contextMenuToken = PubSub.subscribe(EVENT.CONTEXT_MENU, (msg, data) => {
            setLeft(data.x);
            setTop(data.y);
            
            let menuData = [];
            // TODO: Add common context menu items for each type here.
            switch(data.menuType) {
                case CONTEXT_MENU_TYPE.NOMINAL_BAR:
                    menuData.push({ title: data.title });
                    menuData.push({ isSeparator: true })
                    menuData.push(...data.items);
                    break;
                default:
                    break;
            }
            setMenuItemData(menuData);
        });

        return () => {
            PubSub.unsubscribe(contextMenuToken);
        };
    });

    // Function to make an item of ContextMenu
    function ContextMenuItem(props) {
        const {
            isSeparator,
            key,
            icon,
            title,
            action
        } = props;
    
        return (
            isSeparator ? 
                <hr className="chw-context-menu-separator"/>
                : <div className={action ? "chw-context-menu-item" : "chw-context-menu-item-title"} key={key}
                    onClick={action}
                    style={{ 
                        display: "flex", 
                        alignItems: "center" 
                    }}>
                    {icon ?
                        <svg className="chgw-button-sm chgw-search-button chgw-button-static"
                            viewBox={icon.viewBox}>
                            <path d={icon.path} fill="currentColor"/>
                        </svg>
                        : action ? 
                            <svg className="chgw-button-sm chgw-button-static"/>
                            : null
                    }
                    {title}
                </div>
        );
    }

    return (
        <div className="chw-context-menu-container"
            style={{
                visibility: menuItemData.length === 0 ? "hidden" : "visible",
                left, top
            }}>
            {menuItemData.map((d, i) => {
                return ContextMenuItem({...d, key: i})
            })}
        </div>
    )
}