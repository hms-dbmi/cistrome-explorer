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
            switch(data.menuType) {
                case CONTEXT_MENU_TYPE.NOMINAL_BAR:
                    menuData.push({ title: data.title });
                    menuData.push({ isSeparator: true });
                    menuData.push(...data.items);
                    break;
                case CONTEXT_MENU_TYPE.TREE_ANCESTOR:
                    menuData.push({ title: data.title });
                    menuData.push({ isSeparator: true });
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

    // Function to generate an item for ContextMenu
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
                <hr className="hm-context-menu-separator" key={key}/>
                : <div className={action ? "hm-context-menu-item" : "hm-context-menu-item-title"} key={key}
                    onClick={action}
                    style={{ 
                        display: "flex", 
                        alignItems: "center" 
                    }}>
                    {icon ?
                        <svg className="hm-button-sm hm-button-static"
                            viewBox={icon.viewBox}>
                            <path d={icon.path} fill="currentColor"/>
                        </svg>
                        : action ? 
                            <svg className="hm-button-sm hm-button-static"/>
                            : null
                    }
                    {title}
                </div>
        );
    }

    return (
        <div className="hm-context-menu-container"
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