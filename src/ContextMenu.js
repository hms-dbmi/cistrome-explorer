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
    const [field, setField] = useState(null);
    const [value, setValue] = useState(null);
    const [onFilter, setOnFilter] = useState(null);
    const [onSearch, setOnSearch] = useState(null);
    const [menuType, setMenuType] = useState(null);
    const [menuItemData, setMenuItemData] = useState([]);

    function ContextMenuItem(props) {
        const {
            key,
            text,
            icon,
            action
        } = props;
    
        return (
            <div className="chw-context-menu-item" key={key}
                onClick={action}
                 style={{ 
                     display: "flex", 
                     alignItems: "center" 
                }}>
                {icon?
                    <svg className="chgw-button-sm chgw-search-button chgw-button-static"
                        viewBox={icon.viewBox}>
                        <path d={icon.path} fill="currentColor"/>
                    </svg>
                    : <svg className="chgw-button-sm chgw-button-static"/>
                }
                {text}
            </div>
        );
    }

    useEffect(() => {
        const contextMenuToken = PubSub.subscribe(EVENT.CONTEXT_MENU, (msg, data) => {
            setLeft(data.x);
            setTop(data.y);
            setMenuType(data.menuType);
            switch(menuType) {
                case CONTEXT_MENU_TYPE.NOMINAL_BAR:
                    setField(data.field);
                    setValue(data.value);
                    setOnFilter((f,t,v) => data.onFilter(f,t,v));
                    setOnSearch((f,t,v) => data.onSearch(f,t,v));
                    break;
                default:
                    break;
            }
            
        });

        return () => {
            PubSub.unsubscribe(contextMenuToken);
        };
    });

    useEffect(() => {
        let menuData = [];
        switch(menuType) {
            case CONTEXT_MENU_TYPE.NOMINAL_BAR:
                menuData.push({text: "Highlight rows", action: () => { onFilter(field, "nominal", value) }});
                menuData.push({text: "Filter rows", icon: FILTER, action: () => { onSearch(field, "nominal", value) }});
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