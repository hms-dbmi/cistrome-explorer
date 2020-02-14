import React, { createContext, useReducer } from "react";

/*
 * Context mapping tileset UIDs to row infos.
 * Structured as follows:
 * ```
 * {
 *   "my_tileset_id": {
 * 	    rowInfo: [{}, {}],
 * 	    selectedRows: [1, 2, 3],
 * 	    highlitRows: [1, 2]
 *   }
 * }
 * ```
 */
const initialState = {};
const InfoContext = createContext(initialState);

const InfoProvider = ({ children }) => {
    const [state, dispatch] = useReducer((state, action) => {
        switch(action.type) {
            case 'set_row_info':
                return {
                    ...state,
                    [action.tilesetUid]: {
                        rowInfo: action.rowInfo,
                        selectedRows: null,
                        highlitRows: null,
                        rowSort: [],
                    }
                };
            case 'set_selected_rows':
                state[action.tilesetUid].selectedRows = action.selectedRows;
                return state;
            default:
                throw new Error();
        };
    }, initialState);

    return (
        <InfoContext.Provider value={{ state, dispatch }}>
            {children}
        </InfoContext.Provider>
    );
}

export { InfoContext, InfoProvider };