import React, { createContext, useReducer } from "react";

/**
 * Helper function for constructing a reducer 
 * from an object mapping action types to handler functions.
 * See https://redux.js.org/recipes/reducing-boilerplate#reducers
 * @param {object} handlers Keys are action type strings, values are handler functions.
 * @returns {function} The reducer function.
 */
function createReducer(handlers) {
    return function reducer(state, action) {
        if (handlers.hasOwnProperty(action.type)) {
            return handlers[action.type](state, action)
        } else {
            return state
        }
    }
}

const reducer = createReducer({
    'set_row_info': (state, action) => {
        return {
            ...state,
            [action.tilesetUid]: {
                rowInfo: action.rowInfo,
                selectedRows: null,
                highlitRows: null,
                rowSort: [],
            }
        };
    },
    'set_selected_rows': (state, action) => {
        state[action.tilesetUid].selectedRows = action.selectedRows;
        return state;
    }
});

const initialState = {};
export const InfoContext = createContext(initialState);
export const InfoProvider = ({ children }) => {
    const [state, dispatch] = useReducer(reducer, initialState);

    return (
        <InfoContext.Provider value={{ state, dispatch }}>
            {children}
        </InfoContext.Provider>
    );
};