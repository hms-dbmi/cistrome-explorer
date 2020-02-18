import React, { createContext, useReducer } from "react";
import merge from 'lodash/merge';

export const ACTION = Object.freeze({
    SET_ROW_INFO: "set_row_info",
    SELECT_ROWS: "select_rows",
    HIGHLIGHT_ROWS: "highlight_rows"
});

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
    [ACTION.SET_ROW_INFO]: (state, action) => {
        return merge(state,
            {
                [action.viewId]: {
                    [action.trackId]: {
                        rowInfo: action.rowInfo,
                        selectedRows: null,
                        highlitRows: null,
                    }
                }
            }
        );
    },
    [ACTION.SELECT_ROWS]: (state, action) => {
        if(state[action.viewId] && state[action.viewId][action.trackId]) {
            state[action.viewId][action.trackId].selectedRows = action.selectedRows;
        }
        return state;
    },
    [ACTION.HIGHLIGHT_ROWS]: (state, action) => {
        if(state[action.viewId] && state[action.viewId][action.trackId]) {
            state[action.viewId][action.trackId].highlitRows = action.highlitRows;
        }
        return state;
    }
});

/*
 * The following code is loosely based on this article:
 * https://blog.logrocket.com/use-hooks-and-context-not-react-and-redux/
 */
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