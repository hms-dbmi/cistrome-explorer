import cloneDeep from 'lodash/cloneDeep';
import uuidv4 from 'uuid/v4';

import { TRACK_TYPE } from './constants.js';
import { removeItemFromArray } from './array.js';
import { VIEWPORT_OPTIONS } from './viewport.js';

/**
 * Execute a callback function for every view, track, and innerTrack in a view config object.
 * @param {object} viewConf A valid HiGlass viewConfig object.
 * @param {function} callback A function that will be called with an object parameter with the attributes:
 *                              `viewI, viewId, trackI, trackType, trackId, innerTrackI, innerTrackType, innerTrackId`.
 *                            These object attributes will be undefined if the attribute does not apply to a particular item.
 */
export function traverseViewConfig(viewConf, callback) {
    if(viewConf && viewConf.views && Array.isArray(viewConf.views)) {
        for(let [viewI, view] of viewConf.views.entries()) {
            if(view && view.uid) {
                callback({
                    viewI,
                    viewId: view.uid
                });
                for(let [tracksPos, tracks] of Object.entries(view.tracks)) {
                    if(Array.isArray(tracks)) {
                        for(let [trackI, track] of tracks.entries()) {
                            callback({
                                viewI,
                                track,
                                viewId: view.uid,
                                trackPos: tracksPos,
                                trackI,
                                trackType: track.type,
                                trackId: track.uid,
                                trackTilesetId: track.tilesetUid,
                                trackOptions: track.options
                            });
                            if(track.type === TRACK_TYPE.COMBINED && Array.isArray(track.contents)) {
                                for(let [innerTrackI, innerTrack] of track.contents.entries()) {
                                    callback({ 
                                        viewI,
                                        viewId: view.uid,
                                        track,
                                        trackPos: tracksPos,
                                        trackI,
                                        trackType: track.type, 
                                        trackId: track.uid, 
                                        innerTrack,
                                        innerTrackI,
                                        innerTrackType: innerTrack.type, 
                                        innerTrackId: innerTrack.uid,
                                        innerTrackTilesetId: innerTrack.tilesetUid,
                                        innerTrackOptions: innerTrack.options
                                    });
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

/**
 * Get a track definition from a higlass viewConfig object.
 * @param {object} viewConf A valid HiGlass viewConfig object.
 * @param {string} viewId The uid of view containing the `horizontal-multivec` track that was the target of the action.
 * @param {string} trackId The uid of the `horizontal-multivec` track that was the target of the action. 
 * @return {object} A part of HiGlass viewConfig object for a specific track.
 */
export function getTrackDefFromViewConfig(viewConfig, viewId, trackId) {
    let newViewConfig = {};
    traverseViewConfig(cloneDeep(viewConfig), (d) => {
        // The horizontal-multivec track could be standalone, or within a "combined" track.
        if(d.trackType === TRACK_TYPE.HORIZONTAL_MULTIVEC 
            && d.viewId === viewId 
            && d.trackId === trackId
        ) {
            newViewConfig = d.track;
            return;
        } else if(d.trackType === TRACK_TYPE.COMBINED 
            && d.innerTrackType === TRACK_TYPE.HORIZONTAL_MULTIVEC 
            && d.viewId === viewId 
            && d.innerTrackId === trackId
        ) {
            newViewConfig = d.innerTrack;
            return;
        }
    });
    return newViewConfig;
}

/**
 * Add a track definitaion object to a higlass view config object.
 * @param {object} currViewConfig A valid higlass view config object.
 * @param {object} trackDef A track definition object.
 * @param {string} targetViewId The view ID for the track of interest.
 * @param {string} position The target position of the track, such as "top" or "bottom".
 * @returns {object} The new view config. 
 */
export function addTrackDefToViewConfig(currViewConfig, trackDef, targetViewId, position) {
    const newViewConfig = cloneDeep(currViewConfig);
    let viewIndex = -1;
    // Get view index.
    traverseViewConfig(currViewConfig, (d) => {
        // The horizontal-multivec track could be standalone, or within a "combined" track.
        if((d.trackType === TRACK_TYPE.HORIZONTAL_MULTIVEC 
            && d.viewId === targetViewId) || 
            (d.trackType === TRACK_TYPE.COMBINED 
            && d.innerTrackType === TRACK_TYPE.HORIZONTAL_MULTIVEC 
            && d.viewId === targetViewId)
        ) {
            viewIndex = d.viewI;
        }
    });
    if(viewIndex !== -1) {
        newViewConfig.views[viewIndex].tracks[position].push(trackDef);
    } else {
        console.log(`WARNING: The following track is not found (${targetViewId}, ${neighborTrackId}).`);
    }
    return newViewConfig;
}

/**
 * This function finds the horizontal-multivec tracks in the view config.
 * @param {object} viewConf A valid HiGlass viewConfig object.
 * @returns {array} Array containing `{ viewId, trackId, ... }` for each horizontal-multivec track.
 */
export function getHMTrackIdsFromViewConfig(viewConf) {
    const mvTracks = [];
    traverseViewConfig(viewConf, ({ viewId, trackType, trackId, trackTilesetId, innerTrackType, innerTrackId, innerTrackTilesetId }) => {
        // The horizontal-multivec track could be standalone, or within a "combined" track.
        if(trackType === TRACK_TYPE.HORIZONTAL_MULTIVEC) {
            mvTracks.push({ viewId, trackId, trackTilesetId });
        } else if(trackType === TRACK_TYPE.COMBINED && innerTrackType === TRACK_TYPE.HORIZONTAL_MULTIVEC) {
            mvTracks.push({ viewId, trackId: innerTrackId, trackTilesetId: innerTrackTilesetId, combinedTrackId: trackId });
        }
    });
    return mvTracks;
}

/**
 * This function finds the `viewport-projection-horizontal` tracks that are siblings of a particular target track.
 * @param {object} viewConf A valid HiGlass viewConfig object.
 * @param {string} targetViewId The `uid` of the view that is the parent of the track to target.
 * @returns {object[]} Array containing `{ viewId, trackId }` for each sibling `viewport-projection-horizontal` track.
 */
export function getSiblingVPHTrackIdsFromViewConfig(viewConf, targetViewId) {
    const matches = [];
    traverseViewConfig(viewConf, ({ viewId, trackPos, trackType, trackId }) => {
        if(viewId === targetViewId && trackPos === 'whole' && trackType === TRACK_TYPE.VIEWPORT_PROJECTION_HORIZONTAL) {
            matches.push({ viewId, trackId });
        }
    });
    return matches;
}

/**
 * This function removes a viewport horizontal track from a view config.
 * @param {object} viewConfig The current HiGlass view config.
 * @param {string} viewId The uid of view containing the `viewport-horizontal` track that was the target of the action.
 * @param {string} trackId The uid of the `viewport-horizontal` track that was the target of the action.
 * @returns {object} The updated HiGlass view config.
 */
export function removeViewportFromViewConfig(viewConfig, viewId, trackId) {
    const newViewConfig = cloneDeep(viewConfig);

    // Find the view associated with this viewId.
    const foundViewIndex = newViewConfig.views.findIndex(v => v.uid === viewId);
    const foundView = newViewConfig.views[foundViewIndex];

    if(!foundView.tracks['whole']) {
        // There is no viewport projection horizontal track in this view.
        console.log(`There is no ${TRACK_TYPE.VIEWPORT_PROJECTION_HORIZONTAL} in a HiGlass view.`);
        return newViewConfig;
    }
    if(foundView.tracks['whole'].find(d => d.uid === trackId)) {
        const viewportIndex = foundView.tracks['whole'].findIndex(d => d.uid === trackId);
        foundView.tracks['whole'] = removeItemFromArray(foundView.tracks['whole'], viewportIndex);
        newViewConfig.views[foundViewIndex] = foundView;
    }
    return newViewConfig;
}

/**
 * This function updates the view config when the user would like to create a genomic interval selection.
 * @param {object} viewConfig The current HiGlass view config.
 * @param {string} viewId The uid of view containing the `horizontal-multivec` track that was the target of the action.
 * @param {number} startProp The relative starting position (0 ~ 1) of the genomic interval selection from the screen.
 * @param {number} endProp The relative ending position (0 ~ 1) of the genomic interval selection from the screen.
 * @param {string} uid The uid of the parent view.
 * @returns {object} The updated HiGlass view config.
 */
export function updateViewConfigOnSelectGenomicInterval(viewConfig, viewId, startProp, endProp, uid) {
    const newViewConfig = cloneDeep(viewConfig);

    // Find the view associated with this viewId.
    const foundViewIndex = newViewConfig.views.findIndex(v => v.uid === viewId);
    const foundView = newViewConfig.views[foundViewIndex];

    const xRange = foundView.initialXDomain[1] - foundView.initialXDomain[0];
    const xDomain = [
        foundView.initialXDomain[0] + xRange * startProp,
        foundView.initialXDomain[0] + xRange * endProp
    ];

    const newUid = `${foundView.uid}-${uid}`;
    
    if(!foundView.tracks['whole']) {
        foundView.tracks['whole'] = [];
    }
    if(foundView.tracks['whole'].find(d => d.uid === newUid)) {
        const projIndex = foundView.tracks['whole'].findIndex(d => d.uid === newUid);
        foundView.tracks['whole'][projIndex] = {
            ...foundView.tracks['whole'][projIndex],
            projectionXDomain: xDomain
        };
    } else {
        const newProjectionTrackDef = {
            type: TRACK_TYPE.VIEWPORT_PROJECTION_HORIZONTAL,
            uid: newUid,
            fromViewUid: null,
            projectionXDomain: xDomain,
            options: VIEWPORT_OPTIONS["black"]
        }
        foundView.tracks['whole'].push(newProjectionTrackDef);
        newViewConfig.views[foundViewIndex] = foundView;
    }

    return newViewConfig;
}

/**
 * Set the `selectRows` option for a particular track, and return the updated view config object.
 * @param {object} currViewConfig A valid higlass view config object.
 * @param {number[]} selectedRows The array of row indices, which will become the value of the track option.
 * @param {string} targetViewId The view ID for the track of interest.
 * @param {string} targetTrackId The track ID for the track of interest.
 * @returns {object} The new view config. 
 */
export function updateViewConfigOnSelectRowsByTrack(currViewConfig, selectedRows, targetViewId, targetTrackId) {
    const newViewConfig = cloneDeep(currViewConfig);
    traverseViewConfig(currViewConfig, (d) => {
        // The horizontal-multivec track could be standalone, or within a "combined" track.
        if(d.trackType === TRACK_TYPE.HORIZONTAL_MULTIVEC 
            && d.viewId === targetViewId 
            && d.trackId === targetTrackId
        ) {
            newViewConfig.views[d.viewI].tracks[d.trackPos][d.trackI].options.selectRows = selectedRows;
        } else if(d.trackType === TRACK_TYPE.COMBINED 
            && d.innerTrackType === TRACK_TYPE.HORIZONTAL_MULTIVEC 
            && d.viewId === targetViewId 
            && d.innerTrackId === targetTrackId
        ) {
            newViewConfig.views[d.viewI].tracks[d.trackPos][d.trackI].contents[d.innerTrackI].options.selectRows = selectedRows;
        }
    });
    return newViewConfig;
}

/**
 * Get the `selectRows` option for a particular track.
 * @param {object} viewConfig A valid higlass view config object.
 * @param {string} targetViewId The view ID for the track of interest.
 * @param {string} targetTrackId The track ID for the track of interest.
 * @returns {(number[]|null)} The value of the `selectRows` option for the track.
 */
export function getHMSelectedRowsFromViewConfig(viewConfig, targetViewId, targetTrackId) {
    let selectedRows = null;
    traverseViewConfig(viewConfig, (d) => {
        // The horizontal-multivec track could be standalone, or within a "combined" track.
        if(d.trackType === TRACK_TYPE.HORIZONTAL_MULTIVEC 
            && d.viewId === targetViewId 
            && d.trackId === targetTrackId
        ) {
            selectedRows = d.trackOptions.selectRows;
        } else if(d.trackType === TRACK_TYPE.COMBINED 
            && d.innerTrackType === TRACK_TYPE.HORIZONTAL_MULTIVEC 
            && d.viewId === targetViewId 
            && d.innerTrackId === targetTrackId
        ) {
            selectedRows = d.innerTrackOptions.selectRows;
        }
    });
    return selectedRows;
}

/**
 * Get a unique viewId and trackId by referring a viewConfig.
 * @param {object} viewConfig A valid higlass view config object.
 * @param {(string|null)} baseId The base Id to make a new unique ID.
 * @param {(string|null)} idKey Either "viewId" or "trackId".
 * @param {(string|null)} interfix A prefered interfix to use for the new IDs.
 * @return {string} A new unique Id, formated as "{basedId}[-{interfix}]-{number}" or 
 *                  just a random string when baseId is not provided.
 */
export function getUniqueViewOrTrackId(viewConfig, { baseId, idKey, interfix }) {
    let newId = baseId;
    if(!baseId) {
        // Generate UID with a random string.
        newId = uuidv4();
    } else {
        if(interfix) {
            // Append suffix.
            newId = `${newId}-${interfix}`;
        }
        let isNotUnique = true, suffixNum = 1;
        const MAX_ITER = 100;
        while(isNotUnique || suffixNum > MAX_ITER) {
            // Search until we find a unique id.
            isNotUnique = false;
            traverseViewConfig(viewConfig, (d) => {
                if(d[idKey] === `${newId}-${suffixNum}`) {
                    isNotUnique = true;
                }
            })
            suffixNum++;
        }
        newId += `-${(suffixNum - 1)}`;
    }
    return newId;
}

/**
 * Get all view ID and track ID pairs of HiGlass tracks.
 * @param {object} viewConfig A valid higlass view config object.
 * @param {boolean} onlyHorizontalMultivec Only search for HiGlass `horizontal multivec tracks` (default: true).
 * @param {boolean} notSide Exclude tracks on the left or right positions (default: false).
 * @returns {array} An array of objects of {traciId, viewId}.
 */
export function getAllViewAndTrackPairs(viewConfig, options={}) {
    const {
        onlyHorizontalMultivec=true,
        notSide=false
    } = options;

    let pairs = [];
    traverseViewConfig(viewConfig, (d) => {
        // Screening tracks.
        if(!d.trackId) {
            return;
        } else if(notSide && (d.trackPos === "left" || d.trackPos === "right")) {
            return;
        }
        if(!onlyHorizontalMultivec) {
            pairs.push({
                viewId: d.viewId,
                trackId: d.trackId
            });
        } else {
            // The horizontal-multivec track could be standalone, or within a "combined" track.
            if(d.trackType === TRACK_TYPE.HORIZONTAL_MULTIVEC) {
                pairs.push({
                    viewId: d.viewId,
                    trackId: d.trackId
                });
            } else if(d.trackType === TRACK_TYPE.COMBINED && d.innerTrackType === TRACK_TYPE.HORIZONTAL_MULTIVEC) {
                pairs.push({
                    viewId: d.viewId,
                    trackId: d.innerTrackId
                });
            }
        }
    });
    return pairs;
}