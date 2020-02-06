import cloneDeep from 'lodash/cloneDeep';
import uuidv4 from 'uuid/v4';

import { TRACK_TYPE } from '../constants.js';

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
                                viewId: view.uid,
                                trackI,
                                trackType: track.type, 
                                trackId: track.uid
                            });
                            if(track.type === TRACK_TYPE.COMBINED && Array.isArray(track.contents)) {
                                for(let [innerTrackI, innerTrack] of track.contents.entries()) {
                                    callback({ 
                                        viewI,
                                        viewId: view.uid,
                                        trackI,
                                        trackType: track.type, 
                                        trackId: track.uid, 
                                        innerTrackI,
                                        innerTrackType: innerTrack.type, 
                                        innerTrackId: innerTrack.uid
                                    });
                                }
                            }
                        }
                    }
                }
            }
        }
    }
};

/**
 * This function finds the horizontal-multivec tracks in the view config.
 * @param {object} viewConf A valid HiGlass viewConfig object.
 * @returns {array} Array containing `[viewId, trackId, null]` or `[viewId, trackId, combinedTrackId]` for each horizontal-multivec track.
 */
export function getHMTrackIdsFromViewConfig(viewConf) {
    const mvTracks = [];
    traverseViewConfig(viewConf, ({ viewId, trackType, trackId, innerTrackType, innerTrackId }) => {
        // The horizontal-multivec track could be standalone, or within a "combined" track.
        if(trackType === TRACK_TYPE.HORIZONTAL_MULTIVEC && trackId) {
            mvTracks.push([viewId, trackId, null]);
        } else if(trackType === TRACK_TYPE.COMBINED && trackId && innerTrackType === TRACK_TYPE.HORIZONTAL_MULTIVEC && innerTrackId) {
            mvTracks.push([viewId, innerTrackId, trackId]);
        }
    });
    return mvTracks;
};

/**
 * This function finds the `viewport-projection-horizontal` tracks that are siblings of a particular target track.
 * @param {object} viewConf A valid HiGlass viewConfig object.
 * @param {string} targetTrackId The `uid` of the track to target.
 * @returns {array} Array containing `[viewId, trackId, null]` or `[viewId, trackId, combinedTrackId]` for each sibling `viewport-projection-horizontal` track.
 */
export function getSiblingProjectionTracksFromViewConfig(viewConf, targetTrackId) {
    const potentialMatches = {};
    traverseViewConfig(viewConf, ({ viewId, trackType, trackId, innerTrackType, innerTrackId }) => {
        if(trackType === TRACK_TYPE.COMBINED) {
            if(!innerTrackType) {
                potentialMatches[trackId] = {
                    match: false,
                    siblingIds: []
                };
            } else if(innerTrackType === TRACK_TYPE.VIEWPORT_PROJECTION_HORIZONTAL) {
                potentialMatches[trackId].siblingIds.push([viewId, innerTrackId, trackId]);
            } else if(innerTrackId === targetTrackId) {
                potentialMatches[trackId].match = true;
            }
        }
    });
    const matchObj = Object.values(potentialMatches).find(d => d.match);
    if(matchObj) {
        return matchObj.siblingIds;
    }
    return [];
};


/**
 * This function updates the view config when the user would like to create a genomic interval selection.
 * @param {object} currViewConfig The current HiGlass view config.
 * @param {string} viewId The uid of view containing the `horizontal-multivec` track that was the target of the action.
 * @param {string} trackId The uid of the `horizontal-multivec` track that was the target of the action.
 * @param {boolean} [addView=false] Whether or not to add the projection as a new view. If true, also updates view widths and x-offsets.
 * @returns {object} The updated HiGlass view config.
 */
export function updateViewConfigOnSelectGenomicInterval(currViewConfig, viewId, trackId, addView = false) {
    const newViewConfig = cloneDeep(currViewConfig);

    // Find the view associated with this viewId.
    const foundViewIndex = newViewConfig.views.findIndex(v => v.uid === viewId);
    const foundView = newViewConfig.views[foundViewIndex];

    if(addView) {
        foundView.layout.w = Math.max(1, (foundView.layout.w / 2) - 1);
    }
    
    // Find the track object.
    let foundTrack;
    for(let [tracksPos, tracks] of Object.entries(foundView.tracks)) {
        if(Array.isArray(tracks)) {
            for(let [i, track] of tracks.entries()) {
                if(track.type === TRACK_TYPE.HORIZONTAL_MULTIVEC && track.uid === trackId) {
                    foundTrack = track;

                    const newView = cloneDeep(foundView);
                    newView.uid = newView.uid + "-with-projection";
                    if(addView) {
                        // Need view.uid values to be unique.
                        newView.uid += uuidv4();
                    }
                    const newTrackInner = cloneDeep(foundTrack);

                    newView.tracks[tracksPos][i] = {
                        type: TRACK_TYPE.COMBINED,
                        uid: newTrackInner.uid + "-combined",
                        height: newTrackInner.height,
                        width: newTrackInner.width,
                        contents: [
                            newTrackInner,
                            {
                                uid: newTrackInner.uid + "-projection",
                                type: TRACK_TYPE.VIEWPORT_PROJECTION_HORIZONTAL,
                                fromViewUid: viewId,
                                options: {
                                    projectionFillColor: "#777",
                                    projectionStrokeColor: "#777",
                                    projectionFillOpacity: 0.3,
                                    projectionStrokeOpacity: 0.7,
                                    strokeWidth: 1
                                }
                            }
                        ]
                    };

                    newViewConfig.views[foundViewIndex] = newView;

                    if(addView) {
                        foundView.layout.x = newView.layout.w + 1;
                        newViewConfig.views.push(foundView);
                    }
                } else if(track.type === TRACK_TYPE.COMBINED) {
                    // We currently do not need to handle this case, but may want to in the future, 
                    // to allow multiple interval selections per `horizontal-multivec` track.
                }
            }
        }
    }

    return newViewConfig;
};