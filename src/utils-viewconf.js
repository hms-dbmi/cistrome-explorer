import cloneDeep from 'lodash/cloneDeep';
import uuidv4 from 'uuid/v4';

import { TRACK_TYPE } from './constants.js';

/**
 * This function finds the horizontal-multivec tracks in the view config.
 * @param {object} viewConf A valid HiGlass viewConfig object.
 * @returns {array} Array containing `[viewId, trackId, null]` or `[viewId, trackId, combinedTrackId]` for each horizontal-multivec track.
 */
export function getTracksIdsFromViewConfig(viewConf) {
    const mvTracks = [];
    if(viewConf && viewConf.views && Array.isArray(viewConf.views)) {
        for(let view of viewConf.views) {
            const viewId = view.uid;
            if(view && view.tracks) {
                for(let [tracksPos, tracks] of Object.entries(view.tracks)) {
                    if(Array.isArray(tracks)) {
                        for(let track of tracks) {
                            // The horizontal-multivec track could be standalone, or within a "combined" track.
                            if(track.type === TRACK_TYPE.HORIZONTAL_MULTIVEC && track.uid) {
                                mvTracks.push([viewId, track.uid, null]);
                            } else if(track.type === TRACK_TYPE.COMBINED && track.uid) {
                                for(let innerTrack of track.contents) {
                                    if(innerTrack.type === TRACK_TYPE.HORIZONTAL_MULTIVEC && innerTrack.uid) {
                                        mvTracks.push([viewId, innerTrack.uid, track.uid]);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    return mvTracks;
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