/**
 * This function finds the horizontal-multivec tracks in the viewConfig.
 * @param {object} viewConf A valid HiGlass viewConfig object.
 * @returns {array} Array containing `[viewId, trackId, inCombined]` for each horizontal-multivec track.
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
                            if(track.type === "horizontal-multivec" && track.uid) {
                                mvTracks.push([viewId, track.uid, false]);
                            } else if(track.type === "combined" && track.uid) {
                                for(let innerTrack of track.contents) {
                                    if(innerTrack.type === "horizontal-multivec" && innerTrack.uid) {
                                        mvTracks.push([viewId, innerTrack.uid, true]);
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
}