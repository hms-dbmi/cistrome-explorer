/**
 * This function finds the horizontal-multivec tracks in the viewConfig.
 * @param {object} viewConf A valid HiGlass viewConfig object.
 * @returns {array} Array containing [viewId, trackId] for each horizontal-multivec track.
 */
export function getHorizontalMultivecTracksFromViewConfig(viewConf) {
    const mvTracks = [];
    if(viewConf && viewConf.views && Array.isArray(viewConf.views)) {
        for(let view of viewConf.views) {
            const viewId = view.uid;
            if(view && view.tracks) {
                for(let [tracksPos, tracks] of Object.entries(view.tracks)) {
                    if(Array.isArray(tracks)) {
                        for(let track of tracks) {
                            if(track.type === "horizontal-multivec" && track.uid) {
                                const trackId = track.uid;
                                mvTracks.push([viewId, trackId]);
                            }
                        }
                    }
                }
            }
        }
    }
    return mvTracks;
}