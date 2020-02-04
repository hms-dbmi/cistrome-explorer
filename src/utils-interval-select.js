import cloneDeep from 'lodash/cloneDeep';
import uuidv4 from 'uuid/v4';

export function onSelectGenomicInterval(viewId, trackId, hgApi) {
    const newViewConfig = hgApi.getViewConfig();

    // Find the view associated with this viewId.
    const foundView = newViewConfig.views.find(v => v.uid === viewId);

    foundView.layout.w = Math.max(1, (foundView.layout.w / 2) - 1);
    
    // Find the track object.
    let foundTrack;
    for(let [tracksPos, tracks] of Object.entries(foundView.tracks)) {
        if(Array.isArray(tracks)) {
            for(let [i, track] of tracks.entries()) {
                if(track.type === "horizontal-multivec" && track.uid === trackId) {
                    foundTrack = track;

                    const newView = cloneDeep(foundView);
                    newView.uid = newView.uid + "-copy-" + uuidv4();
                    newView.layout.x = foundView.layout.w + 1;
                    const newTrackInner = cloneDeep(foundTrack);

                    newView.tracks[tracksPos][i] = {
                        type: "combined",
                        uid: newTrackInner.uid + "-combined",
                        height: newTrackInner.height,
                        width: newTrackInner.width,
                        contents: [
                            newTrackInner,
                            {
                                uid: newTrackInner + "-projection",
                                type: "viewport-projection-horizontal",
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
                    }

                    newViewConfig.views.push(newView);
                }
            }
        }
    }

    hgApi.setViewConfig(newViewConfig);
};