import React, { useEffect, useState } from 'react';

function getCistromeTrackValue(vc, key) {
    if(vc && vc.views && vc.views.length > 0 && vc.views[0].tracks && vc.views[0].tracks.center && vc.views[0].tracks.center.length > 0 && vc.views[0].tracks.center[0].uid && vc.views[0].tracks.center[0].uid === "cistrome-track" && vc.views[0].tracks.center[0][key]) {
        return vc.views[0].tracks.center[0][key];
    }
    return null;
}

export default function CistromeGroupings(props) {
    
    const { hgApi } = props;

    const [trackServer, setTrackServer] = useState(null);
    const [tilesetUid, setTilesetUid] = useState(null);
    const [trackWidth, setTrackWidth] = useState(0);
    const [trackHeight, setTrackHeight] = useState(0);

    useEffect(() => {
        if(!hgApi) return;
        hgApi.on('viewConfig', (viewConfString) => {
            const viewConf = JSON.parse(viewConfString)
            setTrackServer(getCistromeTrackValue(viewConf, "server"));
            setTilesetUid(getCistromeTrackValue(viewConf, "tilesetUid"));
            setTrackWidth(getCistromeTrackValue(viewConf, "width"));
            setTrackHeight(getCistromeTrackValue(viewConf, "height"));

            console.log(viewConf);
        });
    }, [hgApi]);

    return (
        <div className="cistrome-hgw-groupings">

        </div>
    );
}