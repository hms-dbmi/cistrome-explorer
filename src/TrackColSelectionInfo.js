import React, { useEffect, useState } from 'react';

import { resolveIntervalCoordinates } from './utils-genome.js';

import './TrackColSelectionInfo.scss';

/**
 * Component for rendering information about a particular genomic interval selection.
 * @prop {number} trackX The track horizontal offset.
 * @prop {number} trackY The track vertical offset.
 * @prop {number} trackWidth The track width.
 * @prop {number} trackHeight The track height.
 * @prop {string} colToolsPosition The value of the `colToolsPosition` option.
 */
export default function TrackColSelectionInfo(props) {

    const {
        projectionTrack
    } = props;

    const [chrStartName, setChrStartName] = useState(null);
    const [chrStartPos, setChrStartPos] = useState(null);
    const [chrEndName, setChrEndName] = useState(null);
    const [chrEndPos, setChrEndPos] = useState(null);

    console.log(projectionTrack)

    useEffect(() => {
        const absDomain = projectionTrack.viewportXDomain;
        resolveIntervalCoordinates("hg38", absDomain[0], absDomain[1])
        .then(result => {
            setChrStartName(result[0][0]);
            setChrStartPos(result[0][1]);
            setChrEndName(result[1][0]);
            setChrEndPos(result[1][1]);
        });
    });

    return (
        <div>
            {chrStartName}:{chrStartPos} - {chrEndName}:{chrEndPos}
        </div>
    );
};