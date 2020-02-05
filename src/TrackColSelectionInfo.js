import React, { useEffect, useState } from 'react';

import { resolveIntervalCoordinates } from './utils-genome.js';

import './TrackColSelectionInfo.scss';

const ASSEMBLY = "hg38"; // TODO: obtain this from the tilesetInfo for the horizontal-multivec track.

function makeDbToolkitURL(assembly, chrStartName, chrStartPos, chrEndName, chrEndPos) {
    if(!assembly || !chrStartName || !chrStartPos || !chrEndName || !chrEndPos) {
        // We need all of these values to be able to generate the URL.
        return null;
    }
    if(chrStartName !== chrEndName) {
        // Bail out, interval spans across more than one chromosome.
        return null;
    }
    if(chrEndPos - chrStartPos > 2000000) {
        // Bail out, interval spans across more than 2 Mb, the maximum for dbtoolkit's interval search.
        return null;
    }
    return `http://dbtoolkit.cistrome.org/?specie=${assembly}&factor=tf&interval=${chrStartName}%3A${chrStartPos}-${chrEndPos}`;
}

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

    useEffect(() => {
        let didUnmount = false;
        const absDomain = projectionTrack.viewportXDomain;
        resolveIntervalCoordinates(ASSEMBLY, absDomain[0], absDomain[1])
        .then(result => {
            if(!didUnmount) {
                // Only update state if the component has not yet unmounted.
                // See https://github.com/facebook/react/issues/14369#issuecomment-468267798
                setChrStartName(result[0][0]);
                setChrStartPos(result[0][1]);
                setChrEndName(result[1][0]);
                setChrEndPos(result[1][1]);
            }
        });

        return (() => { didUnmount = true; });
    });

    const dbToolkitURL = makeDbToolkitURL(ASSEMBLY, chrStartName, chrStartPos, chrEndName, chrEndPos);

    return (
        <div>
            <div>
                {chrStartName}:{chrStartPos} - {chrEndName}:{chrEndPos}
            </div>
            {dbToolkitURL ? (
                <a 
                    href={dbToolkitURL}
                    target="_blank"
                >
                    View on DB Toolkit
                </a>
            ) : null}
        </div>
    );
};