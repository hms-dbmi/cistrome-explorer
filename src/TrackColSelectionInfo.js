import React, { useEffect, useState } from 'react';
import { format as d3_format } from 'd3-format';

import { resolveIntervalCoordinates } from './utils-genome.js';
import { CISTROME_DBTOOLKIT_MAX_INTERVAL_SIZE } from './constants.js';

import './TrackColSelectionInfo.scss';

function makeDbToolkitURL(assembly, chrStartName, chrStartPos, chrEndName, chrEndPos) {
    if(chrStartName !== chrEndName) {
        // Bail out, interval spans across more than one chromosome.
        return null;
    }
    if(chrEndPos - chrStartPos > CISTROME_DBTOOLKIT_MAX_INTERVAL_SIZE) {
        // Bail out, interval is too large for dbtoolkit's interval search.
        return null;
    }
    return `http://dbtoolkit.cistrome.org/?specie=${assembly}&factor=tf&interval=${chrStartName}%3A${chrStartPos}-${chrEndPos}`;
}

const numberFormatter = d3_format(",");

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
        projectionTrack,
        trackAssembly
    } = props;

    if(!trackAssembly) {
        return null;
    }

    const [chrStartName, setChrStartName] = useState(null);
    const [chrStartPos, setChrStartPos] = useState(null);
    const [chrEndName, setChrEndName] = useState(null);
    const [chrEndPos, setChrEndPos] = useState(null);

    useEffect(() => {
        let didUnmount = false;
        const absDomain = projectionTrack.viewportXDomain;
        resolveIntervalCoordinates(trackAssembly, absDomain[0], absDomain[1])
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

    if(!chrStartName || !chrStartPos || !chrEndName || !chrEndPos) {
        return null;
    }

    const dbToolkitURL = makeDbToolkitURL(trackAssembly, chrStartName, chrStartPos, chrEndName, chrEndPos);

    return (
        <div>
            <div>
                {chrStartName}:{numberFormatter(chrStartPos)} - {chrEndName}:{numberFormatter(chrEndPos)} (selected)
            </div>
            {dbToolkitURL ? (
                <a 
                    href={dbToolkitURL}
                    target="_blank"
                >
                    Search interval on Cistrome DB Toolkit
                </a>
            ) : (
                <p className="col-selection-info-disabled">
                    Search requires interval &le; 2 Mb
                </p>
            )}
        </div>
    );
};