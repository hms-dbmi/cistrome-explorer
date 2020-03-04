import React, { useContext, useEffect, useState, useCallback } from 'react';

import { InfoContext, ACTION } from "./utils/contexts.js";
import TrackColTools from './TrackColTools.js';
import TrackRowInfo from './TrackRowInfo.js';
import TrackRowHighlight from './TrackRowHighlight.js';
import DataTable from './DataTable.js';

// TODO: remove the below fakedata import.
//       see https://github.com/hms-dbmi/cistrome-higlass-wrapper/issues/26
import fakedata from './demo/fakedata/index.js';

/**
 * Wrapper component associated with a particular HiGlass track.
 * @prop {object} options Options associated with the track. Contains values for all possible options.
 * @prop {object} multivecTrack A `horizontal-multivec` track object returned by `hgc.api.getTrackObject()`.
 * @prop {string} multivecTrackViewId The viewId for the multivecTrack.
 * @prop {string} multivecTrackTrackId The trackId for the multivecTrack.
 * @prop {(object|null)} combinedTrack A `combined` track object returned by `hgc.api.getTrackObject()`.
 *                              If not null, it is the parent track of the `multivecTrack`.
 * @prop {object[]} siblingTracks An array of `viewport-projection-horizontal` track objects, which
 *                                are siblings of `multivecTrack` (children of the same `combined` track).
 * @prop {function} onSelectGenomicInterval The function to call upon selection of a genomic interval.
 *                                          Passed down to the `TrackColTools` component.
 * @prop {function} onSortRows The function to call upon a sort interaction.
 * @prop {function} onSearchRows The function to call upon a search interaction.
 * @prop {function} onFilterRows The function to call upon a filer interaction.
 * @prop {function} drawRegister The function for child components to call to register their draw functions.
 */
export default function TrackWrapper(props) {
    const { 
        options, 
        multivecTrack,
        multivecTrackViewId,
        multivecTrackTrackId,
        combinedTrack,
        siblingTracks,
        onSelectGenomicInterval,
        onSortRows,
        onSearchRows,
        onFilterRows,
        onMetadataLoad,
        drawRegister
    } = props;

    const context = useContext(InfoContext);

    const [shouldCallOnMetadataLoad, setShouldCallOnMetadataLoad] = useState(false);

    useEffect(() => {
        if(shouldCallOnMetadataLoad) {
            onMetadataLoad();
        }
    }, [shouldCallOnMetadataLoad, multivecTrackViewId, multivecTrackTrackId]);

    if(!multivecTrack || !multivecTrack.tilesetInfo || !multivecTrack.tilesetInfo.shape) {
        // The track or track tileset info has not yet loaded.
        return null;
    }
    
    // Attributes to visualize based on the position
    const leftAttrs = options.rowInfoAttributes.filter(d => d.position === "left");
    const rightAttrs = options.rowInfoAttributes.filter(d => d.position === "right");

    const trackX = multivecTrack.position[0];
    const trackY = multivecTrack.position[1];
    const trackWidth = multivecTrack.dimensions[0];
    const trackHeight = multivecTrack.dimensions[1];
    const totalNumRows = multivecTrack.tilesetInfo.shape[1];

    // Attempt to obtain metadata values from the `tilesetInfo` field of the track.
    let rowInfo = [];
    let trackAssembly = null;
    try {
        trackAssembly = multivecTrack.tilesetInfo.coordSystem;
        // TODO: uncomment the below line to use the real metadata coming from the HiGlass Server.
        //       see https://github.com/hms-dbmi/cistrome-higlass-wrapper/issues/26
        // rowInfo = multivecTrack.tilesetInfo.row_infos.map(JSON.parse);

        // TODO: remove the below line.
        //       see https://github.com/hms-dbmi/cistrome-higlass-wrapper/issues/26
        rowInfo = fakedata[multivecTrack.id].tilesetInfo.rowInfo.slice(0, totalNumRows);
        
        if(!context.state[multivecTrackViewId] || !context.state[multivecTrackViewId][multivecTrackTrackId]) {
            context.dispatch({
                type: ACTION.SET_ROW_INFO,
                viewId: multivecTrackViewId,
                trackId: multivecTrackTrackId,
                rowInfo: rowInfo
            });
            setShouldCallOnMetadataLoad(true);
        }
    } catch(e) {
        console.log(e);
    }

    let selectedRows;
    let highlitRows;
    try {
        selectedRows = context.state[multivecTrackViewId][multivecTrackTrackId].selectedRows;
        highlitRows = context.state[multivecTrackViewId][multivecTrackTrackId].highlitRows;
    } catch(e) {
        // pass
        console.log(e);
    }

    const transformedRowInfo = (!selectedRows ? rowInfo : selectedRows.map(i => rowInfo[i]));

    // States and functions related to DataTable.
    // TODO: Move this part when we decide how to properly show the DataTable.
    const [dataTableRows, setDataTableRows] = useState(null);
    const [dataTableColumns, setDataTableColumns] = useState(null);
    function requestIntervalTfs(url) {
        // url example: "http://dbtoolkit.cistrome.org/api_interval?species=hg38&factor=tf&interval=chr6:151690496-152103274"
        const _url = "http://dbtoolkit.cistrome.org/api_interval?species=hg38&factor=tf&interval=chr6:151690496-152103274";
        console.log(_url)
        requestJSON(_url, (error, data) => {
            if(error !== null) {
                console.log("WARNING: URL not accesible " + _url)
            }
            else {
                console.log(data);
            }
        });
        setDataTableRows(rowInfo);
        setDataTableColumns([ "Species", "attr_3", "attr_4", "Cell Type", "Tissue Type", "attr_9" ]);
    }
    function requestJSON(url, callback) {
        var xhr = new XMLHttpRequest();
        xhr.open('GET', url, true);
        xhr.responseType = 'json';
        xhr.onload = function() {
          var status = xhr.status;
          if (status === 200) {
            callback(null, xhr.response);
          } else {
            callback(status, xhr.response);
          }
        };
        xhr.send();
    };
    ///

    console.log("TrackWrapper.render");
    return (
        <div className="cistrome-hgw-track-wrapper">
            {leftAttrs.length !== 0 ? 
                (<TrackRowInfo 
                    transformedRowInfo={transformedRowInfo}
                    viewId={multivecTrackViewId}
                    trackId={multivecTrackTrackId}
                    trackX={trackX}
                    trackY={trackY}
                    trackHeight={trackHeight}
                    trackWidth={trackWidth}
                    rowInfoAttributes={leftAttrs}
                    rowSort={options.rowSort}
                    rowFilter={options.rowFilter}
                    rowHighlight={options.rowHighlight}
                    rowInfoPosition="left"
                    onSortRows={onSortRows}
                    onSearchRows={onSearchRows}
                    onFilterRows={onFilterRows}
                    drawRegister={drawRegister}
                />) : null}
            {rightAttrs.length !== 0 ? 
                (<TrackRowInfo
                    transformedRowInfo={transformedRowInfo}
                    viewId={multivecTrackViewId}
                    trackId={multivecTrackTrackId}
                    trackX={trackX}
                    trackY={trackY}
                    trackHeight={trackHeight}
                    trackWidth={trackWidth}
                    rowInfoAttributes={rightAttrs}
                    rowSort={options.rowSort}
                    rowFilter={options.rowFilter}
                    rowHighlight={options.rowHighlight}
                    rowInfoPosition="right"
                    onSortRows={onSortRows}
                    onSearchRows={onSearchRows}
                    onFilterRows={onFilterRows}
                    drawRegister={drawRegister}
                />) : null}
            {options.colToolsPosition !== "hidden" ? 
                (<TrackColTools
                    trackX={trackX}
                    trackY={trackY}
                    trackHeight={trackHeight}
                    trackWidth={trackWidth}
                    trackAssembly={trackAssembly}
                    combinedTrack={combinedTrack}
                    siblingTracks={siblingTracks}
                    colToolsPosition={options.colToolsPosition}
                    onSelectGenomicInterval={onSelectGenomicInterval}
                    onRequestIntervalTFs={requestIntervalTfs}
                    drawRegister={drawRegister}
                />) : null}
            <TrackRowHighlight 
                trackX={trackX}
                trackY={trackY}
                trackHeight={trackHeight}
                trackWidth={trackWidth}
                totalNumRows={totalNumRows}
                selectedRows={selectedRows}
                highlitRows={highlitRows}
                drawRegister={drawRegister}
            />
            {dataTableRows ? 
                <DataTable
                    left={trackX}
                    top={trackY + trackHeight + 56}
                    width={trackWidth}
                    height={600}
                    rows={dataTableRows}
                    columns={dataTableColumns}
                />
                : null}
        </div>
    );
}