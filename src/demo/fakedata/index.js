/*
 * TODO: Remove this folder of fake track metadata.
 * See [#26](https://github.com/hms-dbmi/cistrome-explorer/issues/26) for more info.
 */
import rowInfo1 from "./cistrome-track-1/rowInfo.json";
import rowInfo2 from "./cistrome-track-2/rowInfo.json";
import rowInfo3 from "./cistrome-track-3/rowInfo.json";
import rowInfo10 from "./cistrome-track-10/rowInfo.json";
import rowInfoAtac from "./cistrome-track-atac/rowInfoWithQc.json";
import rowInfoAtacRevision from "./cistrome-track-atac/rowInfoRevision.json";
import rowInfo3k27 from "./cistrome-track-3k27/rowInfoWithQc.json";
import rowInfo3k27Revision from "./cistrome-track-3k27/rowInfoRevision.json";
import rowInfo3k4 from './cistrome-track-3k4/rowInfo.json';
import rowInfoMiraMouse from './mira-track-mouse/rowInfo.json';

const _fakedata = {
    "default": {
        tilesetInfo: {
            rowInfo: rowInfo1
        }
    },
    "cistrome-track-1": {
        tilesetInfo: {
            rowInfo: rowInfo1
        }
    },
    "cistrome-track-atac": {
        tilesetInfo: {
            rowInfo: rowInfoAtac
        }
    },
    "cistrome-track-atac-revision": {
        tilesetInfo: {
            rowInfo: rowInfoAtacRevision
        }
    },
    "cistrome-track-3k27": {
        tilesetInfo: {
            rowInfo: rowInfo3k27
        }
    },
    "cistrome-track-3k27-revision": {
        tilesetInfo: {
            rowInfo: rowInfo3k27Revision
        }
    },
    "cistrome-track-3k4": {
        tilesetInfo: {
            rowInfo: rowInfo3k4
        }
    },
    "mira-track-mouse": {
        tilesetInfo: {
            // TODO:
            rowInfo: rowInfoMiraMouse
        }
    },
    "cistrome-track-1b": {
        tilesetInfo: {
            rowInfo: rowInfo1
        }
    },
    "cistrome-track-2": {
        tilesetInfo: {
            rowInfo: rowInfo2
        }
    },
    "cistrome-track-2b": {
        tilesetInfo: {
            rowInfo: rowInfo3
        }
    },
    "cistrome-track-4-1": {
        tilesetInfo: {
            rowInfo: rowInfo2
        }
    },
    "cistrome-track-4-2": {
        tilesetInfo: {
            rowInfo: rowInfo2
        }
    },
    "cistrome-track-6-1": {
        tilesetInfo: {
            rowInfo: rowInfo1
        }
    },
    "cistrome-track-6-2": {
        tilesetInfo: {
            rowInfo: rowInfo1
        }
    },
    "cistrome-track-7": {
        tilesetInfo: {
            rowInfo: rowInfo1
        }
    },
    "cistrome-track-8-1": {
        tilesetInfo: {
            rowInfo: rowInfo1
        }
    },
    "cistrome-track-8-2": {
        tilesetInfo: {
            rowInfo: rowInfo2
        }
    },
    "cistrome-track-10": {
        tilesetInfo: {
            rowInfo: rowInfo10
        }
    },
};

export default new Proxy(_fakedata, {
    get: function(data, key) {
        // Get default data when key is not provided.
        return data.hasOwnProperty(key) ? data[key] : data["default"];
    }
});
