/*
 * TODO: Remove this folder of fake track metadata.
 * See [#26](https://github.com/hms-dbmi/cistrome-explorer/issues/26) for more info.
 */
import rowInfo1 from './cistrome-track-1/rowInfo.json';
import rowInfo2 from './cistrome-track-2/rowInfo.json';
import rowInfo3 from './cistrome-track-3/rowInfo.json';
import rowInfo10 from './cistrome-track-10/rowInfo.json';

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