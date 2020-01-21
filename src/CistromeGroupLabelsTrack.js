
class CistromeGroupLabelsTrack {
    constructor(HGC, ...args) {
        class CistromeGroupLabelsTrackClass extends HGC.tracks.HorizontalMultivecTrack {
            constructor(context, options) {
                super(context, options);
                
                console.log(context);
            }

            draw() {
                try {
                    console.log(this.tilesetInfo.row_infos);
                    console.log(this.dimensions);
                } catch(e) {

                }
            }
    
        }
        return new CistromeGroupLabelsTrackClass(...args);
    }
};

const icon = '<svg version="1.0" xmlns="http://www.w3.org/2000/svg" width="20px" height="20px" viewBox="0 0 5640 5420" preserveAspectRatio="xMidYMid meet"> <g></g> </svg>';

// default
CistromeGroupLabelsTrack.config = {
  type: 'cistrome-group-labels',
  datatype: ['multivec'],
  local: false,
  orientation: '1d-horizontal',
  thumbnail: new DOMParser().parseFromString(icon, 'text/xml').documentElement,
  availableOptions: ['labelPosition'],
  defaultOptions: {
    labelPosition: 'outerRight',
  },
  otherOptions: {
    
  }
};


export default CistromeGroupLabelsTrack;
