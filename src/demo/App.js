import React, { useState } from 'react';

import { CistromeHGW } from '../index.js';

import hgDemoViewConfig1 from '../viewconfigs/horizontal-multivec-1.json';
import hgDemoViewConfig2 from '../viewconfigs/horizontal-multivec-2.json';

import './App.scss';

const demos = {
    "Demo 1": {
        viewConfig: hgDemoViewConfig1,
        options: {
            rowInfoPosition: "right",
            rowLinkPosition: "left",
            colToolsPosition: "bottom",
            infoAttrPrimary: "attr_5",
            infoAttrSecondary: "attr_6"
        }
    },
    "Demo 2": {
        viewConfig: hgDemoViewConfig2,
        options: {
            rowInfoPosition: "right",
            rowLinkPosition: "left",
            colToolsPosition: "bottom",
            infoAttrPrimary: "r1",
            infoAttrSecondary: "r2"
        }
    }
};

export default function App() {

    const [selectedDemo, setSelectedDemo] = useState(Object.keys(demos)[0]);

    return (
        <div className="app">
            <div className="header">
                <h4>Cistrome HiGlass Wrapper</h4>
            </div>
            <div className="viewconf-options">
                <form>
                    {Object.keys(demos).map(vcKey => (
                        <span key={vcKey}>
                            <input 
                                type="radio" 
                                name="viewconf"
                                value={vcKey} 
                                id={vcKey} 
                                checked={selectedDemo === vcKey} 
                                onChange={e => setSelectedDemo(e.target.value)}
                            />
                            <label htmlFor={vcKey}>{vcKey}</label>
                        </span>
                    ))}
                </form>
            </div>

            <div className="container">
                <CistromeHGW 
                    viewConfig={demos[selectedDemo].viewConfig}
                    options={demos[selectedDemo].options}
                />
            </div>
        </div>
    );
};
