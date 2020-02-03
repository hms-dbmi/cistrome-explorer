import React, { useState } from 'react';

import CistromeHGW from './CistromeHGW.js';

import hgDemoViewConfig1 from './configs/horizontal-multivec.json';
import hgDemoViewConfig2 from './configs/horizontal-multivec-local.json';

import './App.scss';

const hgViewConfigs = {
    "Old demo (with TSV metadata)": hgDemoViewConfig1,
    "New demo (with JSON metadata, local server)": hgDemoViewConfig2
};

export default function App() {

    const [selectedKey, setSelectedKey] = useState(Object.keys(hgViewConfigs)[0]);

    return (
        <div className="app">
            <div className="header">
                <p>Cistrome HiGlass Wrapper</p>
            </div>
            <div className="viewconf-options">
                <form>
                    {Object.keys(hgViewConfigs).map(vcKey => (
                        <span key={vcKey}>
                            <input 
                                type="radio" 
                                name="viewconf"
                                value={vcKey} 
                                id={vcKey} 
                                checked={selectedKey === vcKey} 
                                onChange={e => setSelectedKey(e.target.value)}
                            />
                            <label htmlFor={vcKey}>{vcKey}</label>
                        </span>
                    ))}
                </form>
            </div>

            <div className="container">
                <CistromeHGW 
                    viewConfig={hgViewConfigs[selectedKey]}
                    options={{
                        rowInfoPosition: "right",
                        rowLinkPosition: "left"
                    }}
                />
            </div>
        </div>
    );
};
