import React, { useState } from 'react';

import { CistromeHGW } from '../index.js';

import hgDemoViewConfig1 from '../viewconfigs/horizontal-multivec.json';
import hgDemoViewConfig2 from '../viewconfigs/horizontal-multivec-local.json';

import './App.scss';

const hgViewConfigs = {
    "Demo 1 (TSV metadata)": hgDemoViewConfig1,
    "Demo 2 (JSON metadata, requires local server)": hgDemoViewConfig2
};

export default function App() {

    const [selectedKey, setSelectedKey] = useState(Object.keys(hgViewConfigs)[0]);

    return (
        <div className="app">
            <div className="header">
                <p><strong>Cistrome HiGlass Wrapper</strong></p>
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
