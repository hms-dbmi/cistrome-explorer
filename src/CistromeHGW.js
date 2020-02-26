import React from 'react';

import { InfoProvider } from './utils/contexts.js';
import CistromeHGWConsumer from "./CistromeHGWConsumer.js";

/**
 * Cistrome HiGlass Wrapper, a React component that wraps around HiGlass 
 * to provide visualization features for cistrome data.
 * @prop {object} viewConfig A HiGlass viewConfig object.
 * @prop {(object|object[])} options Options for the wrapper component.
 * @prop {function} onViewConfigChange A function to call upon change of the HiGlass view config. Optional.
 * @example
 * <CistromeHGW
 *  viewConfig={higlassViewConfig}
 *  options={wrapperOptions}
 * />
 */
export default function CistromeHGW(props) {
    
    console.log("CistromeHGW.render");
    return (
        <InfoProvider>
            <CistromeHGWConsumer {...props} />
        </InfoProvider>
    );
}