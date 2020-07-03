import React from 'react';

import { InfoProvider } from './utils/contexts.js';
import HiGlassWithMetadataConsumer from "./HiGlassWithMetadataConsumer.js";

/**
 * Cistrome HiGlass Wrapper, a React component that wraps around HiGlass 
 * to provide visualization features for cistrome data.
 * @prop {object} viewConfig A HiGlass viewConfig object.
 * @prop {(object|object[])} options Options for the wrapper component.
 * @prop {function} onViewConfigChange A function to call upon change of the HiGlass view config. Optional.
 * @example
 * <HiGlassWithMetadata
 *  viewConfig={higlassViewConfig}
 *  options={wrapperOptions}
 * />
 */
export default function HiGlassWithMetadata(props) {
    
    console.log("HiGlassWithMetadata.render");
    return (
        <InfoProvider>
            <HiGlassWithMetadataConsumer {...props} />
        </InfoProvider>
    );
}