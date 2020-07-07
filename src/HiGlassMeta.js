import React from 'react';

import { InfoProvider } from './utils/contexts.js';
import HiGlassMetaConsumer from "./HiGlassMetaConsumer.js";

/**
 * Cistrome HiGlass Wrapper, a React component that wraps around HiGlass 
 * to provide visualization features for cistrome data.
 * @prop {object} viewConfig A HiGlass viewConfig object.
 * @prop {(object|object[])} options Options for the wrapper component.
 * @prop {function} onViewChanged A function to call upon change of the view config and option. Optional.
 * @prop {function} onGenomicIntervalSearch A function to call upon searching for TFs by using the selected interval. Optional.
 * @example
 * <HiGlassMeta
 *  viewConfig={higlassViewConfig}
 *  options={wrapperOptions}
 * />
 */
export default function HiGlassMeta(props) {
    
    console.log("HiGlassMeta.render");
    return (
        <InfoProvider>
            <HiGlassMetaConsumer {...props} />
        </InfoProvider>
    );
}