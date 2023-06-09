import React, { forwardRef } from 'react';

import { InfoProvider } from './utils/contexts.jsx';
import HiGlassMetaConsumer from './HiGlassMetaConsumer.jsx';

/**
 * HiGlassMeta, a React component that wraps around HiGlass
 * to provide interactive metadata visualizations.
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
const HiGlassMeta = forwardRef((props, ref) => {
	// console.log("HiGlassMeta.render");
	return (
		<InfoProvider>
			<HiGlassMetaConsumer ref={ref} {...props} />
		</InfoProvider>
	);
});

export default HiGlassMeta;
