import React, { useRef, useEffect } from 'react';
import debounce from 'lodash/debounce';
import './TrackRowZoomOverlay.scss';

/**
 * An overlay element to intercept wheel events in horizontal multivec tracks.
 * @prop {number} trackX The x offset of the horizontal multivec track.
 * @prop {number} trackY The y offset of the horizontal multivec track.
 * @prop {number} trackWidth The width of the horizontal multivec track.
 * @prop {number} trackHeight The height of the horizontal multivec track.
 * @prop {function} onZoomRows The function to call on zoom events.
 * @prop {boolean} isWheelListening Should the wheel events be intercepted right now?
 * i.e. is the designated key pressed down?
 */
export default function TrackRowZoomOverlay(props) {
	const { trackX, trackY, trackWidth, trackHeight, onZoomRows, isWheelListening } = props;

	const top = trackY;
	const left = trackX;
	const height = trackHeight;
	const width = trackWidth;

	const overlayRef = useRef();

	useEffect(() => {
		const wheelHandler = debounce(event => {
			const { deltaMode, deltaY, layerY } = event;
			onZoomRows(layerY / trackHeight, deltaY, deltaMode);
		}, 50);
		if (isWheelListening) {
			overlayRef.current.addEventListener('wheel', wheelHandler);
		}
		return () => {
			overlayRef.current.removeEventListener('wheel', wheelHandler);
		};
	}, [overlayRef, isWheelListening, onZoomRows, trackHeight]);

	// Use CSS to show or hide the overlay element
	const style = isWheelListening
		? {
				position: 'absolute',
				top: `${top}px`,
				left: `${left}px`,
				width: `${width}px`,
				height: `${height}px`,
				boxSizing: 'border-box',
				display: 'block',
				lineHeight: 1
		  }
		: {
				display: 'none'
		  };

	return (
		<div ref={overlayRef} style={style}>
			<span className="hm-vertical-zoom-info">Vertical zooming enabled</span>
		</div>
	);
}
