import { theme } from '../viewconfigs/horizontal-multivec-1.js';
import { REMOVE_ALLOWED_TAG_TRACKID } from '../HiGlassMetaConsumer';

export function getExampleCistromeTracks() {
	return [
		{ cid: 85057, gsm: 'GSM2695653', ct: 'Macrophage', f: 'SPI1' },
		{ cid: 73239, gsm: 'GSM1420193', ct: 'B cell precursor', f: 'IKZF1' }
	].map(d => getCistromeTrack(d));
}

export function getCistromeTrack(props) {
	const { cid, gsm, ct, f } = props;
	return {
		data: {
			type: 'cistrome-bigwig',
			cid,
			chromSizesUrl: 'https://aveit.s3.amazonaws.com/higlass/data/sequence/hg38.chrom.sizes'
		},
		uid: 'cistrome-' + cid + REMOVE_ALLOWED_TAG_TRACKID,
		type: 'gosling-track',
		options: {
			showMousePosition: true,
			mousePositionColor: 'black',
			name: `Cistrome ID ${cid} | ${gsm} | ${ct} | ${f}`,
			labelPosition: 'topLeft',
			fontSize: 12,
			labelColor: 'black',
			labelShowResolution: false,
			labelBackgroundColor: '#F6F6F6',
			labelTextOpacity: 0.6,
			labelLeftMargin: 4,
			labelRightMargin: 0,
			labelTopMargin: 2,
			labelBottomMargin: 0,
			backgroundColor: 'transparent',
			theme,
			spec: {
				data: {
					type: 'cistrome-bigwig',
					cid,
					chromSizesUrl: 'https://aveit.s3.amazonaws.com/higlass/data/sequence/hg38.chrom.sizes'
				},
				mark: 'bar',
				x: { field: 'start', type: 'genomic' },
				xe: { field: 'end', type: 'genomic' },
				y: { field: 'value', type: 'quantitative', axis: 'none' },
				color: { value: '#22908D' },
				tooltip: [
					{ field: 'start', type: 'genomic' },
					{ field: 'end', type: 'genomic' },
					{ field: 'value', type: 'quantitative' }
				],
				style: { outlineWidth: 0 },
				width: 100,
				height: 30
			}
		},
		width: 100,
		height: 40
	};
}
