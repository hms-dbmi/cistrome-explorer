import { format } from 'd3-format';

function ScaleLegendTrack(HGC, ...args) {
	if (!new.target) {
		throw new Error('Uncaught TypeError: Class constructor cannot be invoked without "new"');
	}

	// Utils
	const { colorToHex } = HGC.utils;

	class ScaleLegendTrackClass extends HGC.tracks.BarTrack {
		constructor(context, options) {
			context.dataFetcher = new EmptyDataFetcher(context.dataConfig);

			super(context, options);

			const { animate } = context;
			this.animate = animate;

			this.options = options;
			this.initOptions();
		}

		initOptions() {
			this.fontSize = +this.options.fontSize;

			this.color = colorToHex(this.options.color);
			this.background = colorToHex(this.options.background);
			this.trackStrokeWidth = this.options.trackStrokeWidth;
			this.text = this.options.text;
			this.lineStyle = this.options.lineStyle;

			this.textOptions = {
				fontSize: `${this.fontSize}px`,
				fontFamily: this.options.fontFamily,
				fontWeight: this.options.fontWeight,
				fill: this.color
			};
		}

		initTile(tile) {}
		updateTile() {}
		drawTile(tile) {}
		renderTile(tile) {}

		/*
		 * Redraw the track because the options changed
		 */
		rerender(options) {
			this.options = options;
			this.draw();
		}

		zoomed(newXScale, newYScale) {
			this.xScale(newXScale);
			this.yScale(newYScale);
			this.draw();
		}

		draw() {
			this.renderArrowAndText();
		}

		renderArrowAndText() {
			const graphics = this.pForeground;

			graphics.clear();
			graphics.removeChildren();

			// track size
			const [trackWidth, trackHeight] = this.dimensions;

			/* Track Background */
			graphics.lineStyle(this.trackStrokeWidth, this.color, 1, 0.5);
			graphics.beginFill(this.background, 1);
			graphics.drawRect(0, 0, trackWidth, trackHeight);

			/* Arrows */
			// line
			graphics.lineStyle(
				1,
				this.color,
				1,
				0.5 // alignment of the line to draw, (0 = inner, 0.5 = middle, 1 = outter)
			);

			if (this.lineStyle === 'solid') {
				graphics.moveTo(0, trackHeight / 2.0);
				graphics.lineTo(trackWidth, trackHeight / 2.0);
			} else if (this.lineStyle === 'dotted') {
				const [dash, gap] = [3, 3];
				let curX = 0;
				do {
					graphics.moveTo(curX, trackHeight / 2.0);
					graphics.lineTo(Math.min(curX + dash, trackWidth), trackHeight / 2.0);
					curX += dash + gap;
				} while (curX < trackWidth);
			} else {
				console.warn(`Not supported line style for a Legend Scale track (${this.lineStyle})`);
			}

			// triangles
			const triH = Math.min(10, trackHeight);
			const triW = triH * Math.cos((30 * Math.PI) / 180);
			const triTopY = (trackHeight - triH) / 2.0;
			const triBotY = (trackHeight + triH) / 2.0;

			// left triangle
			graphics.lineStyle(0, this.color, 1, 0.5);
			graphics.beginFill(this.color, 1);
			graphics.drawPolygon([0, trackHeight / 2.0, triW, triTopY, triW, triBotY, 0, trackHeight / 2.0]);
			graphics.endFill();

			// right triangle
			graphics.lineStyle(0, this.color, 1, 0.5);
			graphics.beginFill(this.color, 1);
			graphics.drawPolygon([
				trackWidth,
				trackHeight / 2.0,
				trackWidth - triW,
				triTopY,
				trackWidth - triW,
				triBotY,
				trackWidth,
				trackHeight / 2.0
			]);
			graphics.endFill();

			/* Text */
			// text
			const f = format('.2~s');
			const scaleNumber = this._xScale.invert(this.dimensions[0]) - this._xScale.invert(0);
			const scaleText = `${f(scaleNumber)}b`;
			const text = new HGC.libraries.PIXI.Text(scaleText, this.textOptions);
			text.interactive = true;
			text.visible = true;
			text.anchor.x = 0.5;
			text.anchor.y = 0.5;
			text.y = trackHeight / 2;
			text.x = trackWidth / 2;
			graphics.addChild(text);

			// text background
			const textStyleObj = new HGC.libraries.PIXI.TextStyle(this.textOptions);
			const { width: textWidth, height: textHeight } = HGC.libraries.PIXI.TextMetrics.measureText(
				scaleText,
				textStyleObj
			);
			const bgMargin = 4;
			graphics.lineStyle(0, this.color, 1, 0.5);
			graphics.beginFill(this.background, 1);
			graphics.drawRect(
				(trackWidth - textWidth) / 2.0 - bgMargin,
				(trackHeight - textHeight) / 2.0,
				textWidth + bgMargin * 2,
				textHeight
			);
		}

		destroyTile(tile) {
			tile.graphics.destroy();
			tile = null;
		}

		zoomedY() {}
		movedY() {}
		getMouseOverHtml() {}
		exportSVG() {}
	}
	return new ScaleLegendTrackClass(...args);
}

const icon =
	'<svg width="20" height="20" xmlns="http://www.w3.org/2000/svg"><path fill="#fff" d="M-1-1h22v22H-1z"/><g><path stroke="#007fff" stroke-width="1.5" fill="#007fff" d="M-.667-.091h5v20.167h-5z"/><path stroke-width="1.5" stroke="#e8e500" fill="#e8e500" d="M5.667.242h5v20.167h-5z"/><path stroke-width="1.5" stroke="#ff0038" fill="#ff0038" d="M15.833.076h5v20.167h-5z"/><path stroke="green" stroke-width="1.5" fill="green" d="M10.833-.258H14.5v20.167h-3.667z"/></g></svg>';

// default
ScaleLegendTrack.config = {
	type: 'scale-legend',
	datatype: ['bedlike'],
	local: false,
	orientation: '1d-horizontal',
	thumbnail: new DOMParser().parseFromString(icon, 'text/xml').documentElement,
	availableOptions: [
		'text',
		'color',
		'background',
		'trackStroke',
		'fontSize',
		'fontFamily',
		'fontWeight',
		'lineStyle'
	],
	defaultOptions: {
		text: '',
		color: '#5C5C5C',
		background: '#F6F6F6',
		trackStroke: 'black',
		trackStrokeWidth: 0,
		fontSize: 16,
		fontFamily: 'Arial',
		fontWeight: 'normal',
		lineStyle: ['solid', 'dotted'][0]
	}
};

class EmptyDataFetcher {
	constructor(dataConfig) {
		this.dataConfig = dataConfig;
	}

	tilesetInfo(callback) {
		this.tilesetInfoLoading = false;

		// Dummy values - not actually used
		const TILE_SIZE = 1024;
		const MAX_ZOOM = 22;

		const retVal = {
			tile_size: TILE_SIZE,
			bins_per_dimension: TILE_SIZE,
			max_zoom: MAX_ZOOM,
			max_width: TILE_SIZE * 2 ** MAX_ZOOM,
			min_pos: [0],
			max_pos: [3000000000]
		};
		if (callback) {
			callback(retVal);
		}
		return retVal;
	}

	fetchTilesDebounced(receivedTiles, tileIds) {
		const tiles = {};
		return tiles;
	}

	tile(z, x) {
		return this.tilesetInfo().then(tsInfo => {
			return [];
		});
	}
}

export default ScaleLegendTrack;
