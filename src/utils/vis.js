import { TwoText, TwoRectangle } from './two.js';

/**
 * Common function for rendering the title text for a vertical visualization component.
 * @param {string} text The title text content.
 * @param {object} options Options for drawing title.
 * @param {Two} options.two The object of Two.js class.
 * @param {boolean} options.isLeft In this view placed on the left side of the track?
 * @param {number} options.width Width of this view.
 * @param {string} options.titleSuffix The suffix of a title, information about sorting and filtering status. Optional.
 * @param {string} options.backgroundColor If defined, creates a rect behind the text element. Optional. By default, "#FFF".
 */
export function drawVisTitle(text, options) {
	const { two, isLeft, width, height, titleSuffix = '', backgroundColor = '#F5F5F5' } = options;

	const margin = 4;
	const titleFontSize = 14;
	const titleLeftInitial = isLeft ? margin : width - margin;
	const titleRotate = isLeft ? -Math.PI / 2 : Math.PI / 2;
	const fullText = `${text}${titleSuffix}`;

	const title = new TwoText(titleLeftInitial, 0, width, height, fullText);
	title.fill = '#434343';
	title.fontsize = titleFontSize;
	title.align = isLeft ? 'end' : 'start';
	title.baseline = isLeft ? 'top' : 'bottom';
	title.rotation = titleRotate;

	let titleDims = two.measureText(title);
	let titleLeft = isLeft ? titleLeftInitial : titleLeftInitial - titleDims.height;
	title.x = titleLeft;

	// Experimental, place the title on the top in the horizontal direction.
	const IS_TITLE_HORIZONTAL = true;
	if (IS_TITLE_HORIZONTAL) {
		title.text = text;
		title.align = 'start';
		title.baseline = 'bottom';
		title.rotation = null;
		titleDims = two.measureText(title);
		title.x = 2;
		title.y = 30 - 6;
		title.overflow = 'ellipsis';
	}

	if (backgroundColor && !IS_TITLE_HORIZONTAL) {
		const rect = new TwoRectangle(titleLeft - 2, 0, titleDims.height + 6, titleDims.width);
		rect.fill = backgroundColor;
		rect.stroke = null;
		rect.opacity = 0.5;

		two.append(rect);
	}

	two.append(title);
}
