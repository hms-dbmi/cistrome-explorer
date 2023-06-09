/**
 * Change number (0 ~ 255) to hex string.
 * @prop {number} c Color intensity ranges from 0 to 255.
 * @returns {string} The hex string of a given intensity.
 */
export function componentToHex(c) {
	var hex = c.toString(16);
	return hex.length == 1 ? '0' + hex : hex;
}

/**
 * A function to convert from rgb to hex color.
 * @prop {number} r Red color intensity (0 ~ 255).
 * @prop {number} g Green color intensity (0 ~ 255).
 * @prop {number} b Blue color intensity (0 ~ 255).
 * @returns {string} Color in a hex string format.
 */
export function rgbToHex([r, g, b]) {
	return '#' + componentToHex(r) + componentToHex(g) + componentToHex(b);
}

/**
 * Get unique color with a number
 * @prop {number} i Index of a unique color.
 * @returns {string} The hex string of unique color
 */
export function generateNextUniqueColor(i) {
	if (i < 16777215) {
		// https://stackoverflow.com/questions/15804149/rgb-color-permutation/15804183#15804183
		return rgbToHex([i & 0xff, (i & 0xff00) >> 8, (i & 0xff0000) >> 16]);
	} else {
		console.warn('Unique colors are out of range generateNextUniqueColor().');
		return '#000000';
	}
}

// Viridis
export const DEFAULT_COLOR_RANGE = [
	'rgba(68,1,84,1)',
	'rgba(68,2,85,1)',
	'rgba(69,3,87,1)',
	'rgba(69,5,88,1)',
	'rgba(69,6,90,1)',
	'rgba(70,8,91,1)',
	'rgba(70,9,93,1)',
	'rgba(70,11,94,1)',
	'rgba(70,12,96,1)',
	'rgba(71,14,97,1)',
	'rgba(71,15,98,1)',
	'rgba(71,17,100,1)',
	'rgba(71,18,101,1)',
	'rgba(71,20,102,1)',
	'rgba(72,21,104,1)',
	'rgba(72,22,105,1)',
	'rgba(72,24,106,1)',
	'rgba(72,25,108,1)',
	'rgba(72,26,109,1)',
	'rgba(72,28,110,1)',
	'rgba(72,29,111,1)',
	'rgba(72,30,112,1)',
	'rgba(72,32,113,1)',
	'rgba(72,33,115,1)',
	'rgba(72,34,116,1)',
	'rgba(72,36,117,1)',
	'rgba(72,37,118,1)',
	'rgba(72,38,119,1)',
	'rgba(72,39,120,1)',
	'rgba(71,41,121,1)',
	'rgba(71,42,121,1)',
	'rgba(71,43,122,1)',
	'rgba(71,44,123,1)',
	'rgba(71,46,124,1)',
	'rgba(70,47,125,1)',
	'rgba(70,48,126,1)',
	'rgba(70,49,126,1)',
	'rgba(70,51,127,1)',
	'rgba(69,52,128,1)',
	'rgba(69,53,129,1)',
	'rgba(69,54,129,1)',
	'rgba(68,56,130,1)',
	'rgba(68,57,131,1)',
	'rgba(68,58,131,1)',
	'rgba(67,59,132,1)',
	'rgba(67,60,132,1)',
	'rgba(67,62,133,1)',
	'rgba(66,63,133,1)',
	'rgba(66,64,134,1)',
	'rgba(65,65,134,1)',
	'rgba(65,66,135,1)',
	'rgba(65,67,135,1)',
	'rgba(64,69,136,1)',
	'rgba(64,70,136,1)',
	'rgba(63,71,136,1)',
	'rgba(63,72,137,1)',
	'rgba(62,73,137,1)',
	'rgba(62,74,137,1)',
	'rgba(61,75,138,1)',
	'rgba(61,77,138,1)',
	'rgba(60,78,138,1)',
	'rgba(60,79,138,1)',
	'rgba(59,80,139,1)',
	'rgba(59,81,139,1)',
	'rgba(58,82,139,1)',
	'rgba(58,83,139,1)',
	'rgba(57,84,140,1)',
	'rgba(57,85,140,1)',
	'rgba(56,86,140,1)',
	'rgba(56,87,140,1)',
	'rgba(55,88,140,1)',
	'rgba(55,89,140,1)',
	'rgba(54,91,141,1)',
	'rgba(54,92,141,1)',
	'rgba(53,93,141,1)',
	'rgba(53,94,141,1)',
	'rgba(52,95,141,1)',
	'rgba(52,96,141,1)',
	'rgba(51,97,141,1)',
	'rgba(51,98,141,1)',
	'rgba(51,99,141,1)',
	'rgba(50,100,142,1)',
	'rgba(50,101,142,1)',
	'rgba(49,102,142,1)',
	'rgba(49,103,142,1)',
	'rgba(48,104,142,1)',
	'rgba(48,105,142,1)',
	'rgba(47,106,142,1)',
	'rgba(47,107,142,1)',
	'rgba(47,108,142,1)',
	'rgba(46,109,142,1)',
	'rgba(46,110,142,1)',
	'rgba(45,111,142,1)',
	'rgba(45,112,142,1)',
	'rgba(45,112,142,1)',
	'rgba(44,113,142,1)',
	'rgba(44,114,142,1)',
	'rgba(43,115,142,1)',
	'rgba(43,116,142,1)',
	'rgba(43,117,142,1)',
	'rgba(42,118,142,1)',
	'rgba(42,119,142,1)',
	'rgba(41,120,142,1)',
	'rgba(41,121,142,1)',
	'rgba(41,122,142,1)',
	'rgba(40,123,142,1)',
	'rgba(40,124,142,1)',
	'rgba(40,125,142,1)',
	'rgba(39,126,142,1)',
	'rgba(39,127,142,1)',
	'rgba(38,128,142,1)',
	'rgba(38,129,142,1)',
	'rgba(38,130,142,1)',
	'rgba(37,131,142,1)',
	'rgba(37,131,142,1)',
	'rgba(37,132,142,1)',
	'rgba(36,133,142,1)',
	'rgba(36,134,142,1)',
	'rgba(35,135,142,1)',
	'rgba(35,136,142,1)',
	'rgba(35,137,142,1)',
	'rgba(34,138,141,1)',
	'rgba(34,139,141,1)',
	'rgba(34,140,141,1)',
	'rgba(33,141,141,1)',
	'rgba(33,142,141,1)',
	'rgba(33,143,141,1)',
	'rgba(32,144,141,1)',
	'rgba(32,145,140,1)',
	'rgba(32,146,140,1)',
	'rgba(32,147,140,1)',
	'rgba(31,147,140,1)',
	'rgba(31,148,140,1)',
	'rgba(31,149,139,1)',
	'rgba(31,150,139,1)',
	'rgba(31,151,139,1)',
	'rgba(30,152,139,1)',
	'rgba(30,153,138,1)',
	'rgba(30,154,138,1)',
	'rgba(30,155,138,1)',
	'rgba(30,156,137,1)',
	'rgba(30,157,137,1)',
	'rgba(30,158,137,1)',
	'rgba(30,159,136,1)',
	'rgba(30,160,136,1)',
	'rgba(31,161,136,1)',
	'rgba(31,162,135,1)',
	'rgba(31,163,135,1)',
	'rgba(31,163,134,1)',
	'rgba(32,164,134,1)',
	'rgba(32,165,134,1)',
	'rgba(33,166,133,1)',
	'rgba(33,167,133,1)',
	'rgba(34,168,132,1)',
	'rgba(35,169,131,1)',
	'rgba(35,170,131,1)',
	'rgba(36,171,130,1)',
	'rgba(37,172,130,1)',
	'rgba(38,173,129,1)',
	'rgba(39,174,129,1)',
	'rgba(40,175,128,1)',
	'rgba(41,175,127,1)',
	'rgba(42,176,127,1)',
	'rgba(43,177,126,1)',
	'rgba(44,178,125,1)',
	'rgba(46,179,124,1)',
	'rgba(47,180,124,1)',
	'rgba(48,181,123,1)',
	'rgba(50,182,122,1)',
	'rgba(51,183,121,1)',
	'rgba(53,183,121,1)',
	'rgba(54,184,120,1)',
	'rgba(56,185,119,1)',
	'rgba(57,186,118,1)',
	'rgba(59,187,117,1)',
	'rgba(61,188,116,1)',
	'rgba(62,189,115,1)',
	'rgba(64,190,114,1)',
	'rgba(66,190,113,1)',
	'rgba(68,191,112,1)',
	'rgba(70,192,111,1)',
	'rgba(72,193,110,1)',
	'rgba(73,194,109,1)',
	'rgba(75,194,108,1)',
	'rgba(77,195,107,1)',
	'rgba(79,196,106,1)',
	'rgba(81,197,105,1)',
	'rgba(83,198,104,1)',
	'rgba(85,198,102,1)',
	'rgba(88,199,101,1)',
	'rgba(90,200,100,1)',
	'rgba(92,201,99,1)',
	'rgba(94,201,98,1)',
	'rgba(96,202,96,1)',
	'rgba(98,203,95,1)',
	'rgba(101,204,94,1)',
	'rgba(103,204,92,1)',
	'rgba(105,205,91,1)',
	'rgba(108,206,90,1)',
	'rgba(110,206,88,1)',
	'rgba(112,207,87,1)',
	'rgba(115,208,85,1)',
	'rgba(117,208,84,1)',
	'rgba(119,209,82,1)',
	'rgba(122,210,81,1)',
	'rgba(124,210,79,1)',
	'rgba(127,211,78,1)',
	'rgba(129,212,76,1)',
	'rgba(132,212,75,1)',
	'rgba(134,213,73,1)',
	'rgba(137,213,72,1)',
	'rgba(139,214,70,1)',
	'rgba(142,215,68,1)',
	'rgba(144,215,67,1)',
	'rgba(147,216,65,1)',
	'rgba(149,216,63,1)',
	'rgba(152,217,62,1)',
	'rgba(155,217,60,1)',
	'rgba(157,218,58,1)',
	'rgba(160,218,57,1)',
	'rgba(163,219,55,1)',
	'rgba(165,219,53,1)',
	'rgba(168,220,51,1)',
	'rgba(171,220,50,1)',
	'rgba(173,221,48,1)',
	'rgba(176,221,46,1)',
	'rgba(179,221,45,1)',
	'rgba(181,222,43,1)',
	'rgba(184,222,41,1)',
	'rgba(187,223,39,1)',
	'rgba(189,223,38,1)',
	'rgba(192,223,36,1)',
	'rgba(195,224,35,1)',
	'rgba(197,224,33,1)',
	'rgba(200,225,32,1)',
	'rgba(203,225,30,1)',
	'rgba(205,225,29,1)',
	'rgba(208,226,28,1)',
	'rgba(211,226,27,1)',
	'rgba(213,226,26,1)',
	'rgba(216,227,25,1)',
	'rgba(219,227,24,1)',
	'rgba(221,227,24,1)',
	'rgba(224,228,24,1)',
	'rgba(226,228,24,1)',
	'rgba(229,228,24,1)',
	'rgba(232,229,25,1)',
	'rgba(234,229,25,1)',
	'rgba(237,229,26,1)',
	'rgba(239,230,27,1)',
	'rgba(242,230,28,1)',
	'rgba(244,230,30,1)',
	'rgba(247,230,31,1)',
	'rgba(249,231,33,1)',
	'rgba(251,231,35,1)',
	'rgba(254,231,36,1)'
];
