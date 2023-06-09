/**
 * Simple calculation on how many tracks should be rendered, including band-connection tracks
 * @prop {boolean} isLeft Is the track placed on the left side of the HiGlass?
 * @prop {boolean[]} resolveYScales An array of boolean value that represent `resolveYScale` value of a track
 */
export function getNumOfTracks(isLeft, resolveYScales) {
	// count for regular tracks
	let num = resolveYScales.length;

	// A single track can have maximum two adjacent band-connection tracks
	num += resolveYScales.filter(d => d).length * 2;

	// Remove duplicated band-connection tracks
	// (i.e., when two adjacent regular tracks uses independent y scales)
	for (let i = 0; i < resolveYScales.length - 1; i++) {
		if (resolveYScales[i] && resolveYScales[i + 1]) {
			num--;
		}
	}

	// Remove a band-connection track faraway from HiGlass because we do not want to draw that
	if ((isLeft && resolveYScales[0]) || (!isLeft && resolveYScales[resolveYScales.length - 1])) {
		num--;
	}

	return num;
}
