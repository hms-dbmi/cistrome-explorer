/**
 * Get start and end positions inside of a predetermined range.
 * @prop {number} pos1 First position to consider.
 * @prop {number} pos2 Second position to consider.
 * @prop {number} left The position of the left side.
 * @prop {number} right The position of the right side.
 * @returns {object} Corrected positions of start and end.
 */
export function getRange(pos1, pos2, left, right) {
    let start = pos1;
    let end = pos2;

    if(end < start) {
        // Handle the oposite direction of dragging along x-axis.
        const temp = start;
        start = end;
        end = temp;
    }

    if(end < left || start > right) {
        // Suggested positions are both out of range.
        return { start: null, end: null };
    }

    // Handle when the viewport goes outside of the current view range.
    start = Math.max(start, left);
    end = Math.min(end, right);
    
    return { start, end };
}