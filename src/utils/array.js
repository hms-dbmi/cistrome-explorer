/**
 * Insert item to an array and return it.
 * @prop {array} array Array to be updated.
 * @prop {number} index Index of array to insert new item.
 * @prop {any} item Item to be inserted.
 */
export function insertItemToArray(array, index, item) {
    return [
        ...array.slice(0, index),
        item,
        ...array.slice(index)
    ];
}

/**
 * Insert item to an array and return it.
 * @prop {array} array Array to be updated.
 * @prop {number} index Index of array to change item.
 * @prop {any} item Item to be inserted.
 */
export function modifyItemInArray(array, index, item) {
    return [
        ...array.slice(0, index),
        item,
        ...array.slice(index + 1)
    ];
}