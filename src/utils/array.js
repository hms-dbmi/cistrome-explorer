/**
 * Insert item to an array and return it.
 * @prop {array} array Array to be updated.
 * @prop {number} index Index of array to insert new item.
 * @prop {any} item Item to be inserted.
 * @returns Updated array.
 */
export function insertItemToArray(array, index, item) {
    if(!array) array = [];
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
 * @returns Updated array.
 */
export function modifyItemInArray(array, index, item) {
    return [
        ...array.slice(0, index),
        item,
        ...array.slice(index + 1)
    ];
}

/**
 * Remove item from an array stored in a certain index.
 * @prop {array} array Array to be updated.
 * @prop {number} index Index of an item to be removed.
 * @returns Updated array.
 */
export function removeItemFromArray(array, index) {
    return {
      item: array[index],
      array: [
        ...array.slice(0, index),
        ...array.slice(index + 1)
      ]
    };
}