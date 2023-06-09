import isEqual from 'lodash/isEqual';
import { processWrapperOptions } from './options';

/**
 * Find the difference between two view options.
 * @param {Object} o1 A JSON object that contains view options of `HiGlassMeta`.
 * @param {Object} o2 A JSON object that contains view options of `HiGlassMeta`.
 * @returns {boolean} Returns `true` if two view options are different, otherwise `false`.
 * TODO: Update this function to return the type of difference.
 * For example, this function can return "filter" when a filter is newly applied.
 */
export function diffViewOptions(o1, o2) {
	return !isEqual(processWrapperOptions(o1), processWrapperOptions(o2));
}
