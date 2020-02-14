import PubSub from 'pubsub-js';
import { EVENT } from './../constants.js';

export function destroyTooltip() {
    PubSub.publish(EVENT.TOOLTIP, {
        x: null,
        y: null,
        content: null
    });
}