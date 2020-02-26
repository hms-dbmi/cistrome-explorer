/* eslint-env node */

import React from "react";
import { render, unmountComponentAtNode } from "react-dom";
import { act } from "react-dom/test-utils";

/*
import Enzyme from 'enzyme';
import Adapter from 'enzyme-adapter-react-16';

Enzyme.configure({ adapter: new Adapter() });
*/
import CistromeHGW from "./CistromeHGW";
import hgDemoViewConfig1 from "./viewconfigs/horizontal-multivec-1.json";


let container = null;
beforeEach(() => {
  // setup a DOM element as a render target
  container = document.createElement("div");
  document.body.appendChild(container);
});

afterEach(() => {
  // cleanup on exiting
  unmountComponentAtNode(container);
  container.remove();
  container = null;
});

it("renders the CistromeHGW component", () => {
  act(() => {
    render(
        <CistromeHGW
            viewConfig={hgDemoViewConfig1}
        />,
        container
    );
  });
  expect(container.textContent).toBe("Hey, stranger");
});