import React from 'react';
import { render } from '@testing-library/react';
import App from './App';

test('renders title text', () => {
  const { getByText } = render(<App />);
  const pElement = getByText(/HiGlass/i);
  expect(pElement).toBeInTheDocument();
});
