import { simulate } from './model.js';

self.onmessage = (ev) => {
  const { id, params } = ev.data;
  try {
    const R = simulate(params, (f) => self.postMessage({ id, type: 'progress', f }));
    const buffers = [R.T.buffer, ...R.Y.map(a => a.buffer), ...Object.values(R.D).map(a => a.buffer)];
    self.postMessage({ id, type: 'done', result: R }, buffers);
  } catch (err) {
    self.postMessage({ id, type: 'error', message: err.message || String(err) });
  }
};
