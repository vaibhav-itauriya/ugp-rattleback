import { simulate } from './model.js';

// 'run' returns full time histories; 'sweep' runs the model once per entry of
// m.patches (each a set of parameter overrides) and returns only the summaries.
// Stopping is done by the page terminating this worker.
self.onmessage = (ev) => {
  const m = ev.data;
  try {
    if (m.type === 'run') {
      const R = simulate(m.params, (f) => self.postMessage({ id: m.id, type: 'progress', f }));
      const buffers = [R.T.buffer, ...R.Y.map(a => a.buffer), ...Object.values(R.D).map(a => a.buffer)];
      self.postMessage({ id: m.id, type: 'done', result: R }, buffers);
    } else if (m.type === 'sweep') {
      const n = m.patches.length, t0 = Date.now();
      for (let i = 0; i < n; i++) {
        let row;
        try {
          const R = simulate({ ...m.base, ...m.patches[i] }, (f) => self.postMessage({ id: m.id, type: 'progress', f: (i + f) / n }));
          row = { x: m.values[i], ok: true, S: R.S, steps: R.stats.nSteps, ms: R.stats.ms };
        } catch (e) {
          row = { x: m.values[i], ok: false, error: e.message };
        }
        self.postMessage({ id: m.id, type: 'row', row, i, n });
      }
      self.postMessage({ id: m.id, type: 'sweepDone', ms: Date.now() - t0 });
    }
  } catch (err) {
    self.postMessage({ id: m.id, type: 'error', message: err.message || String(err) });
  }
};
