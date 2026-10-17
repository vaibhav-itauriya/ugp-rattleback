import uPlot from '../vendor/uPlot.esm.js';
import { DEFAULTS } from './model.js';

const $ = (id) => document.getElementById(id);
const css = (name) => getComputedStyle(document.documentElement).getPropertyValue(name).trim();
const reduceMotion = matchMedia('(prefers-reduced-motion: reduce)').matches;

// ------------------------------------------------------------------ presets
const PRESETS = [
  { id: 'coulomb', name: 'Coulomb friction, μ = 0.15', note: 'The paper\'s main case: reverses at 4.22 s and loses 12% of its energy.',
    p: { law: 'coulomb', mu: 0.15 } },
  { id: 'rolling', name: 'Rolling without slipping', note: 'No energy is lost: the spin comes back at full speed after reversing at 2.83 s.',
    p: { law: 'rolling' } },
  { id: 'drill', name: 'Coulomb + drilling, ρ = 0.5 mm', note: 'Reverses at 3.98 s, then the drilling torque brings it nearly to rest.',
    p: { law: 'coulomb', mu: 0.15, drill: true, rho: 0.5 } },
  { id: 'viscous', name: 'Viscous friction, c = 500 kg/s', note: 'Reverses at 2.90 s, close to rolling.',
    p: { law: 'viscous', cv: 500 } },
  { id: 'weakvisc', name: 'Weak viscous friction, c = 1.5 kg/s', note: 'Same nominal size as μ = 0.15, yet the spin never reverses.',
    p: { law: 'viscous', cv: 1.5 } },
  { id: 'bigpatch', name: 'Oversized contact patch, ρ = 2 mm', note: 'Drilling stops the body at about 3 s, before the rocking can grow: a spin-down, not a reversal.',
    p: { law: 'coulomb', mu: 0.15, drill: true, rho: 2 } },
  { id: 'fixn', name: 'Shortcut N = mg, μ = 1', note: 'Reverses at 2.43 s instead of the correct 2.83 s. Takes several seconds to compute.',
    p: { law: 'coulomb', mu: 1, nmode: 'fixed' } },
  { id: 'coarse', name: 'Unconverged ε = 10⁻³ m/s', note: 'Reverses at 4.85 s, about 15% later than the converged 4.22 s.',
    p: { law: 'coulomb', mu: 0.15, eps: '1e-3' } },
  { id: 'pref', name: 'Preferred sense, Coulomb', note: 'Spun the other way: no rocking, no slip, no reversal.',
    p: { law: 'coulomb', mu: 0.15, sense: '1' } },
];
const BASE_FORM = { law: 'coulomb', mu: 0.15, eps: '3e-5', drill: false, rho: 0.5, cv: 500, nmode: 'solved',
  sense: '-1', w0: 5, tEnd: 8, tilt0: 0, m: 1, a1: 0.2, a2: 0.03, a3: 0.02, b: 0.005, I11: 2, I22: 16, I33: 17.5, I12: -0.5 };
const NUM = ['mu', 'rho', 'cv', 'w0', 'tEnd', 'tilt0', 'm', 'a1', 'a2', 'a3', 'b', 'I11', 'I22', 'I33', 'I12'];

// ------------------------------------------------------------------ form
function readForm() {
  const f = {};
  f.law = document.querySelector('input[name=law]:checked').value;
  f.sense = document.querySelector('input[name=sense]:checked').value;
  f.eps = $('eps').value; f.nmode = $('nmode').value; f.drill = $('drill').checked;
  for (const k of NUM) f[k] = parseFloat($(k).value);
  return f;
}
function writeForm(f) {
  document.querySelector(`input[name=law][value=${f.law}]`).checked = true;
  document.querySelector(`input[name=sense][value="${f.sense}"]`).checked = true;
  $('eps').value = f.eps; $('nmode').value = f.nmode; $('drill').checked = !!f.drill;
  for (const k of NUM) $(k).value = f[k];
  syncSliders(); updateVisibility();
}
function toParams(f) {
  for (const k of NUM) if (!Number.isFinite(f[k])) throw new Error(`Enter a number for ${$(k).labels[0].textContent.trim()}.`);
  if (f.mu <= 0 && f.law === 'coulomb') throw new Error('The friction coefficient must be above zero.');
  if (f.tEnd <= 0 || f.tEnd > 30) throw new Error('Choose a simulated time between 1 and 30 s.');
  if (f.a1 <= 0 || f.a2 <= 0 || f.a3 <= 0 || f.m <= 0) throw new Error('Mass and semi-axes must be positive.');
  if (f.b >= f.a3) throw new Error('The centre offset b must be smaller than a₃, or the body is not stable upright.');
  // the preferred sense of this body is +w3 when I12 < 0 (Kane and Levinson's sign)
  const pref = f.I12 < 0 ? 1 : f.I12 > 0 ? -1 : 1;
  // tilt 0 = the paper's start: 1e-6 rad on each Euler angle
  const tilt = f.tilt0 > 0 ? f.tilt0 * Math.PI / 180 / Math.SQRT2 : 1e-6;
  return {
    law: f.law, mu: f.mu, eps: parseFloat(f.eps), rho: f.law === 'coulomb' && f.drill ? f.rho * 1e-3 : 0,
    cv: f.cv, fixN: f.law !== 'rolling' && f.nmode === 'fixed',
    w0: Math.abs(f.w0) * (f.sense === '1' ? pref : -pref), tEnd: f.tEnd,
    phi0: tilt, theta0: tilt, psi0: 1e-6,
    m: f.m, a1: f.a1, a2: f.a2, a3: f.a3, b: f.b,
    I11: f.I11 * 1e-4, I22: f.I22 * 1e-4, I33: f.I33 * 1e-4, I12: f.I12 * 1e-4,
  };
}
function describe(f) {
  const s = f.sense === '1' ? 'preferred sense' : 'unstable sense';
  const w = `|ω₀| = ${f.w0} rad/s`;
  if (f.law === 'rolling') return `Rolling, ${s}, ${w}`;
  const nm = f.nmode === 'fixed' ? ', N = mg' : '';
  if (f.law === 'viscous') return `Viscous c = ${f.cv} kg/s, ${s}, ${w}${nm}`;
  const dr = f.drill ? ` + drilling ρ = ${f.rho} mm` : '';
  return `Coulomb μ = ${f.mu}${dr}, ${s}, ${w}${nm}`;
}

function updateVisibility() {
  const f = { law: document.querySelector('input[name=law]:checked').value, drill: $('drill').checked };
  document.querySelectorAll('[data-for]').forEach(el => {
    const k = el.dataset.for;
    el.hidden = !(k === f.law || (k === 'drill' && f.law === 'coulomb' && f.drill) || (k === 'slip' && f.law !== 'rolling'));
  });
}
// paired number + slider inputs (cv slider is logarithmic)
const PAIRS = [['mu', 'mu-r', 'lin'], ['rho', 'rho-r', 'lin'], ['cv', 'cv-r', 'log'], ['w0', 'w0-r', 'lin']];
function syncSliders() {
  for (const [n, r, kind] of PAIRS) { const v = parseFloat($(n).value); if (Number.isFinite(v)) $(r).value = kind === 'log' ? Math.log10(Math.max(v, 1)) : v; }
}
for (const [n, r, kind] of PAIRS) {
  $(r).addEventListener('input', () => {
    const v = parseFloat($(r).value);
    $(n).value = kind === 'log' ? +Math.pow(10, v).toPrecision(2) : v;
    markCustom();
  });
  $(n).addEventListener('input', () => { syncSliders(); });
}
$('controls').addEventListener('input', (e) => { if (e.target.id !== 'preset') markCustom(); updateVisibility(); });

// presets
const presetSel = $('preset');
presetSel.innerHTML = PRESETS.map(p => `<option value="${p.id}">${p.name}</option>`).join('') + '<option value="custom">Your own settings</option>';
function applyPreset(id) {
  const P = PRESETS.find(p => p.id === id); if (!P) return;
  writeForm({ ...BASE_FORM, ...P.p });
  presetSel.value = id; $('presetNote').textContent = P.note;
}
function markCustom() { if (presetSel.value !== 'custom') { presetSel.value = 'custom'; $('presetNote').textContent = 'Press Run simulation to see the effect of your changes.'; } }
presetSel.addEventListener('change', () => { if (presetSel.value !== 'custom') { applyPreset(presetSel.value); run(); } });

// share link
function toHash(f) { return '#' + new URLSearchParams(Object.entries(f).map(([k, v]) => [k, String(v)])).toString(); }
function fromHash() {
  if (location.hash.length < 2) return null;
  const q = new URLSearchParams(location.hash.slice(1)), f = { ...BASE_FORM };
  for (const [k, v] of q) if (k in f) f[k] = NUM.includes(k) ? parseFloat(v) : k === 'drill' ? v === 'true' : v;
  return f;
}
$('linkBtn').addEventListener('click', async () => {
  const url = location.href.split('#')[0] + toHash(readForm());
  history.replaceState(null, '', url);
  try { await navigator.clipboard.writeText(url); setStatus('Link copied. It opens this page with the same settings.'); }
  catch { setStatus('Link is in the address bar; copy it from there.'); }
});

// theme
$('themeBtn').addEventListener('click', () => {
  const cur = document.documentElement.dataset.theme || (matchMedia('(prefers-color-scheme: dark)').matches ? 'dark' : 'light');
  const next = cur === 'dark' ? 'light' : 'dark';
  document.documentElement.dataset.theme = next;
  try { localStorage.setItem('rb-theme', next); } catch {}
  stage && stage.setTheme(next === 'dark'); rebuildCharts();
});
try { const t = localStorage.getItem('rb-theme'); if (t) document.documentElement.dataset.theme = t; } catch {}

// ------------------------------------------------------------------ worker
let worker = null, jobId = 0, running = false;
function newWorker() {
  worker = new Worker(new URL('./worker.js', import.meta.url), { type: 'module' });
  worker.onmessage = onWorker;
  worker.onerror = (e) => { finishRun(); setStatus('The simulation could not start in this browser: ' + (e.message || 'worker error'), true); };
}
function setStatus(msg, err = false) { const s = $('status'); s.textContent = msg; s.classList.toggle('error', err); }
let pendingForm = null, t0 = 0;
function run() {
  if (running) { cancel(); return; }
  let f, p;
  try { f = readForm(); p = toParams(f); } catch (e) { setStatus(e.message, true); return; }
  if (!worker) newWorker();
  running = true; pendingForm = f; t0 = performance.now();
  $('runBtn').textContent = 'Cancel'; $('bar').style.width = '0%';
  setStatus('Integrating the equations of motion…');
  worker.postMessage({ id: ++jobId, params: p });
}
function cancel() { worker.terminate(); worker = null; finishRun(); setStatus('Cancelled. Nothing changed.'); }
function finishRun() { running = false; $('runBtn').textContent = 'Run simulation'; }
function onWorker(ev) {
  const m = ev.data; if (m.id !== jobId) return;
  if (m.type === 'progress') { $('bar').style.width = (100 * m.f).toFixed(1) + '%'; return; }
  finishRun();
  if (m.type === 'error') { $('bar').style.width = '0%'; setStatus(m.message, true); return; }
  $('bar').style.width = '100%';
  const R = m.result; R.form = pendingForm; R.label = describe(pendingForm);
  setStatus(`Done in ${((performance.now() - t0) / 1000).toFixed(1)} s, ${R.stats.nSteps.toLocaleString()} integration steps.`);
  showResult(R);
}
$('controls').addEventListener('submit', (e) => { e.preventDefault(); run(); });

// ------------------------------------------------------------------ results
let current = null;
const fmt = (x, d = 2) => (x >= 0 ? '+' : '−') + Math.abs(x).toFixed(d);
function showResult(R) {
  current = R;
  const S = R.S, p = R.params;
  $('runLabel').textContent = R.label;
  $('r-trev').innerHTML = S.reversal ? `${S.tCross.toFixed(2)} s` : S.spinDown ? `none <small>spin-down at ${S.tCross.toFixed(2)} s</small>` : `none <small>within ${p.tEnd} s</small>`;
  $('r-wend').innerHTML = `${fmt(S.w3End)} <small>rad/s</small>`;
  $('r-th').innerHTML = `${S.thetaMax.toFixed(2)}°`;
  $('r-de').innerHTML = `${Math.abs(100 * S.dEfrac) < 0.005 ? '0.00' : (-100 * S.dEfrac).toFixed(2)}%`;
  $('r-n').innerHTML = `${S.NminMg.toFixed(2)} to ${S.NmaxMg.toFixed(2)}`;
  $('r-slip').innerHTML = p.law === 'rolling' ? '0 <small>by definition</small>' : `${(1000 * S.slipMax).toFixed(1)} <small>mm/s</small>`;
  $('r-ratio').innerHTML = S.ratioMax.toFixed(2) + (p.law === 'rolling' ? ' <small>needed to roll</small>' : '');

  const v = $('verdict'); v.hidden = false;
  const pref = (p.w0 > 0) === (p.I12 < 0);
  if (S.reversal) v.textContent = `The spin reversed at ${S.tCross.toFixed(2)} s and ended at ${fmt(S.w3End)} rad/s` +
    (p.law === 'rolling' ? ', its full starting speed: rolling loses no energy.' : `, after losing ${(-100 * S.dEfrac).toFixed(1)}% of its energy.`);
  else if (S.spinDown) v.textContent = `The spin passed through zero at ${S.tCross.toFixed(2)} s but never picked up again: friction stopped the body. That is a spin-down, not a reversal.`;
  else if (pref) v.textContent = `Spun in its preferred sense, the body keeps turning the same way: ${S.thetaMax < 0.05 ? 'no rocking develops' : 'it rocks a little'} and the spin ends at ${fmt(S.w3End)} rad/s.`;
  else v.textContent = `No reversal within ${p.tEnd} s. The spin is ${fmt(S.w3End)} rad/s at the end` + (Math.abs(S.w3End) > 0.9 * Math.abs(p.w0) ? ': the contact is too slippery to couple the rocking to the spin.' : '.');
  if (S.NminMg <= 0) v.textContent += ' Warning: N became negative, so the body would leave the table and the model no longer applies.';

  const C = S.checks, rolling = p.law === 'rolling';
  const items = [
    [rolling || C.residual < 1e-5, 'Contact maintained', rolling ? 'exact for rolling' : `max |ġ| ${C.residual.toExponential(0)} m/s`],
    [C.energyOk, 'Energy never increases', ''],
    [C.Npositive, 'Normal reaction stays positive', `min ${S.NminMg.toFixed(2)} mg`],
    [rolling || C.denMin > 0, 'N formula well posed', rolling ? '' : `min denominator ${C.denMin.toFixed(2)}`],
  ];
  $('checks').innerHTML = items.map(([ok, t, val]) => `<li><span class="${ok ? 'ok' : 'no'}" aria-hidden="true">${ok ? '✓' : '✕'}</span><span>${t}<span class="sr-only">: ${ok ? 'passed' : 'failed'}</span></span><span class="val">${val}</span></li>`).join('');

  $('keepBtn').disabled = kept.length >= 6; $('csvBtn').disabled = false;
  $('keepBtn').textContent = kept.length >= 6 ? 'Six runs kept' : 'Keep for comparison';
  stage && stage.load(R);
  $('stageEmpty').hidden = true; $('hud').hidden = false;
  $('playBtn').disabled = false; $('scrub').disabled = false; $('scrub').max = R.T.length - 1;
  playhead = 0; setFrame(0);
  rebuildCharts();
  if (!reduceMotion && autoplayOnce) { autoplayOnce = false; play(); }
}
let autoplayOnce = true;

// CSV
$('csvBtn').addEventListener('click', () => {
  const R = current; if (!R) return;
  const head = 't_s,w1,w2,w3,phi_rad,theta_rad,psi_rad,v1,v2,v3,xG_m,yG_m,N_over_mg,F_N,slip_m_s,E_J,tilt_deg';
  const mg = R.params.m * R.params.g, lines = [`# ${R.label}`, head];
  for (let i = 0; i < R.T.length; i++) {
    const row = [R.T[i], ...R.Y.map(a => a[i]), R.D.N[i] / mg, R.D.F[i], R.D.slip[i], R.D.E[i], R.D.tilt[i]];
    lines.push(row.map(x => +x.toPrecision(8)).join(','));
  }
  const a = document.createElement('a');
  a.href = URL.createObjectURL(new Blob([lines.join('\n')], { type: 'text/csv' }));
  a.download = 'rattleback_' + R.label.replace(/[^a-z0-9]+/gi, '_').toLowerCase().slice(0, 60) + '.csv';
  a.click(); setTimeout(() => URL.revokeObjectURL(a.href), 2000);
});

// ------------------------------------------------------------------ comparison
const kept = [];
const SLOT = ['--s1', '--s2', '--s3', '--s4', '--s5', '--s6'];
$('keepBtn').addEventListener('click', () => {
  if (!current || kept.length >= 6 || kept.includes(current)) return;
  const used = new Set(kept.map(k => k.slot));
  current.slot = SLOT.find(s => !used.has(s));
  kept.push(current);
  $('keepBtn').disabled = true; $('keepBtn').textContent = 'Kept';
  rebuildCharts();
});
function renderCompareList() {
  $('compareEmpty').hidden = kept.length > 0;
  $('compareList').innerHTML = kept.map((k, i) =>
    `<li><span class="sw" style="background:var(${k.slot})"></span>${k.label}<button type="button" data-i="${i}" aria-label="Remove ${k.label}">remove</button></li>`).join('') +
    (current && !kept.includes(current) ? `<li><span class="sw" style="background:var(--ink)"></span>Current: ${current.label}</li>` : '');
}
$('compareList').addEventListener('click', (e) => {
  const i = e.target.dataset.i; if (i == null) return;
  kept.splice(+i, 1); if (current) { $('keepBtn').disabled = kept.includes(current) || kept.length >= 6; $('keepBtn').textContent = kept.includes(current) ? 'Kept' : 'Keep for comparison'; }
  rebuildCharts();
});

// ------------------------------------------------------------------ charts
let plots = [];
const sync = uPlot.sync('rb');
function baseOpts(el, series, yLabel, extra = {}) {
  const ink2 = css('--ink-2'), grid = css('--grid'), axis = css('--axis');
  const ax = { stroke: ink2, grid: { stroke: grid, width: 1 }, ticks: { stroke: axis, width: 1, size: 4 }, font: '12px Barlow, system-ui, sans-serif', labelFont: '13px Barlow, system-ui, sans-serif' };
  return {
    width: Math.max(280, el.clientWidth), height: 210, series: [{ label: 't (s)', value: (u, v) => v == null ? '–' : v.toFixed(3) + ' s' }, ...series],
    axes: [{ ...ax, label: 't (s)' }, { ...ax, label: yLabel, size: 52 }],
    scales: { x: { time: false } }, legend: { live: true },
    cursor: { sync: { key: sync.key }, drag: { x: true, y: false }, points: { size: 7 } },
    ...extra,
  };
}
function line(label, color, unit, extra = {}) {
  return { label, stroke: color, width: 1.6, points: { show: false }, value: (u, v) => v == null ? '–' : `${v.toFixed(3)} ${unit}`, ...extra };
}
function decimate(R, maxPts = 4000) { const n = R.T.length, step = Math.max(1, Math.ceil(n / maxPts)); const idx = []; for (let i = 0; i < n; i += step) idx.push(i); if (idx[idx.length - 1] !== n - 1) idx.push(n - 1); return idx; }
function pick(a, idx, f = (x) => x) { const o = new Float64Array(idx.length); for (let j = 0; j < idx.length; j++) o[j] = f(a[idx[j]]); return o; }

function rebuildCharts() {
  plots.forEach(p => p.u.destroy()); plots = [];
  renderCompareList();
  const R = current; if (!R) return;
  const idx = R.T.length > 4000 ? decimate(R, 4000) : [...Array(R.T.length).keys()];
  const t = pick(R.T, idx), mg = R.params.m * R.params.g;
  const ink = css('--ink'), s1 = css('--s1'), s2 = css('--s2'), s3 = css('--s3'), muted = css('--muted');
  const mk = (id, series, data, yl, extra) => {
    const el = $(id); const u = new uPlot(baseOpts(el, series, yl, extra), data, el);
    const mark = document.createElement('div');
    mark.style.cssText = `position:absolute;top:0;bottom:0;width:1px;background:${muted};pointer-events:none;left:0`;
    u.over.appendChild(mark); plots.push({ u, el, mark });
  };
  mk('p-w', [line('ω₁', s1, 'rad/s', { width: 1 }), line('ω₂', s2, 'rad/s', { width: 1 }), line('ω₃', s3, 'rad/s', { width: 2 })],
    [t, pick(R.Y[0], idx), pick(R.Y[1], idx), pick(R.Y[2], idx)], 'rad/s');
  mk('p-tilt', [line('tilt', s1, '°')], [t, pick(R.D.tilt, idx)], 'degrees');
  mk('p-n', [line('N/mg', s1, ''), { label: 'N = mg', stroke: muted, width: 1, dash: [5, 4], points: { show: false }, value: () => '1' }],
    [t, pick(R.D.N, idx, x => x / mg), t.map(() => 1)], 'N/mg');
  mk('p-f', [line('|F|', s1, 'N')], [t, pick(R.D.F, idx)], 'newtons');
  mk('p-slip', [line('slip speed', s1, 'mm/s')], [t, pick(R.D.slip, idx, x => 1000 * x)], 'mm/s');
  mk('p-e', [line('E', s1, 'J', { value: (u, v) => v == null ? '–' : v.toFixed(5) + ' J' })], [t, pick(R.D.E, idx)], 'joules');

  // comparison: kept runs in their own slots, current run in ink
  const runs = [...kept]; if (!kept.includes(R)) runs.push(R);
  const longest = runs.reduce((a, b) => (b.T.length > a.T.length ? b : a));
  const lidx = decimate(longest, 3000), tl = pick(longest.T, lidx), dt = longest.params.dtOut;
  const cdata = [tl], cser = [];
  for (const k of runs) {
    const col = k === R && !kept.includes(R) ? ink : css(k.slot);
    cdata.push(Array.from(tl, (tt) => { const i = Math.round(tt / dt); return i < k.T.length ? k.Y[2][i] : null; }));
    cser.push(line(k.label.length > 40 ? k.label.slice(0, 38) + '…' : k.label, col, 'rad/s', { width: k === R ? 2 : 1.6, spanGaps: false }));
  }
  mk('p-compare', cser, cdata, 'ω₃ (rad/s)', { legend: { show: false } });
  setMarks();
}
new ResizeObserver(() => plots.forEach(p => p.u.setSize({ width: Math.max(280, p.el.clientWidth), height: 210 }))).observe($('charts'));
function setMarks() {
  if (!current) return; const tt = current.T[playhead] ?? 0;
  for (const p of plots) { const x = p.u.valToPos(tt, 'x'); p.mark.style.transform = `translateX(${x}px)`; p.mark.style.display = Number.isFinite(x) ? '' : 'none'; }
}

// ------------------------------------------------------------------ playback
let stage = null, playhead = 0, playing = false, simT = 0, lastTs = 0;
function setFrame(i) {
  const R = current; if (!R) return;
  playhead = Math.max(0, Math.min(R.T.length - 1, i));
  stage && stage.show(playhead);
  const w3 = R.Y[2][playhead];
  $('hudW').textContent = `ω₃ ${fmt(w3)} rad/s`;
  $('hudT').textContent = `t = ${R.T[playhead].toFixed(2)} s, tilt ${R.D.tilt[playhead].toFixed(2)}°`;
  const ex = stage ? stage.exag : 1;
  $('hudX').textContent = ex !== 1 ? `Tilt drawn ×${ex}` : '';
  $('clock').textContent = `${R.T[playhead].toFixed(2)} s`;
  $('scrub').value = playhead;
  setMarks();
}
function play() {
  if (!current) return;
  if (playhead >= current.T.length - 1) playhead = 0;
  playing = true; simT = current.T[playhead]; lastTs = 0; $('playBtn').textContent = 'Pause';
  requestAnimationFrame(tick);
}
function pause() { playing = false; $('playBtn').textContent = 'Play'; }
function tick(ts) {
  if (!playing || !current) return;
  if (lastTs) simT += (ts - lastTs) / 1000 * parseFloat($('speed').value);
  lastTs = ts;
  const i = Math.round(simT / current.params.dtOut);
  if (i >= current.T.length - 1) { setFrame(current.T.length - 1); pause(); return; }
  setFrame(i); requestAnimationFrame(tick);
}
$('playBtn').addEventListener('click', () => (playing ? pause() : play()));
$('scrub').addEventListener('input', () => { pause(); setFrame(parseInt($('scrub').value, 10)); });
$('exag').addEventListener('change', () => { if (stage) { stage.exag = parseFloat($('exag').value); setFrame(playhead); } });
$('trail').addEventListener('change', () => stage && stage.setTrail($('trail').checked));
$('follow').addEventListener('change', () => { if (stage) stage.follow = $('follow').checked; });

// ------------------------------------------------------------------ start
(async function start() {
  const fromLink = fromHash();
  if (fromLink) { writeForm(fromLink); presetSel.value = 'custom'; $('presetNote').textContent = 'Settings from a shared link.'; }
  else applyPreset('coulomb');
  try {
    const { Stage } = await import('./scene.js');
    stage = new Stage($('stage'));
    $('stageEmpty').textContent = 'Press Run simulation to see the motion.';
  } catch (e) {
    $('stageEmpty').textContent = 'The 3D view needs WebGL, which this browser has turned off. The results and charts still work.';
  }
  run();
})();
