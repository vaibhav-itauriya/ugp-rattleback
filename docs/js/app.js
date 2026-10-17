import uPlot from '../vendor/uPlot.esm.js';

const $ = (id) => document.getElementById(id);
const css = (name) => getComputedStyle(document.documentElement).getPropertyValue(name).trim();
const reduceMotion = matchMedia('(prefers-reduced-motion: reduce)').matches;
const sgn = (x, d = 2) => (x >= 0 ? '+' : '−') + Math.abs(x).toFixed(d);
const fmtInt = (n) => String(Math.round(n)).replace(/\B(?=(\d{3})+(?!\d))/g, ',');
const esc = (s) => String(s).replace(/[&<>"]/g, (c) => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;' }[c]));

// =====================================================================  setup fields
const EPS_OPTS = [['1e-3', '1×10⁻³'], ['3e-4', '3×10⁻⁴'], ['1e-4', '1×10⁻⁴'], ['3e-5', '3×10⁻⁵'], ['1e-5', '1×10⁻⁵']];
const FIELDS = {
  contact: [
    { id: 'mu', label: 'Friction coefficient μ', min: 0.001, max: 5, step: 0.01, slider: [0.01, 2, 0.01], show: f => f.law === 'coulomb' },
    { id: 'eps', label: 'Slip regularization ε', unit: 'm/s', select: EPS_OPTS, show: f => f.law === 'coulomb',
      hint: 'Smooths the Coulomb law at zero slip. Converged at 3×10⁻⁵; 10⁻³ reverses about 15% late.' },
    { id: 'drill', label: 'Drilling (spin) friction at the contact patch', check: true, show: f => f.law === 'coulomb' },
    { id: 'rho', label: 'Drilling lever arm ρ', unit: 'mm', min: 0, max: 5, step: 0.05, slider: [0, 2.5, 0.05], show: f => f.law === 'coulomb' && f.drill,
      hint: 'Hertzian estimate 0.16–0.52 mm. Near 2 mm the body stops before it can rock.' },
    { id: 'cv', label: 'Viscous coefficient c', unit: 'kg/s', min: 0.01, max: 1e5, step: 'any', slider: [0, 3.5, 0.01], log: true, show: f => f.law === 'viscous',
      hint: 'Slider is logarithmic, 1 to 3000 kg/s. Below about 200 there is no reversal.' },
    { id: 'nmode', label: 'Normal reaction N', select: [['solved', 'Solved'], ['fixed', 'Fixed at mg']], show: f => f.law !== 'rolling',
      hint: 'Fixed at mg is the common shortcut; it is 14% off near the rolling limit.' },
  ],
  launch: [
    { id: 'w0', label: 'Initial spin |ω₀|', unit: 'rad/s', min: 0.1, max: 40, step: 0.5, slider: [0.5, 15, 0.5] },
    { id: 'tEnd', label: 'Simulated time', unit: 's', min: 0.5, max: 30, step: 1 },
    { id: 'tilt0', label: 'Initial tilt', unit: 'deg', min: 0, max: 15, step: 0.1, hint: '0 starts from the paper\'s 10⁻⁶ rad nudge. A larger tilt is like a tap on one end.' },
  ],
  body: [
    { id: 'm', label: 'Mass m', unit: 'kg', min: 0.001, max: 100, step: 0.1 },
    { id: 'a1', label: 'Semi-axis a₁', unit: 'm', min: 0.001, max: 2, step: 0.01 },
    { id: 'a2', label: 'Semi-axis a₂', unit: 'm', min: 0.001, max: 2, step: 0.005 },
    { id: 'a3', label: 'Semi-axis a₃', unit: 'm', min: 0.001, max: 2, step: 0.005 },
    { id: 'b', label: 'Offset of centre above G, b', unit: 'm', min: 0, max: 1, step: 0.001 },
    { id: 'I11', label: 'I₁₁', unit: '×10⁻⁴ kg m²', min: 1e-6, max: 1e6, step: 0.1 },
    { id: 'I22', label: 'I₂₂', unit: '×10⁻⁴ kg m²', min: 1e-6, max: 1e6, step: 0.1 },
    { id: 'I33', label: 'I₃₃', unit: '×10⁻⁴ kg m²', min: 1e-6, max: 1e6, step: 0.1 },
    { id: 'I12', label: 'I₁₂', unit: '×10⁻⁴ kg m²', min: -1e6, max: 1e6, step: 0.1, hint: 'The chirality: its sign picks the preferred sense; at 0 nothing reverses.' },
  ],
  solver: [
    { id: 'tol', label: 'Tolerance', select: [['std', 'Standard, 10⁻⁸'], ['fast', 'Fast, 10⁻⁶'], ['strict', 'Strict, 10⁻¹⁰']] },
    { id: 'dtOut', label: 'Output step', select: [['0.001', '1 ms'], ['0.002', '2 ms'], ['0.005', '5 ms']], hint: 'Spacing of the stored samples. The solver picks its own internal steps.' },
  ],
};
const ALL = Object.values(FIELDS).flat();
const NUMERIC = ALL.filter(f => !f.select && !f.check).map(f => f.id);
const BASE = { law: 'coulomb', sense: '-1', mu: 0.15, eps: '3e-5', drill: false, rho: 0.5, cv: 500, nmode: 'solved',
  w0: 5, tEnd: 8, tilt0: 0, m: 1, a1: 0.2, a2: 0.03, a3: 0.02, b: 0.005, I11: 2, I22: 16, I33: 17.5, I12: -0.5, tol: 'std', dtOut: '0.001' };

function fieldHTML(f) {
  if (f.check) return `<label class="chk" data-fid="${f.id}"><input type="checkbox" id="${f.id}"> ${f.label}</label>`;
  const unit = f.unit ? ` <span class="u">${f.unit}</span>` : '';
  const ctl = f.select
    ? `<select id="${f.id}">${f.select.map(([v, t]) => `<option value="${v}">${t}</option>`).join('')}</select>`
    : `<input type="number" id="${f.id}" step="${f.step}" min="${f.min}" max="${f.max}" inputmode="decimal">`;
  const sl = f.slider ? `<input type="range" id="${f.id}-r" min="${f.slider[0]}" max="${f.slider[1]}" step="${f.slider[2]}" aria-label="${f.label} slider" tabindex="-1">` : '';
  const hint = f.hint ? `<span class="hint">${f.hint}</span>` : '';
  return `<div class="f" data-fid="${f.id}"><label for="${f.id}">${f.label}${unit}</label>${ctl}${sl}<span class="err" id="${f.id}-err" hidden></span>${hint}</div>`;
}
for (const [g, list] of Object.entries(FIELDS)) $('fields-' + g).innerHTML = list.map(fieldHTML).join('');

// =====================================================================  presets
const PRESETS = [
  { id: 'coulomb', name: 'Coulomb friction, μ = 0.15', p: {}, note: 'The paper\'s main case: reverses at 4.22 s and loses 12% of its energy.' },
  { id: 'rolling', name: 'Rolling without slipping', p: { law: 'rolling' }, note: 'Conservative: the spin returns at full speed after reversing at 2.83 s.' },
  { id: 'drill', name: 'Coulomb + drilling, ρ = 0.5 mm', p: { drill: true, rho: 0.5 }, note: 'Reverses at 3.98 s, then the drilling torque brings it nearly to rest.' },
  { id: 'viscous', name: 'Viscous friction, c = 500 kg/s', p: { law: 'viscous', cv: 500 }, note: 'Reverses at 2.90 s, close to rolling.' },
  { id: 'weakvisc', name: 'Weak viscous friction, c = 1.5 kg/s', p: { law: 'viscous', cv: 1.5 }, note: 'Same nominal size as μ = 0.15, yet the spin never reverses.' },
  { id: 'bigpatch', name: 'Oversized contact patch, ρ = 2 mm', p: { drill: true, rho: 2 }, note: 'Drilling stops the body at about 3 s: a spin-down, not a reversal.' },
  { id: 'fixn', name: 'Shortcut N = mg, μ = 1', p: { mu: 1, nmode: 'fixed' }, note: 'Reverses at 2.43 s instead of the correct 2.83 s.' },
  { id: 'coarse', name: 'Unconverged ε = 10⁻³ m/s', p: { eps: '1e-3' }, note: 'Reverses at 4.85 s, about 15% later than the converged 4.22 s.' },
  { id: 'pref', name: 'Preferred sense, Coulomb', p: { sense: '1' }, note: 'Spun the other way: no rocking, no slip, no reversal.' },
];
const presetSel = $('preset');
presetSel.innerHTML = PRESETS.map(p => `<option value="${p.id}">${p.name}</option>`).join('') + '<option value="custom">Custom setup</option>';

// =====================================================================  form state
function readForm() {
  const f = { law: document.querySelector('input[name=law]:checked').value, sense: document.querySelector('input[name=sense]:checked').value };
  for (const d of ALL) f[d.id] = d.check ? $(d.id).checked : d.select ? $(d.id).value : parseFloat($(d.id).value);
  return f;
}
function writeForm(f) {
  document.querySelector(`input[name=law][value=${f.law}]`).checked = true;
  document.querySelector(`input[name=sense][value="${f.sense}"]`).checked = true;
  for (const d of ALL) { if (d.check) $(d.id).checked = !!f[d.id]; else $(d.id).value = f[d.id]; }
  syncSliders(); refreshSetup();
}
function validate(f) {
  const e = {};
  for (const d of ALL) {
    if (d.select || d.check) continue;
    const v = f[d.id];
    if (!Number.isFinite(v)) e[d.id] = 'Enter a number.';
    else if (v < d.min || v > d.max) e[d.id] = `Use a value from ${d.min} to ${d.max}.`;
  }
  if (!e.b && !e.a3 && f.b >= f.a3) e.b = 'Must be smaller than a₃, or the body falls over.';
  if (!e.I12 && !e.I11 && !e.I22 && f.I11 * f.I22 <= f.I12 * f.I12) e.I12 = 'Too large: needs I₁₂² < I₁₁ I₂₂.';
  if (f.law === 'coulomb' && !e.mu && f.mu <= 0) e.mu = 'Must be above zero.';
  return e;
}
const shown = (d, f) => !d.show || d.show(f);
function showErrors(errs, f) {
  for (const d of ALL) {
    if (d.select || d.check) continue;
    const box = document.querySelector(`.f[data-fid="${d.id}"]`), msg = $(d.id + '-err');
    const bad = errs[d.id] && shown(d, f);
    box.classList.toggle('invalid', !!bad); msg.hidden = !bad; msg.textContent = bad ? errs[d.id] : '';
    $(d.id).setAttribute('aria-invalid', bad ? 'true' : 'false');
  }
}
const prefSign = (f) => (f.I12 > 0 ? -1 : 1);       // +w3 is preferred when I12 < 0
function toParams(f) {
  const tilt = f.tilt0 > 0 ? f.tilt0 * Math.PI / 180 / Math.SQRT2 : 1e-6;
  const tol = { std: [1e-8, 1e-10], fast: [1e-6, 1e-8], strict: [1e-10, 1e-12] }[f.tol];
  return {
    law: f.law, mu: f.mu, eps: parseFloat(f.eps), rho: f.law === 'coulomb' && f.drill ? f.rho * 1e-3 : 0, cv: f.cv,
    fixN: f.law !== 'rolling' && f.nmode === 'fixed',
    w0: Math.abs(f.w0) * (f.sense === '1' ? prefSign(f) : -prefSign(f)), tEnd: f.tEnd,
    phi0: tilt, theta0: tilt, psi0: 1e-6,
    m: f.m, a1: f.a1, a2: f.a2, a3: f.a3, b: f.b,
    I11: f.I11 * 1e-4, I22: f.I22 * 1e-4, I33: f.I33 * 1e-4, I12: f.I12 * 1e-4,
    rtol: tol[0], atol: tol[1], dtOut: parseFloat(f.dtOut),
  };
}
function describe(f) {
  const s = f.sense === '1' ? 'preferred' : 'unstable';
  const nm = f.law !== 'rolling' && f.nmode === 'fixed' ? ', N = mg' : '';
  const w = `, ${s} sense, |ω₀| ${f.w0}`;
  if (f.law === 'rolling') return `Rolling${w}`;
  if (f.law === 'viscous') return `Viscous c ${f.cv}${nm}${w}`;
  return `Coulomb μ ${f.mu}${f.drill ? `, ρ ${f.rho} mm` : ''}${f.eps !== '3e-5' ? `, ε ${f.eps}` : ''}${nm}${w}`;
}
function refreshSetup() {
  const f = readForm();
  document.querySelectorAll('[data-fid]').forEach(el => { const d = ALL.find(x => x.id === el.dataset.fid); el.hidden = !shown(d, f); });
  showErrors(validate(f), f);
  $('sum-contact').textContent = f.law === 'rolling' ? 'no slip' : f.law === 'viscous' ? `c = ${f.cv} kg/s` : `μ = ${f.mu}${f.drill ? `, ρ = ${f.rho} mm` : ''}`;
  $('sum-launch').textContent = `${f.sense === '1' ? '+' : '−'}${f.w0} rad/s, ${f.tEnd} s`;
  $('sum-solver').textContent = { std: 'standard', fast: 'fast', strict: 'strict' }[f.tol] + `, ${parseFloat(f.dtOut) * 1000} ms`;
  const pref = prefSign(f);
  $('derived').innerHTML = f.I12 === 0 ? 'With I₁₂ = 0 the body has no preferred sense.'
    : `Preferred sense: <b>${pref > 0 ? 'positive' : 'negative'} ω₃</b> (${pref > 0 ? 'anticlockwise' : 'clockwise'} seen from above). Resting height of G: <b>${((f.a3 - f.b) * 1000).toFixed(1)} mm</b>.`;
  $('dirty').hidden = !lastRunForm || JSON.stringify(f) === JSON.stringify(lastRunForm);
}
// paired number and slider inputs
function syncSliders() {
  for (const d of ALL) if (d.slider) { const v = parseFloat($(d.id).value); if (Number.isFinite(v)) $(d.id + '-r').value = d.log ? Math.log10(Math.max(v, 1)) : v; }
}
for (const d of ALL) if (d.slider) {
  $(d.id + '-r').addEventListener('input', () => { const v = parseFloat($(d.id + '-r').value); $(d.id).value = d.log ? +Math.pow(10, v).toPrecision(2) : v; markCustom(); refreshSetup(); });
  $(d.id).addEventListener('input', syncSliders);
}
$('setup').addEventListener('input', () => { markCustom(); refreshSetup(); });
$('setup').addEventListener('submit', (e) => { e.preventDefault(); run(); });
function applyPreset(id, quiet) {
  const P = PRESETS.find(p => p.id === id); if (!P) return;
  writeForm({ ...BASE, ...P.p }); presetSel.value = id; if (!quiet) log(`Case loaded: ${P.name}. ${P.note}`);
}
function markCustom() { presetSel.value = 'custom'; }
presetSel.addEventListener('change', () => { if (presetSel.value !== 'custom') { applyPreset(presetSel.value); run(); } });

// =====================================================================  log, status, toast
function log(msg, isErr = false) {
  const li = document.createElement('li');
  const d = new Date(), p2 = (x) => String(x).padStart(2, '0');
  li.innerHTML = `<time>${p2(d.getHours())}:${p2(d.getMinutes())}:${p2(d.getSeconds())}</time><span class="${isErr ? 'e' : ''}">${esc(msg)}</span>`;
  $('log').prepend(li);
}
function status(state, msg = '', kind = 'ok') {
  $('sState').textContent = state; $('sMsg').textContent = msg;
  $('sDot').className = 'dot' + (kind === 'busy' ? ' busy' : kind === 'err' ? ' err' : '');
}
let toastTimer;
function toast(msg) {
  let t = document.querySelector('.toast'); if (!t) { t = document.createElement('div'); t.className = 'toast'; t.setAttribute('role', 'status'); document.body.append(t); }
  t.textContent = msg; t.hidden = false; clearTimeout(toastTimer); toastTimer = setTimeout(() => (t.hidden = true), 2400);
}

// =====================================================================  worker
let worker = null, job = 0, busy = null, lastRunForm = null, runStart = 0;
function getWorker() {
  if (worker) return worker;
  worker = new Worker(new URL('./worker.js', import.meta.url), { type: 'module' });
  worker.onmessage = onWorker;
  worker.onerror = (e) => { setBusy(null); status('Error', 'The solver could not start in this browser.', 'err'); log('Worker error: ' + (e.message || 'unknown'), true); };
  return worker;
}
function setBusy(kind) {
  busy = kind;
  $('stopBtn').disabled = !kind; $('runBtn').disabled = false; $('swRun').disabled = kind === 'sweep';
  $('swRun').textContent = kind === 'sweep' ? 'Sweeping…' : 'Run sweep';
  if (!kind) $('sBar').style.width = '0%';
}
function stop() {
  if (!busy) return;
  worker.terminate(); worker = null;
  log(busy === 'sweep' ? 'Sweep stopped. Rows computed so far are kept.' : 'Run stopped. The previous result is still shown.');
  status('Stopped', ''); setBusy(null);
}
function run() {
  const f = readForm(), errs = validate(f);
  showErrors(errs, f);
  const bad = Object.keys(errs).find(k => shown(ALL.find(d => d.id === k), f));
  if (bad) {
    const d = ALL.find(x => x.id === bad);
    const sec = $(bad).closest('details'); if (sec) sec.open = true;
    $(bad).focus(); status('Check the setup', `${d.label}: ${errs[bad]}`, 'err'); return;
  }
  if (busy) stop();
  const id = ++job; runStart = performance.now(); setBusy('run');
  status('Computing', describe(f), 'busy');
  getWorker().postMessage({ type: 'run', id, params: toParams(f) });
  pendingForm = f;
}
let pendingForm = null;
function onWorker(ev) {
  const m = ev.data; if (m.id !== job) return;
  if (m.type === 'progress') { $('sBar').style.width = (100 * m.f).toFixed(1) + '%'; if (busy === 'run') $('sState').textContent = `Computing ${Math.round(100 * m.f)}%`; else $('sState').textContent = `Sweeping ${Math.round(100 * m.f)}%`; return; }
  if (m.type === 'error') { setBusy(null); status('Error', m.message, 'err'); log(m.message, true); return; }
  if (m.type === 'done') {
    setBusy(null);
    const R = m.result, secs = (performance.now() - runStart) / 1000;
    R.form = pendingForm; R.label = describe(pendingForm); R.id = ++runCounter; R.secs = secs;
    lastRunForm = pendingForm;
    status('Done', R.label);
    $('sStats').textContent = `${fmtInt(R.stats.nSteps)} steps, ${secs.toFixed(2)} s`;
    log(`Run ${R.id}: ${R.label}. ${R.S.reversal ? `Reversal at ${R.S.tCross.toFixed(3)} s` : 'No reversal'}, ω₃ ends at ${sgn(R.S.w3End, 3)} rad/s. ${fmtInt(R.stats.nSteps)} steps in ${secs.toFixed(2)} s.`);
    addHistory(R); show(R); return;
  }
  if (m.type === 'row') { sweep.rows[m.i] = m.row; renderSweep(); if (!m.row.ok) log(`Sweep point ${fmtX(m.row.x)}: ${m.row.error}`, true); return; }
  if (m.type === 'sweepDone') { setBusy(null); status('Done', `Sweep of ${sweep.def.label} finished`); $('sStats').textContent = `${sweep.rows.length} runs, ${(m.ms / 1000).toFixed(1)} s`; log(`Sweep of ${sweep.def.label}: ${sweep.rows.length} runs in ${(m.ms / 1000).toFixed(1)} s.`); $('swCsv').disabled = false; }
}
$('runBtn').addEventListener('click', run);
$('stopBtn').addEventListener('click', stop);

// =====================================================================  results
let current = null, runCounter = 0;
function show(R) {
  current = R; const S = R.S, p = R.params;
  $('resBadge').textContent = `Run ${R.id}`;
  $('runLabel').textContent = R.label;
  $('r-trev').innerHTML = S.reversal ? `${S.tCross.toFixed(2)}<small>s</small>` : S.spinDown ? `none<small>spin-down</small>` : `none<small>in ${p.tEnd} s</small>`;
  $('r-wend').innerHTML = `${sgn(S.w3End)}<small>rad/s</small>`;
  $('r-th').innerHTML = `${S.thetaMax.toFixed(2)}°`;
  $('r-de').innerHTML = `${Math.abs(S.dEfrac) < 5e-5 ? '0.00' : (-100 * S.dEfrac).toFixed(2)}%`;
  $('r-n').innerHTML = `${S.NminMg.toFixed(2)} to ${S.NmaxMg.toFixed(2)}`;
  $('r-slip').innerHTML = p.law === 'rolling' ? '0<small>no slip</small>' : `${(1000 * S.slipMax).toFixed(1)}<small>mm/s</small>`;
  $('r-ratio').innerHTML = `${S.ratioMax.toFixed(2)}${p.law === 'rolling' ? '<small>needed</small>' : ''}`;
  const v = $('verdict'); v.hidden = false;
  const pref = Math.sign(p.w0) === (p.I12 > 0 ? -1 : 1);
  if (S.reversal) v.textContent = `Reversed at ${S.tCross.toFixed(2)} s and ended at ${sgn(S.w3End)} rad/s` + (p.law === 'rolling' ? ', its full starting speed.' : `, after losing ${(-100 * S.dEfrac).toFixed(1)}% of its energy.`);
  else if (S.spinDown) v.textContent = `Passed through zero spin at ${S.tCross.toFixed(2)} s without picking up again: friction stopped the body. A spin-down, not a reversal.`;
  else if (pref) v.textContent = `Preferred sense: it keeps turning the same way, ${S.thetaMax < 0.05 ? 'with no rocking' : 'rocking slightly'}, and ends at ${sgn(S.w3End)} rad/s.`;
  else v.textContent = `No reversal within ${p.tEnd} s; the spin is ${sgn(S.w3End)} rad/s at the end` + (Math.abs(S.w3End) > 0.9 * Math.abs(p.w0) ? '. The contact is too slippery to couple the rocking to the spin.' : '.');
  if (S.NminMg <= 0) v.textContent += ' N became negative: the body would leave the table, so the result is not valid.';
  const C = S.checks, rolling = p.law === 'rolling';
  const items = [
    [rolling || C.residual < 1e-5, 'Contact maintained', rolling ? 'exact for rolling' : `largest gap rate ${C.residual.toExponential(0)} m/s`],
    [C.energyOk, 'Energy never increases', `${S.E0.toFixed(4)} J to ${S.Eend.toFixed(4)} J`],
    [C.Npositive, 'Body stays on the table', `N ≥ ${S.NminMg.toFixed(2)} mg`],
    [rolling || C.denMin > 0, 'N formula well posed', rolling ? 'not needed for rolling' : `denominator ≥ ${C.denMin.toFixed(3)}`],
  ];
  $('checks').innerHTML = items.map(([ok, t, val]) => `<li><span class="${ok ? 'ok' : 'no'}" aria-hidden="true">${ok ? '✓' : '✕'}</span><span>${t}<span class="sr-only">, ${ok ? 'passed' : 'failed'}</span></span><span class="val">${val}</span></li>`).join('');
  for (const b of ['bStart', 'bBack', 'bPlay', 'bFwd', 'bEnd', 'scrub']) $(b).disabled = false;
  $('bRev').disabled = !S.reversal;
  $('scrub').max = R.T.length - 1;
  const mk = $('revMark');
  if (S.tCross != null) { mk.hidden = false; mk.style.left = `calc(${(100 * S.tCross / p.tEnd).toFixed(3)}% - 1px)`; mk.dataset.label = S.reversal ? 'reversal' : 'zero spin'; } else mk.hidden = true;
  document.querySelectorAll('#exportMenu [data-act]').forEach(b => (b.disabled = false));
  if (stage) { stage.load(R); $('vpEmpty').hidden = true; }
  $('hud').hidden = false;
  pause(); setFrame(0);
  dirtyCharts(); renderHistory(); refreshSetup();
  if (firstAuto && !reduceMotion) { firstAuto = false; play(); }
}
let firstAuto = true;

// =====================================================================  history & compare
const runHist = [];             // newest first, max 12
const SLOTS = ['--s1', '--s2', '--s3', '--s4', '--s5', '--s6'];
function addHistory(R) {
  runHist.unshift(R); R.shown = false;
  while (runHist.length > 12) runHist.pop();
  $('histCount').textContent = runHist.length;
}
function renderHistory() {
  const tb = $('histTbl').tBodies[0];
  tb.innerHTML = runHist.map((R, i) => {
    const S = R.S, sw = R.slot ? `<span class="sw" style="background:var(${R.slot})"></span> ` : '';
    return `<tr class="${R === current ? 'current' : ''}">
      <td><input type="checkbox" data-i="${i}" ${R.shown ? 'checked' : ''} aria-label="Overlay run ${R.id}"></td>
      <td class="label" title="${esc(R.label)}">${sw}<b>${R.id}</b> ${esc(R.label)}</td>
      <td class="n">${S.reversal ? S.tCross.toFixed(2) + ' s' : '–'}</td>
      <td class="n">${sgn(S.w3End)}</td><td class="n">${S.thetaMax.toFixed(2)}°</td>
      <td class="n">${(-100 * S.dEfrac).toFixed(1)}%</td>
      <td><button class="link" type="button" data-open="${i}">Show</button><button class="link" type="button" data-load="${i}">Use setup</button></td></tr>`;
  }).join('');
}
$('histTbl').addEventListener('change', (e) => {
  const i = e.target.dataset.i; if (i == null) return;
  const R = runHist[+i];
  if (e.target.checked) {
    const used = new Set(runHist.filter(h => h.shown).map(h => h.slot));
    const free = SLOTS.find(s => !used.has(s));
    if (!free) { e.target.checked = false; toast('Six runs are already overlaid. Untick one first.'); return; }
    R.shown = true; R.slot = free;
  } else { R.shown = false; R.slot = null; }
  renderHistory(); buildCompare();
});
$('histTbl').addEventListener('click', (e) => {
  const o = e.target.dataset.open, l = e.target.dataset.load;
  if (o != null) show(runHist[+o]);
  if (l != null) { writeForm(runHist[+l].form); markCustom(); refreshSetup(); toast(`Setup of run ${runHist[+l].id} loaded. Press Run to recompute.`); }
});

// =====================================================================  charts
let plots = {};
const syncKey = uPlot.sync('rb').key;
function axisStyle() {
  return { stroke: css('--ink-2'), grid: { stroke: css('--grid'), width: 1 }, ticks: { stroke: css('--axis'), width: 1, size: 4 },
    font: '12px Barlow, system-ui, sans-serif', labelFont: '12px Barlow, system-ui, sans-serif' };
}
function lineS(label, color, unit, extra = {}) {
  return { label, stroke: color, width: 1.5, points: { show: false }, value: (u, v) => (v == null ? '–' : `${v.toFixed(3)} ${unit}`), ...extra };
}
function makePlot(key, el, series, data, yLabel, extra = {}) {
  if (plots[key]) plots[key].u.destroy();
  const ax = axisStyle();
  const opts = {
    width: Math.max(260, el.clientWidth), height: extra.height || 190,
    series: [{ label: extra.xLabel || 't', value: (u, v) => (v == null ? '–' : extra.xFmt ? extra.xFmt(v) : v.toFixed(3) + ' s') }, ...series],
    axes: [{ ...ax, label: extra.xLabel || 't (s)' }, { ...ax, label: yLabel, size: 50 }],
    scales: { x: { time: false, ...(extra.xLog ? { distr: 3 } : {}) } },
    cursor: { sync: extra.noSync ? undefined : { key: syncKey }, drag: { x: true, y: false }, points: { size: 6 } },
    legend: { show: extra.legend !== false },
  };
  const u = new uPlot(opts, data, el);
  const mark = document.createElement('div');
  mark.style.cssText = `position:absolute;top:0;bottom:0;width:1px;background:${css('--muted')};pointer-events:none;left:0`;
  if (!extra.noMark) u.over.appendChild(mark);
  if (!extra.noSeek) {
    let down = null;
    u.over.addEventListener('mousedown', (e) => (down = e.clientX));
    u.over.addEventListener('click', (e) => { if (down != null && Math.abs(e.clientX - down) > 3) return; seekTime(u.posToVal(e.offsetX, 'x')); });
  }
  plots[key] = { u, el, mark: extra.noMark ? null : mark };
  return u;
}
function decimate(n, maxPts) { const step = Math.max(1, Math.ceil(n / maxPts)), idx = []; for (let i = 0; i < n; i += step) idx.push(i); if (idx[idx.length - 1] !== n - 1) idx.push(n - 1); return idx; }
const pick = (a, idx, f = (x) => x) => { const o = new Float64Array(idx.length); for (let j = 0; j < idx.length; j++) o[j] = f(a[idx[j]]); return o; };
let chartsStale = { hist: true, cmp: true };
function dirtyCharts() { chartsStale.hist = chartsStale.cmp = true; buildVisible(); }
function buildHistories() {
  const R = current; if (!R) return; chartsStale.hist = false;
  const idx = decimate(R.T.length, 4000), t = pick(R.T, idx), mg = R.params.m * R.params.g;
  const s1 = css('--s1'), s2 = css('--s2'), s3 = css('--s3'), muted = css('--muted');
  makePlot('w', $('p-w'), [lineS('ω₁', s1, 'rad/s', { width: 1 }), lineS('ω₂', s2, 'rad/s', { width: 1 }), lineS('ω₃', s3, 'rad/s', { width: 2 })],
    [t, pick(R.Y[0], idx), pick(R.Y[1], idx), pick(R.Y[2], idx)], 'rad/s');
  makePlot('tilt', $('p-tilt'), [lineS('tilt', s1, '°')], [t, pick(R.D.tilt, idx)], 'degrees');
  makePlot('n', $('p-n'), [lineS('N/mg', s1, ''), { label: 'N = mg', stroke: muted, width: 1, dash: [5, 4], points: { show: false }, value: () => '1' }],
    [t, pick(R.D.N, idx, x => x / mg), t.map(() => 1)], 'N/mg');
  makePlot('f', $('p-f'), [lineS('|F|', s1, 'N')], [t, pick(R.D.F, idx)], 'newtons');
  makePlot('slip', $('p-slip'), [lineS('slip', s1, 'mm/s')], [t, pick(R.D.slip, idx, x => 1000 * x)], 'mm/s');
  makePlot('e', $('p-e'), [lineS('E', s1, 'J', { value: (u, v) => (v == null ? '–' : v.toFixed(5) + ' J') })], [t, pick(R.D.E, idx)], 'joules');
  setMarks();
}
function buildCompare() {
  chartsStale.cmp = false; const R = current; if (!R) return;
  const runs = runHist.filter(h => h.shown); if (!runs.includes(R)) runs.push(R);
  const longest = runs.reduce((a, b) => (b.params.tEnd > a.params.tEnd ? b : a));
  const tl = Array.from({ length: 2001 }, (_, i) => (longest.params.tEnd * i) / 2000);
  const data = [tl], ser = [];
  for (const k of runs) {
    const dt = k.params.dtOut;
    data.push(tl.map(tt => { const i = Math.round(tt / dt); return i < k.T.length ? k.Y[2][i] : null; }));
    ser.push(lineS(`${k.id}`, k.shown ? css(k.slot) : css('--ink'), 'rad/s', { width: k === R ? 2 : 1.5 }));
  }
  makePlot('cmp', $('p-cmp'), ser, data, 'ω₃ (rad/s)', { height: 230 });
  setMarks();
}
function buildVisible() {
  const tab = document.querySelector('.tabs [aria-selected=true]').id;
  if (tab === 't-hist' && chartsStale.hist) buildHistories();
  if (tab === 't-cmp') { renderHistory(); if (chartsStale.cmp) buildCompare(); }
  if (tab === 't-sweep') renderSweep();
}
function setMarks() {
  if (!current) return; const tt = current.T[playhead] ?? 0;
  for (const p of Object.values(plots)) if (p.mark) { const x = p.u.valToPos(tt, 'x'); p.mark.style.transform = `translateX(${x}px)`; p.mark.style.display = Number.isFinite(x) && x >= 0 ? '' : 'none'; }
}
new ResizeObserver(() => { for (const p of Object.values(plots)) if (p.el.clientWidth) p.u.setSize({ width: Math.max(260, p.el.clientWidth), height: p.u.height }); setMarks(); }).observe(document.querySelector('.dock'));

// =====================================================================  playback
let stage = null, playhead = 0, playing = false, simT = 0, lastTs = 0;
function setFrame(i) {
  const R = current; if (!R) return;
  playhead = Math.max(0, Math.min(R.T.length - 1, Math.round(i)));
  const k = playhead, t = R.T[k], mg = R.params.m * R.params.g;
  if (stage) stage.show(k);
  $('hudW').textContent = `ω₃ ${sgn(R.Y[2][k])} rad/s`;
  $('hudT').textContent = `t = ${t.toFixed(2)} s`;
  $('hudTilt').textContent = `tilt ${R.D.tilt[k].toFixed(2)}°`;
  $('hudX').textContent = stage && stage.exag !== 1 ? `Tilt drawn ×${stage.exag}` : '';
  $('clock').textContent = `${t.toFixed(2)} / ${R.params.tEnd.toFixed(2)} s`;
  $('scrub').value = k;
  $('liveT').textContent = `t = ${t.toFixed(3)} s`;
  $('l-w12').textContent = `${sgn(R.Y[0][k])}, ${sgn(R.Y[1][k])}`;
  $('l-w3').textContent = `${sgn(R.Y[2][k], 3)} rad/s`;
  $('l-tilt').textContent = `${R.D.tilt[k].toFixed(3)}°`;
  $('l-n').textContent = (R.D.N[k] / mg).toFixed(3);
  $('l-f').textContent = `${R.D.F[k].toFixed(3)} N`;
  $('l-slip').textContent = `${(1000 * R.D.slip[k]).toFixed(2)} mm/s`;
  $('l-e').textContent = `${R.D.E[k].toFixed(5)} J`;
  setMarks();
}
const seekTime = (t) => { if (!current || !Number.isFinite(t)) return; pause(); setFrame(t / current.params.dtOut); };
function play() {
  if (!current) return;
  if (playhead >= current.T.length - 1) playhead = 0;
  playing = true; simT = current.T[playhead]; lastTs = 0;
  $('playTxt').textContent = 'Pause'; $('playIcon').innerHTML = '<rect x="3.5" y="2.5" width="3" height="11" fill="currentColor"/><rect x="9.5" y="2.5" width="3" height="11" fill="currentColor"/>';
  requestAnimationFrame(tick);
}
function pause() {
  playing = false; $('playTxt').textContent = 'Play'; $('playIcon').innerHTML = '<path d="M4 2.5v11l9-5.5z" fill="currentColor"/>';
}
function tick(ts) {
  if (!playing || !current) return;
  if (lastTs) simT += ((ts - lastTs) / 1000) * parseFloat($('speed').value);
  lastTs = ts;
  let i = Math.round(simT / current.params.dtOut);
  if (i >= current.T.length - 1) {
    if ($('bLoop').getAttribute('aria-pressed') === 'true') { simT = 0; i = 0; }
    else { setFrame(current.T.length - 1); pause(); return; }
  }
  setFrame(i); requestAnimationFrame(tick);
}
const step = (ms) => { if (!current) return; pause(); setFrame(playhead + ms / 1000 / current.params.dtOut); };
$('bPlay').addEventListener('click', () => (playing ? pause() : play()));
$('bBack').addEventListener('click', () => step(-10));
$('bFwd').addEventListener('click', () => step(10));
$('bStart').addEventListener('click', () => { pause(); setFrame(0); });
$('bEnd').addEventListener('click', () => { pause(); setFrame(1e9); });
$('bRev').addEventListener('click', () => current && current.S.tCross != null && seekTime(current.S.tCross));
$('bLoop').addEventListener('click', (e) => e.currentTarget.setAttribute('aria-pressed', e.currentTarget.getAttribute('aria-pressed') !== 'true'));
$('scrub').addEventListener('input', () => { pause(); setFrame(parseInt($('scrub').value, 10)); });
document.querySelectorAll('[data-view]').forEach(b => b.addEventListener('click', () => stage && stage.setView(b.dataset.view)));
function toggle(id, fn) { $(id).addEventListener('click', () => { const on = $(id).getAttribute('aria-pressed') !== 'true'; $(id).setAttribute('aria-pressed', on); fn(on); }); }
toggle('tFollow', on => stage && (stage.follow = on));
toggle('tTrail', on => stage && stage.setTrail(on));
toggle('tContact', on => { if (stage) { stage.setContact(on); setFrame(playhead); } });
$('exag').addEventListener('change', () => { if (stage) { stage.exag = parseFloat($('exag').value); setFrame(playhead); } });

// =====================================================================  keyboard
document.addEventListener('keydown', (e) => {
  const typing = /INPUT|SELECT|TEXTAREA/.test(e.target.tagName) && e.target.type !== 'range' && e.target.type !== 'checkbox';
  if ((e.ctrlKey || e.metaKey) && e.key === 'Enter') { e.preventDefault(); run(); return; }
  if (e.key === 'Escape' && busy) { stop(); return; }
  if (typing || e.ctrlKey || e.metaKey || e.altKey || $('helpDlg').open) return;
  if (e.target.type === 'range' && e.target.id === 'scrub' && /Arrow|Home|End/.test(e.key)) return;
  const k = e.key;
  if (k === ' ' && !/BUTTON|SUMMARY|A/.test(e.target.tagName)) { e.preventDefault(); playing ? pause() : play(); }
  else if (k === 'ArrowRight') { e.preventDefault(); step(e.shiftKey ? 100 : 10); }
  else if (k === 'ArrowLeft') { e.preventDefault(); step(e.shiftKey ? -100 : -10); }
  else if (k === 'Home') { pause(); setFrame(0); }
  else if (k === 'End') { pause(); setFrame(1e9); }
  else if (k === 'r' || k === 'R') $('bRev').click();
  else if (k === '?') $('helpDlg').showModal();
  else if (stage && '1234'.includes(k)) stage.setView(['iso', 'front', 'side', 'top'][+k - 1]);
});
$('helpBtn').addEventListener('click', () => $('helpDlg').showModal());
$('helpClose').addEventListener('click', () => $('helpDlg').close());

// =====================================================================  tabs, splitter, menus, theme
const tabs = [...document.querySelectorAll('.tabs [role=tab]')];
function selectTab(t) {
  tabs.forEach(b => { const on = b === t; b.setAttribute('aria-selected', on); b.tabIndex = on ? 0 : -1; $(b.getAttribute('aria-controls')).hidden = !on; });
  requestAnimationFrame(buildVisible);
}
tabs.forEach(b => b.addEventListener('click', () => selectTab(b)));
document.querySelector('.tabs').addEventListener('keydown', (e) => {
  const i = tabs.indexOf(document.activeElement); if (i < 0) return;
  if (e.key === 'ArrowRight' || e.key === 'ArrowLeft') { e.preventDefault(); const n = tabs[(i + (e.key === 'ArrowRight' ? 1 : tabs.length - 1)) % tabs.length]; n.focus(); selectTab(n); }
});
(function splitter() {
  const sp = $('splitter'), center = document.querySelector('.center');
  const setH = (h) => { const max = center.clientHeight - 220; center.style.setProperty('--dock', Math.max(120, Math.min(max, h)) + 'px'); };
  sp.addEventListener('pointerdown', (e) => {
    sp.setPointerCapture(e.pointerId);
    const move = (ev) => setH(center.getBoundingClientRect().bottom - ev.clientY - 3);
    const up = () => { sp.removeEventListener('pointermove', move); sp.removeEventListener('pointerup', up); };
    sp.addEventListener('pointermove', move); sp.addEventListener('pointerup', up);
  });
  sp.addEventListener('keydown', (e) => {
    const cur = document.querySelector('.dock').getBoundingClientRect().height;
    if (e.key === 'ArrowUp') { e.preventDefault(); setH(cur + 30); } if (e.key === 'ArrowDown') { e.preventDefault(); setH(cur - 30); }
  });
})();
document.addEventListener('click', (e) => { const m = $('exportMenu'); if (m.open && !m.contains(e.target)) m.open = false; });
$('themeBtn').addEventListener('click', () => {
  const cur = document.documentElement.dataset.theme || (matchMedia('(prefers-color-scheme: dark)').matches ? 'dark' : 'light');
  const next = cur === 'dark' ? 'light' : 'dark';
  document.documentElement.dataset.theme = next;
  try { localStorage.setItem('rb-theme', next); } catch {}
  if (stage) stage.setTheme(next === 'dark');
  chartsStale.hist = chartsStale.cmp = true; plots = Object.fromEntries(Object.entries(plots).filter(([k, p]) => (p.u.destroy(), false)));
  buildVisible();
});

// =====================================================================  export
function download(name, data, type) {
  const a = document.createElement('a');
  a.href = typeof data === 'string' && data.startsWith('data:') ? data : URL.createObjectURL(new Blob([data], { type }));
  a.download = name; document.body.append(a); a.click(); a.remove();
  if (!a.href.startsWith('data:')) setTimeout(() => URL.revokeObjectURL(a.href), 3000);
}
const slug = (s) => s.replace(/[^a-z0-9]+/gi, '_').replace(/^_|_$/g, '').toLowerCase().slice(0, 60);
function hashOf(f) { return '#' + new URLSearchParams(Object.entries(f).map(([k, v]) => [k, String(v)])).toString(); }
function fromHash() {
  if (location.hash.length < 2) return null;
  const q = new URLSearchParams(location.hash.slice(1)), f = { ...BASE };
  for (const [k, v] of q) if (k in f) f[k] = NUMERIC.includes(k) ? parseFloat(v) : k === 'drill' ? v === 'true' : v;
  return f;
}
$('exportMenu').addEventListener('click', async (e) => {
  const act = e.target.closest('[data-act]')?.dataset.act; if (!act) return;
  $('exportMenu').open = false;
  const R = current;
  if (act === 'csv' && R) {
    const mg = R.params.m * R.params.g, lines = [`# Rattleback simulator, run ${R.id}: ${R.label}`, 't_s,w1_rad_s,w2_rad_s,w3_rad_s,phi_rad,theta_rad,psi_rad,v1_m_s,v2_m_s,v3_m_s,xG_m,yG_m,N_over_mg,F_N,slip_m_s,E_J,tilt_deg'];
    for (let i = 0; i < R.T.length; i++) lines.push([R.T[i], ...R.Y.map(a => a[i]), R.D.N[i] / mg, R.D.F[i], R.D.slip[i], R.D.E[i], R.D.tilt[i]].map(x => +x.toPrecision(8)).join(','));
    download(`rattleback_run${R.id}_${slug(R.label)}.csv`, lines.join('\n'), 'text/csv'); toast('Time histories downloaded.');
  } else if (act === 'summary') {
    const head = 'run,setup,reversal_s,spin_down,w3_end_rad_s,theta_max_deg,energy_lost_pct,N_min_mg,N_max_mg,slip_max_mm_s,F_over_N_max,steps,seconds';
    const rows = runHist.map(h => [h.id, `"${h.label}"`, h.S.reversal ? h.S.tCross.toFixed(4) : '', h.S.spinDown, h.S.w3End.toFixed(5), h.S.thetaMax.toFixed(4), (-100 * h.S.dEfrac).toFixed(3), h.S.NminMg.toFixed(4), h.S.NmaxMg.toFixed(4), (1000 * h.S.slipMax).toFixed(3), h.S.ratioMax.toFixed(4), h.stats.nSteps, h.secs.toFixed(3)].join(','));
    download('rattleback_runs.csv', [head, ...rows].join('\n'), 'text/csv'); toast('Run history downloaded.');
  } else if (act === 'png' && stage) {
    download(`rattleback_run${R ? R.id : ''}_t${R ? R.T[playhead].toFixed(2) : 0}s.png`, stage.snapshot(), 'image/png'); toast('Snapshot downloaded.');
  } else if (act === 'json') {
    download('rattleback_setup.json', JSON.stringify({ setup: readForm(), parameters_SI: toParams(readForm()) }, null, 2), 'application/json'); toast('Settings downloaded.');
  } else if (act === 'link') {
    const url = location.href.split('#')[0] + hashOf(readForm()); window.history.replaceState(null, '', url);
    try { await navigator.clipboard.writeText(url); toast('Link copied. It opens the simulator with these settings.'); }
    catch { toast('The link is in the address bar.'); }
  }
});

// =====================================================================  sweep
const SWEEPS = [
  { key: 'mu', label: 'Friction coefficient μ', unit: '', from: 0.05, to: 1, n: 8, scale: 'lin', patch: (v) => ({ law: 'coulomb', mu: v }), note: 'Uses the Coulomb law.' },
  { key: 'rho', label: 'Drilling lever arm ρ', unit: 'mm', from: 0, to: 2, n: 9, scale: 'lin', patch: (v) => ({ law: 'coulomb', rho: v * 1e-3 }), note: 'Uses the Coulomb law with drilling friction.' },
  { key: 'cv', label: 'Viscous coefficient c', unit: 'kg/s', from: 10, to: 3000, n: 8, scale: 'log', patch: (v) => ({ law: 'viscous', cv: v }), note: 'Uses the viscous law.' },
  { key: 'w0', label: 'Initial spin |ω₀|', unit: 'rad/s', from: 2, to: 10, n: 9, scale: 'lin', patch: (v, base) => ({ w0: Math.sign(base.w0) * v }), note: '' },
  { key: 'eps', label: 'Slip regularization ε', unit: 'm/s', from: 1e-5, to: 1e-3, n: 7, scale: 'log', patch: (v) => ({ law: 'coulomb', eps: v }), note: 'Uses the Coulomb law. A convergence check.' },
  { key: 'I12', label: 'Product of inertia I₁₂', unit: '×10⁻⁴ kg m²', from: -1, to: -0.1, n: 10, scale: 'lin', patch: (v) => ({ I12: v * 1e-4 }), note: 'The spin sense is kept as set.' },
];
const METRICS = {
  tCross: { label: 'Reversal time', unit: 's', get: (S) => (S.reversal ? S.tCross : null) },
  w3End: { label: 'Spin at the end', unit: 'rad/s', get: (S) => S.w3End },
  thetaMax: { label: 'Largest rocking angle', unit: 'deg', get: (S) => S.thetaMax },
  dE: { label: 'Energy lost', unit: '%', get: (S) => -100 * S.dEfrac },
  NmaxMg: { label: 'Largest N/mg', unit: '', get: (S) => S.NmaxMg },
  slipMax: { label: 'Largest slip speed', unit: 'mm/s', get: (S) => 1000 * S.slipMax },
};
const sweep = { def: null, rows: [], values: [] };
$('swParam').innerHTML = SWEEPS.map(s => `<option value="${s.key}">${s.label}${s.unit ? ` (${s.unit})` : ''}</option>`).join('');
function loadSweepDef() {
  const d = SWEEPS.find(s => s.key === $('swParam').value);
  $('swFrom').value = d.from; $('swTo').value = d.to; $('swN').value = d.n; $('swScale').value = d.scale;
  $('swNote').textContent = d.note;
}
$('swParam').addEventListener('change', loadSweepDef); loadSweepDef();
const fmtX = (x) => (Math.abs(x) < 1e-3 && x !== 0 ? x.toExponential(1) : +x.toPrecision(4) + '');
$('sweepForm').addEventListener('submit', (e) => {
  e.preventDefault();
  const f = readForm(), errs = validate(f);
  if (Object.keys(errs).some(k => shown(ALL.find(d => d.id === k), f))) { status('Check the setup', 'Fix the highlighted setup fields first.', 'err'); return; }
  const d = SWEEPS.find(s => s.key === $('swParam').value);
  const a = parseFloat($('swFrom').value), b = parseFloat($('swTo').value), n = parseInt($('swN').value, 10), isLog = $('swScale').value === 'log';
  if (!Number.isFinite(a) || !Number.isFinite(b) || a === b) { toast('Enter two different end values.'); return; }
  if (!(n >= 2 && n <= 40)) { toast('Use 2 to 40 points.'); return; }
  if (isLog && (a <= 0 || b <= 0)) { toast('Logarithmic spacing needs positive end values.'); return; }
  const values = Array.from({ length: n }, (_, i) => (isLog ? Math.exp(Math.log(a) + (Math.log(b) - Math.log(a)) * i / (n - 1)) : a + (b - a) * i / (n - 1)));
  const base = { ...toParams(f), dtOut: Math.max(0.002, parseFloat(f.dtOut)) };
  if (d.key === 'rho' && !base.rho) base.rho = 0;
  const patches = values.map(v => d.patch(v, base));
  if (busy) stop();
  sweep.def = d; sweep.values = values; sweep.rows = []; sweep.log = isLog;
  $('swCsv').disabled = true;
  const id = ++job; setBusy('sweep'); status('Sweeping', `${d.label}, ${n} runs`, 'busy');
  log(`Sweep of ${d.label} from ${fmtX(a)} to ${fmtX(b)}, ${n} points${isLog ? ', logarithmic' : ''}.`);
  getWorker().postMessage({ type: 'sweep', id, base, patches, values });
  renderSweep();
});
$('swMetric').addEventListener('change', renderSweep);
function renderSweep() {
  const d = sweep.def; if (!d) return;
  const M = METRICS[$('swMetric').value];
  const unit = d.unit ? ` (${d.unit})` : '';
  $('swTh').textContent = d.label.split(' ').slice(-1)[0] + unit;
  $('swTitle').textContent = `${M.label} against ${d.label.charAt(0).toLowerCase() + d.label.slice(1)}`;
  const done = sweep.rows.filter(Boolean).length;
  $('swSub').textContent = `${done} of ${sweep.values.length} runs${M.unit ? `; ${M.label.toLowerCase()} in ${M.unit}` : ''}${$('swMetric').value === 'tCross' ? '; gaps mean no reversal' : ''}`;
  $('swTbl').tBodies[0].innerHTML = sweep.values.map((x, i) => {
    const r = sweep.rows[i];
    if (!r) return `<tr><td class="n">${fmtX(x)}</td><td class="n" colspan="5">${i === done && busy === 'sweep' ? 'computing…' : 'waiting'}</td></tr>`;
    if (!r.ok) return `<tr><td class="n">${fmtX(x)}</td><td colspan="5">${esc(r.error)}</td></tr>`;
    const S = r.S;
    return `<tr><td class="n">${fmtX(x)}</td><td class="n">${S.reversal ? S.tCross.toFixed(3) + ' s' : S.spinDown ? 'spin-down' : 'none'}</td><td class="n">${sgn(S.w3End, 3)}</td><td class="n">${S.thetaMax.toFixed(3)}°</td><td class="n">${(-100 * S.dEfrac).toFixed(2)}%</td><td class="n">${S.NmaxMg.toFixed(3)}</td></tr>`;
  }).join('');
  const el = $('p-sweep'); if (!el.clientWidth) return;
  const xs = [], ys = [];
  sweep.values.forEach((x, i) => { const r = sweep.rows[i]; xs.push(x); ys.push(r && r.ok ? M.get(r.S) : null); });
  makePlot('sweep', el, [{ label: M.label, stroke: css('--s1'), width: 1.5, points: { show: true, size: 7, fill: css('--field') }, value: (u, v) => (v == null ? '–' : v.toFixed(3) + (M.unit ? ' ' + M.unit : '')) }],
    [xs, ys], `${M.label}${M.unit ? ` (${M.unit})` : ''}`, { xLabel: d.label + unit, xFmt: fmtX, xLog: sweep.log, noSync: true, noMark: true, noSeek: true, height: 230 });
}
$('swCsv').addEventListener('click', () => {
  const d = sweep.def; if (!d) return;
  const head = `${d.key}${d.unit ? '_' + slug(d.unit) : ''},reversal_s,spin_down,w3_end_rad_s,theta_max_deg,energy_lost_pct,N_max_mg,slip_max_mm_s,error`;
  const rows = sweep.values.map((x, i) => { const r = sweep.rows[i]; if (!r) return `${x},,,,,,,,not computed`; if (!r.ok) return `${x},,,,,,,,"${r.error}"`; const S = r.S;
    return [x, S.reversal ? S.tCross.toFixed(4) : '', S.spinDown, S.w3End.toFixed(5), S.thetaMax.toFixed(4), (-100 * S.dEfrac).toFixed(3), S.NmaxMg.toFixed(4), (1000 * S.slipMax).toFixed(3), ''].join(','); });
  download(`rattleback_sweep_${d.key}.csv`, [`# base setup: ${describe(readForm())}`, head, ...rows].join('\n'), 'text/csv'); toast('Sweep results downloaded.');
});

// =====================================================================  start
(async function start() {
  const linked = fromHash();
  if (linked) { writeForm(linked); presetSel.value = 'custom'; log('Setup loaded from a shared link.'); }
  else applyPreset('coulomb', true);
  log('Ready. Equations as in the paper; normal reaction solved from the contact condition.');
  try {
    const { Stage } = await import('./scene.js');
    stage = new Stage($('viewport'));
    $('vpEmpty').textContent = 'Press Run to compute the motion.';
  } catch (err) {
    $('vpEmpty').textContent = 'The 3D view needs WebGL, which is turned off in this browser. Results and charts still work.';
    log('3D view unavailable: ' + err.message, true);
  }
  run();
})();
