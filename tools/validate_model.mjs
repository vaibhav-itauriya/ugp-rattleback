// Checks the browser model against the paper's MATLAB/Python results.
// Run from the repo root:  node tools/validate_model.mjs
import { simulate } from '../docs/js/model.js';
const cases = [
  ['Rolling',              { law: 'rolling' },                 { t: 2.831, w: 5.000, th: 6.48, Nmax: 5.86 }],
  ['Coulomb mu=0.15',      { law: 'coulomb' },                 { t: 4.223, w: 1.012, th: 0.836, Nmax: 1.301 }],
  ['Coulomb + rho=0.5mm',  { law: 'coulomb', rho: 0.5e-3 },    { t: 3.977, w: 0.018, th: 0.840 }],
  ['Viscous c=500',        { law: 'viscous' },                 { t: 2.898, w: 1.396, th: 3.09, Nmax: 2.63 }],
  ['Coulomb mu=1, N=mg',   { law: 'coulomb', mu: 1, fixN: true }, { t: 2.428 }],
  ['Coulomb pref',         { law: 'coulomb', w0: 5 },          { t: null, w: 5.0 }],
];
let bad = 0;
for (const [name, p, ref] of cases) {
  const R = simulate(p);
  const S = R.S;
  const row = { t: S.tCross, w: S.w3End, th: S.thetaMax, Nmax: S.NmaxMg };
  const fmt = x => x == null ? 'none' : x.toFixed(3);
  const ok = Object.entries(ref).every(([k, v]) => v == null ? row[k] == null : Math.abs(row[k] - v) <= Math.max(0.01, 0.01 * Math.abs(v)));
  if (!ok) bad++;
  console.log(`${ok ? 'ok  ' : 'FAIL'} ${name.padEnd(22)} t_rev ${fmt(row.t)} (${fmt(ref.t)})  w3 ${fmt(row.w)} (${fmt(ref.w)})  thmax ${fmt(row.th)} (${fmt(ref.th)})  Nmax/mg ${fmt(row.Nmax)}  dE ${(100*S.dEfrac).toFixed(2)}%  ${R.stats.ms.toFixed(0)} ms, ${R.stats.nSteps} steps`);
}
process.exit(bad ? 1 : 0);
