// Rattleback on a plane: rolling, or slipping with Coulomb / drilling / viscous friction.
// JavaScript port of code/matlab/rattleback_simulator_123H.m (same equations, same
// 1-2-3 Euler angles). State y = [w1 w2 w3 phi theta psi v1 v2 v3 xG yG],
// angular velocity and CM velocity in body components, xG yG in the fixed frame.
// The normal reaction is solved from the contact condition (paper, Eq. 18):
//   N = m (g + c + d.A0) / (1 - m d.B - m d.A1)

export const DEFAULTS = {
  m: 1.0, g: 9.81,
  a1: 0.20, a2: 0.03, a3: 0.02, b: 0.005,
  I11: 2.0e-4, I22: 1.6e-3, I33: 1.75e-3, I12: -0.5e-4,
  law: 'coulomb',          // 'rolling' | 'coulomb' | 'viscous'
  mu: 0.15, eps: 3e-5, rho: 0.0, cv: 500, fixN: false,
  w0: -5.0, phi0: 1e-6, theta0: 1e-6, psi0: 1e-6,
  tEnd: 8.0, dtOut: 0.001, rtol: 1e-8, atol: 1e-10,
};

// ---------------------------------------------------------------- body ----
export function makeBody(p) {
  const I = [[p.I11, p.I12, 0], [p.I12, p.I22, 0], [0, 0, p.I33]];
  const det = I[0][0] * I[1][1] - I[0][1] * I[1][0];
  if (!(det > 0) || !(p.I33 > 0)) throw new Error('The inertia tensor must be positive definite (I11*I22 > I12^2).');
  const Ii = [[I[1][1] / det, -I[0][1] / det, 0], [-I[1][0] / det, I[0][0] / det, 0], [0, 0, 1 / I[2][2]]];
  return { ...p, I, Ii, A2: [p.a1 * p.a1, p.a2 * p.a2, p.a3 * p.a3] };
}

const cross = (a, b) => [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]];
const dot = (a, b) => a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
const mv = (M, v) => [M[0][0] * v[0] + M[0][1] * v[1] + M[0][2] * v[2],
                      M[1][0] * v[0] + M[1][1] * v[1] + M[1][2] * v[2],
                      M[2][0] * v[0] + M[2][1] * v[1] + M[2][2] * v[2]];

export function rotmat(phi, th, ps) {
  const c1 = Math.cos(phi), s1 = Math.sin(phi), c2 = Math.cos(th), s2 = Math.sin(th), c3 = Math.cos(ps), s3 = Math.sin(ps);
  // R1(phi) R2(theta) R3(psi): body -> space
  return [
    [c2 * c3, -c2 * s3, s2],
    [c1 * s3 + s1 * s2 * c3, c1 * c3 - s1 * s2 * s3, -s1 * c2],
    [s1 * s3 - c1 * s2 * c3, s1 * c3 + c1 * s2 * s3, c1 * c2],
  ];
}

// geometry: vertical n in body frame, contact point r = r_{P/G}, their rates, height h
export function geom(B, phi, th, ps, w) {
  const n = [Math.sin(phi) * Math.sin(ps) - Math.cos(ps) * Math.sin(th) * Math.cos(phi),
             Math.cos(ps) * Math.sin(phi) + Math.sin(ps) * Math.sin(th) * Math.cos(phi),
             Math.cos(th) * Math.cos(phi)];
  const nd = cross(w, n).map(x => -x);
  const A2 = B.A2;
  const S = A2[0] * n[0] * n[0] + A2[1] * n[1] * n[1] + A2[2] * n[2] * n[2];
  const sq = Math.sqrt(S);
  const r = [-A2[0] * n[0] / sq, -A2[1] * n[1] / sq, -A2[2] * n[2] / sq + B.b];
  const Sd = A2[0] * n[0] * nd[0] + A2[1] * n[1] * nd[1] + A2[2] * n[2] * nd[2];
  const S15 = S * sq;
  const rd = [0, 1, 2].map(i => -A2[i] * nd[i] / sq + A2[i] * n[i] * Sd / S15);
  return { n, nd, r, rd, h: sq - B.b * n[2] };
}

function eulerRates(w, th, ps) {
  const ct = Math.cos(th), st = Math.sin(th), cp = Math.cos(ps), sp = Math.sin(ps);
  const phid = (w[0] * cp - w[1] * sp) / ct;
  return [phid, w[0] * sp + w[1] * cp, w[2] - st * phid];
}

// --------------------------------------------------- slipping contact ----
export function contact(B, w, ang, v) {
  const { n, nd, r, rd } = geom(B, ang[0], ang[1], ang[2], w);
  let vp = [v[0] + w[1] * r[2] - w[2] * r[1], v[1] + w[2] * r[0] - w[0] * r[2], v[2] + w[0] * r[1] - w[1] * r[0]];
  const vn = dot(vp, n); vp = [vp[0] - vn * n[0], vp[1] - vn * n[1], vp[2] - vn * n[2]];
  const d = cross(n, r);
  const c = -dot(rd, nd) + dot(r, cross(w, nd));
  const Iw = mv(B.I, w), gyro = cross(w, Iw);
  let A0 = mv(B.Ii, [-gyro[0], -gyro[1], -gyro[2]]);
  const Bv = mv(B.Ii, cross(r, n));
  let A1 = [0, 0, 0], Fv = null, uhat = [0, 0, 0], sgn = 0;
  const vp2 = dot(vp, vp);
  if (B.law === 'viscous') {
    Fv = [-B.cv * vp[0], -B.cv * vp[1], -B.cv * vp[2]];
    const t = cross(r, Fv);
    A0 = mv(B.Ii, [t[0] - gyro[0], t[1] - gyro[1], t[2] - gyro[2]]);
  } else {
    const s = Math.sqrt(vp2 + B.eps * B.eps);
    uhat = [vp[0] / s, vp[1] / s, vp[2] / s];
    let tau = cross(r, [-B.mu * uhat[0], -B.mu * uhat[1], -B.mu * uhat[2]]);
    if (B.rho > 0) {
      const sc = B.rho * dot(w, n);
      sgn = sc / Math.sqrt(sc * sc + vp2 + B.eps * B.eps);
      tau = [tau[0] - B.mu * B.rho * sgn * n[0], tau[1] - B.mu * B.rho * sgn * n[1], tau[2] - B.mu * B.rho * sgn * n[2]];
    }
    A1 = mv(B.Ii, tau);
  }
  const den = 1 - B.m * dot(d, Bv) - B.m * dot(d, A1);
  let N = B.m * (B.g + c + dot(d, A0)) / den;
  if (B.fixN) N = B.m * B.g;
  const F = Fv || [-B.mu * N * uhat[0], -B.mu * N * uhat[1], -B.mu * N * uhat[2]];
  const Tsp = B.rho > 0 && !Fv ? B.mu * N * B.rho * Math.abs(sgn) : 0;
  const wd = [A0[0] + N * (Bv[0] + A1[0]), A0[1] + N * (Bv[1] + A1[1]), A0[2] + N * (Bv[2] + A1[2])];
  const vd = [(N * n[0] + F[0]) / B.m - B.g * n[0] - (w[1] * v[2] - w[2] * v[1]),
              (N * n[1] + F[1]) / B.m - B.g * n[1] - (w[2] * v[0] - w[0] * v[2]),
              (N * n[2] + F[2]) / B.m - B.g * n[2] - (w[0] * v[1] - w[1] * v[0])];
  const gv = dot(v, n) + dot(r, nd);
  return { wd, vd, N, F, vp, den, gv, Tsp };
}

// ------------------------------------------------------ rolling model ----
function solve3(M, b) {
  const [a, bb, c] = M[0], [d, e, f] = M[1], [g, h, i] = M[2];
  const A = e * i - f * h, Bc = -(d * i - f * g), C = d * h - e * g;
  const det = a * A + bb * Bc + c * C;
  const inv = [[A, -(bb * i - c * h), bb * f - c * e], [Bc, a * i - c * g, -(a * f - c * d)], [C, -(a * h - bb * g), a * e - bb * d]];
  return mv(inv, b).map(x => x / det);
}
export function rollingWd(B, w, ang) {
  const { n, r, rd } = geom(B, ang[0], ang[1], ang[2], w);
  const X = [[0, -r[2], r[1]], [r[2], 0, -r[0]], [-r[1], r[0], 0]];
  const XX = [0, 1, 2].map(i => [0, 1, 2].map(j => X[i][0] * X[0][j] + X[i][1] * X[1][j] + X[i][2] * X[2][j]));
  const Ir = [0, 1, 2].map(i => [0, 1, 2].map(j => B.I[i][j] - B.m * XX[i][j]));
  const Iw = mv(B.I, w);
  const g1 = cross(w, Iw), wr = cross(w, r), wwr = cross(w, wr), wrd = cross(w, rd);
  const t1 = cross(r, wwr), t2 = cross(r, wrd), t3 = cross(r, n);
  const rhs = [0, 1, 2].map(i => -g1[i] - B.m * t1[i] - B.m * t2[i] + B.m * B.g * t3[i]);
  return solve3(Ir, rhs);
}
function rollingContact(B, w, ang) {
  const wd = rollingWd(B, w, ang);
  const { n, r, rd } = geom(B, ang[0], ang[1], ang[2], w);
  const aG = [0, 1, 2].map(i => -cross(wd, r)[i] - cross(w, cross(w, r))[i] - cross(w, rd)[i]);
  const fP = [0, 1, 2].map(i => B.m * aG[i] + B.m * B.g * n[i]);
  const N = dot(fP, n);
  const F = [fP[0] - N * n[0], fP[1] - N * n[1], fP[2] - N * n[2]];
  return { wd, N, F, vp: [0, 0, 0], gv: 0, den: 1, Tsp: 0 };
}

// --------------------------------------------------------------- rhs ----
export function rhs(B, y) {
  const w = [y[0], y[1], y[2]], ang = [y[3], y[4], y[5]];
  let wd, vd, v;
  if (B.law === 'rolling') {
    wd = rollingWd(B, w, ang);
    const { r } = geom(B, ang[0], ang[1], ang[2], w);
    v = cross(w, r).map(x => -x); vd = [0, 0, 0];
  } else {
    v = [y[6], y[7], y[8]];
    ({ wd, vd } = contact(B, w, ang, v));
  }
  const er = eulerRates(w, ang[1], ang[2]);
  const Q = rotmat(ang[0], ang[1], ang[2]);
  const vs = mv(Q, v);
  return [wd[0], wd[1], wd[2], er[0], er[1], er[2], vd[0], vd[1], vd[2], vs[0], vs[1]];
}

export function initialState(B) {
  const w = [0, 0, B.w0];
  const { r } = geom(B, B.phi0, B.theta0, B.psi0, w);
  const v = cross(w, r).map(x => -x);          // start rolling: v_P(0) = 0
  return [w[0], w[1], w[2], B.phi0, B.theta0, B.psi0, v[0], v[1], v[2], 0, 0];
}

// ------------------------------------------- Dormand-Prince 5(4) ----
const C = [0, 1 / 5, 3 / 10, 4 / 5, 8 / 9, 1, 1];
const A = [[], [1 / 5], [3 / 40, 9 / 40], [44 / 45, -56 / 15, 32 / 9],
  [19372 / 6561, -25360 / 2187, 64448 / 6561, -212 / 729],
  [9017 / 3168, -355 / 33, 46732 / 5247, 49 / 176, -5103 / 18656],
  [35 / 384, 0, 500 / 1113, 125 / 192, -2187 / 6784, 11 / 84]];
const B5 = [35 / 384, 0, 500 / 1113, 125 / 192, -2187 / 6784, 11 / 84, 0];
const E = [71 / 57600, 0, -71 / 16695, 71 / 1920, -17253 / 339200, 22 / 525, -1 / 40];

// Integrate to tEnd and return samples every dtOut. The regularized Coulomb law is
// stiff (Jacobian ~ mu N / eps); the adaptive step simply shrinks to stay stable.
export function integrate(B, onProgress) {
  const nY = 11, y0 = initialState(B);
  const nOut = Math.round(B.tEnd / B.dtOut) + 1;
  const T = new Float64Array(nOut), Y = Array.from({ length: nY }, () => new Float64Array(nOut));
  let y = y0.slice(), t = 0, h = 1e-5, k1 = rhs(B, y), nSteps = 0, nRej = 0, iOut = 0;
  const store = () => { T[iOut] = t; for (let j = 0; j < nY; j++) Y[j][iOut] = y[j]; iOut++; };
  store();
  const k = new Array(7);
  let lastReport = 0;
  while (iOut < nOut) {
    const tNext = iOut * B.dtOut;
    if (t + h > tNext) h = tNext - t;
    k[0] = k1;
    for (let s = 1; s < 7; s++) {
      const ys = new Array(nY);
      for (let j = 0; j < nY; j++) { let acc = y[j]; for (let q = 0; q < s; q++) acc += h * A[s][q] * k[q][j]; ys[j] = acc; }
      k[s] = rhs(B, ys);
    }
    let err = 0; const yNew = new Array(nY);
    for (let j = 0; j < nY; j++) {
      let acc = y[j]; for (let q = 0; q < 6; q++) acc += h * B5[q] * k[q][j];
      yNew[j] = acc;
      let e = 0; for (let q = 0; q < 7; q++) e += h * E[q] * k[q][j];
      const sc = B.atol + B.rtol * Math.max(Math.abs(y[j]), Math.abs(acc));
      err = Math.max(err, Math.abs(e) / sc);
    }
    if (!isFinite(err)) throw new Error('The integration diverged. Try a larger regularization velocity or smaller friction.');
    if (err <= 1) {
      t += h; y = yNew; k1 = k[6]; nSteps++;
      if (Math.abs(t - tNext) < 1e-12) { t = tNext; store(); }
      if (onProgress && t - lastReport > 0.05) { lastReport = t; onProgress(t / B.tEnd); }
    } else nRej++;
    const fac = err === 0 ? 5 : Math.min(5, Math.max(0.2, 0.9 * Math.pow(err, -0.2)));
    h = Math.max(h * fac, 1e-12);
    if (nSteps + nRej > 5e7) throw new Error('Too many steps. Reduce the simulated time or increase the regularization velocity.');
  }
  return { T, Y, nSteps, nRej };
}

// ------------------------------------------------------- diagnostics ----
export function diagnostics(B, T, Y) {
  const n = T.length, out = {
    N: new Float64Array(n), F: new Float64Array(n), slip: new Float64Array(n), E: new Float64Array(n),
    gv: new Float64Array(n), den: new Float64Array(n), tilt: new Float64Array(n), h: new Float64Array(n),
  };
  for (let i = 0; i < n; i++) {
    const w = [Y[0][i], Y[1][i], Y[2][i]], ang = [Y[3][i], Y[4][i], Y[5][i]];
    const g = geom(B, ang[0], ang[1], ang[2], w);
    let v, cr;
    if (B.law === 'rolling') { v = cross(w, g.r).map(x => -x); cr = rollingContact(B, w, ang); }
    else { v = [Y[6][i], Y[7][i], Y[8][i]]; cr = contact(B, w, ang, v); }
    out.N[i] = cr.N; out.F[i] = Math.hypot(cr.F[0], cr.F[1], cr.F[2]);
    out.slip[i] = Math.hypot(cr.vp[0], cr.vp[1], cr.vp[2]);
    out.gv[i] = cr.gv; out.den[i] = cr.den; out.h[i] = g.h;
    out.tilt[i] = Math.acos(Math.min(1, g.n[2])) * 180 / Math.PI;   // angle of e3 from vertical
    const Iw = mv(B.I, w);
    out.E[i] = 0.5 * dot(w, Iw) + 0.5 * B.m * dot(v, v) + B.m * B.g * g.h;
  }
  return out;
}

export function summarize(B, T, Y, D) {
  const w3 = Y[2], n = T.length, s0 = Math.sign(B.w0) || 1;
  let iCross = -1;
  for (let i = 0; i < n; i++) if (T[i] > 0.05 && Math.sign(w3[i]) === -s0) { iCross = i; break; }
  let settled = 0; if (iCross >= 0) for (let i = iCross; i < n; i++) settled = Math.max(settled, Math.abs(w3[i]));
  let thMax = 0, Nmin = Infinity, Nmax = -Infinity, slipMax = 0, ratioMax = 0, gvMax = 0, denMin = Infinity, eMono = true;
  for (let i = 0; i < n; i++) {
    thMax = Math.max(thMax, Math.abs(Y[4][i]) * 180 / Math.PI);
    Nmin = Math.min(Nmin, D.N[i]); Nmax = Math.max(Nmax, D.N[i]);
    slipMax = Math.max(slipMax, D.slip[i]); gvMax = Math.max(gvMax, Math.abs(D.gv[i]));
    denMin = Math.min(denMin, D.den[i]);
    if (D.N[i] > 0) ratioMax = Math.max(ratioMax, D.F[i] / D.N[i]);
    if (i > 0 && D.E[i] - D.E[i - 1] > 1e-9) eMono = false;
  }
  const mg = B.m * B.g;
  return {
    tCross: iCross >= 0 ? T[iCross] : null,
    reversal: iCross >= 0 && settled > 0.05,
    spinDown: iCross >= 0 && settled <= 0.05,
    w3End: w3[n - 1], thetaMax: thMax,
    NminMg: Nmin / mg, NmaxMg: Nmax / mg, slipMax, ratioMax,
    E0: D.E[0], Eend: D.E[n - 1], dEfrac: (D.E[n - 1] - D.E[0]) / D.E[0],
    checks: { residual: gvMax, energyOk: eMono, Npositive: Nmin > 0, denMin },
  };
}

export function simulate(params, onProgress) {
  const B = makeBody({ ...DEFAULTS, ...params });
  const t0 = (typeof performance !== 'undefined' ? performance : Date).now();
  const { T, Y, nSteps, nRej } = integrate(B, onProgress);
  const D = diagnostics(B, T, Y);
  const S = summarize(B, T, Y, D);
  const ms = (typeof performance !== 'undefined' ? performance : Date).now() - t0;
  return { params: { ...DEFAULTS, ...params }, T, Y, D, S, stats: { nSteps, nRej, ms } };
}
