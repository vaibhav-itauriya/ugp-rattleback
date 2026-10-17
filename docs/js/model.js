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

// Allocation-free right-hand side for the slipping models (the hot loop of the
// stiff Coulomb case). Same equations as contact() + rhs() above, written out in
// scalars; tools/validate_model.mjs checks the two agree.
function rhsSlip(B, y, out) {
  const w0 = y[0], w1 = y[1], w2 = y[2], v0 = y[6], v1 = y[7], v2 = y[8];
  const sf = Math.sin(y[3]), cf = Math.cos(y[3]), st = Math.sin(y[4]), ct = Math.cos(y[4]), sp = Math.sin(y[5]), cp = Math.cos(y[5]);
  const n0 = sf * sp - cp * st * cf, n1 = cp * sf + sp * st * cf, n2 = ct * cf;
  const nd0 = -(w1 * n2 - w2 * n1), nd1 = -(w2 * n0 - w0 * n2), nd2 = -(w0 * n1 - w1 * n0);
  const q0 = B.A2[0], q1 = B.A2[1], q2 = B.A2[2];
  const S = q0 * n0 * n0 + q1 * n1 * n1 + q2 * n2 * n2, sq = Math.sqrt(S), S15 = S * sq;
  const r0 = -q0 * n0 / sq, r1 = -q1 * n1 / sq, r2 = -q2 * n2 / sq + B.b;
  const Sd = q0 * n0 * nd0 + q1 * n1 * nd1 + q2 * n2 * nd2;
  const rd0 = -q0 * nd0 / sq + q0 * n0 * Sd / S15, rd1 = -q1 * nd1 / sq + q1 * n1 * Sd / S15, rd2 = -q2 * nd2 / sq + q2 * n2 * Sd / S15;
  let p0 = v0 + w1 * r2 - w2 * r1, p1 = v1 + w2 * r0 - w0 * r2, p2 = v2 + w0 * r1 - w1 * r0;
  const vn = p0 * n0 + p1 * n1 + p2 * n2; p0 -= vn * n0; p1 -= vn * n1; p2 -= vn * n2;
  const d0 = n1 * r2 - n2 * r1, d1 = n2 * r0 - n0 * r2, d2 = n0 * r1 - n1 * r0;
  const c = -(rd0 * nd0 + rd1 * nd1 + rd2 * nd2) + r0 * (w1 * nd2 - w2 * nd1) + r1 * (w2 * nd0 - w0 * nd2) + r2 * (w0 * nd1 - w1 * nd0);
  const I = B.I, Ii = B.Ii;
  const Iw0 = I[0][0] * w0 + I[0][1] * w1, Iw1 = I[1][0] * w0 + I[1][1] * w1, Iw2 = I[2][2] * w2;
  const g0 = w1 * Iw2 - w2 * Iw1, g1 = w2 * Iw0 - w0 * Iw2, g2 = w0 * Iw1 - w1 * Iw0;
  const a00 = Ii[0][0], a01 = Ii[0][1], a11 = Ii[1][1], a22 = Ii[2][2];
  const rn0 = r1 * n2 - r2 * n1, rn1 = r2 * n0 - r0 * n2, rn2 = r0 * n1 - r1 * n0;
  const Bv0 = a00 * rn0 + a01 * rn1, Bv1 = a01 * rn0 + a11 * rn1, Bv2 = a22 * rn2;
  const pp = p0 * p0 + p1 * p1 + p2 * p2;
  let A00, A01, A02, A10 = 0, A11 = 0, A12 = 0, F0, F1, F2, u0 = 0, u1 = 0, u2 = 0;
  const visc = B.law === 'viscous';
  if (visc) {
    F0 = -B.cv * p0; F1 = -B.cv * p1; F2 = -B.cv * p2;
    const t0 = r1 * F2 - r2 * F1 - g0, t1 = r2 * F0 - r0 * F2 - g1, t2 = r0 * F1 - r1 * F0 - g2;
    A00 = a00 * t0 + a01 * t1; A01 = a01 * t0 + a11 * t1; A02 = a22 * t2;
  } else {
    A00 = -(a00 * g0 + a01 * g1); A01 = -(a01 * g0 + a11 * g1); A02 = -a22 * g2;
    const s = Math.sqrt(pp + B.eps * B.eps); u0 = p0 / s; u1 = p1 / s; u2 = p2 / s;
    const mu = B.mu;
    let t0 = -mu * (r1 * u2 - r2 * u1), t1 = -mu * (r2 * u0 - r0 * u2), t2 = -mu * (r0 * u1 - r1 * u0);
    if (B.rho > 0) {
      const sc = B.rho * (w0 * n0 + w1 * n1 + w2 * n2);
      const sg = sc / Math.sqrt(sc * sc + pp + B.eps * B.eps), k = mu * B.rho * sg;
      t0 -= k * n0; t1 -= k * n1; t2 -= k * n2;
    }
    A10 = a00 * t0 + a01 * t1; A11 = a01 * t0 + a11 * t1; A12 = a22 * t2;
  }
  const m = B.m, g = B.g;
  const den = 1 - m * (d0 * Bv0 + d1 * Bv1 + d2 * Bv2) - m * (d0 * A10 + d1 * A11 + d2 * A12);
  let N = m * (g + c + d0 * A00 + d1 * A01 + d2 * A02) / den;
  if (B.fixN) N = m * g;
  if (!visc) { F0 = -B.mu * N * u0; F1 = -B.mu * N * u1; F2 = -B.mu * N * u2; }
  out[0] = A00 + N * (Bv0 + A10); out[1] = A01 + N * (Bv1 + A11); out[2] = A02 + N * (Bv2 + A12);
  const phid = (w0 * cp - w1 * sp) / ct;
  out[3] = phid; out[4] = w0 * sp + w1 * cp; out[5] = w2 - st * phid;
  out[6] = (N * n0 + F0) / m - g * n0 - (w1 * v2 - w2 * v1);
  out[7] = (N * n1 + F1) / m - g * n1 - (w2 * v0 - w0 * v2);
  out[8] = (N * n2 + F2) / m - g * n2 - (w0 * v1 - w1 * v0);
  out[9] = ct * cp * v0 - ct * sp * v1 + st * v2;
  out[10] = (cf * sp + sf * st * cp) * v0 + (cf * cp - sf * st * sp) * v1 - sf * ct * v2;
}
export function rhsInto(B, y, out) {
  if (B.law === 'rolling') { const f = rhs(B, y); for (let j = 0; j < 11; j++) out[j] = f[j]; }
  else rhsSlip(B, y, out);
}

// Integrate to tEnd and return samples every dtOut. The regularized Coulomb law is
// stiff (Jacobian ~ mu N / eps); the adaptive step simply shrinks to stay stable.
export function integrate(B, onProgress, shouldStop) {
  const nY = 11, nOut = Math.round(B.tEnd / B.dtOut) + 1;
  const T = new Float64Array(nOut), Y = Array.from({ length: nY }, () => new Float64Array(nOut));
  let y = Float64Array.from(initialState(B)), yNew = new Float64Array(nY);
  const ys = new Float64Array(nY), k = Array.from({ length: 7 }, () => new Float64Array(nY));
  rhsInto(B, y, k[0]);
  let t = 0, h = 1e-5, nSteps = 0, nRej = 0, iOut = 0, lastReport = 0;
  const store = () => { T[iOut] = t; for (let j = 0; j < nY; j++) Y[j][iOut] = y[j]; iOut++; };
  store();
  while (iOut < nOut) {
    const tNext = iOut * B.dtOut;
    if (t + h > tNext) h = tNext - t;
    for (let s = 1; s < 7; s++) {
      const As = A[s];
      for (let j = 0; j < nY; j++) { let acc = y[j]; for (let q = 0; q < s; q++) acc += h * As[q] * k[q][j]; ys[j] = acc; }
      rhsInto(B, ys, k[s]);
    }
    let err = 0;
    for (let j = 0; j < nY; j++) {
      let acc = y[j], e = 0;
      for (let q = 0; q < 6; q++) acc += h * B5[q] * k[q][j];
      for (let q = 0; q < 7; q++) e += h * E[q] * k[q][j];
      yNew[j] = acc;
      const sc = B.atol + B.rtol * Math.max(Math.abs(y[j]), Math.abs(acc));
      const r = Math.abs(e) / sc; if (r > err) err = r;
    }
    if (!(err < Infinity)) throw new Error('The integration diverged. Try a larger regularization velocity ε or a smaller friction coefficient.');
    if (err <= 1) {
      t += h; const tmp = y; y = yNew; yNew = tmp;
      const k6 = k[6]; k[6] = k[0]; k[0] = k6;           // FSAL: last stage is the next first stage
      nSteps++;
      if (Math.abs(t - tNext) < 1e-12) { t = tNext; store(); }
      if (t - lastReport > 0.05) {
        lastReport = t;
        if (onProgress) onProgress(t / B.tEnd);
        if (shouldStop && shouldStop()) throw new Error('Stopped.');
      }
    } else nRej++;
    const fac = err === 0 ? 5 : Math.min(5, Math.max(0.2, 0.9 * Math.pow(err, -0.2)));
    h = Math.max(h * fac, 1e-12);
    if (nSteps + nRej > 5e7) throw new Error('Too many steps. Shorten the simulated time or increase the regularization velocity ε.');
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

export function simulate(params, onProgress, shouldStop) {
  const B = makeBody({ ...DEFAULTS, ...params });
  const t0 = (typeof performance !== 'undefined' ? performance : Date).now();
  const { T, Y, nSteps, nRej } = integrate(B, onProgress, shouldStop);
  const D = diagnostics(B, T, Y);
  const S = summarize(B, T, Y, D);
  const ms = (typeof performance !== 'undefined' ? performance : Date).now() - t0;
  return { params: { ...DEFAULTS, ...params }, T, Y, D, S, stats: { nSteps, nRej, ms } };
}
