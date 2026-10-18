// Flory-Huggins ternary solver: binodal, spinodal, tie-lines and cloud points.
// Pure functions, no DOM. Used by the page (inlined by build.py) and by checks.js.

'use strict';

// Components: non-solvent (1), solvent (2), polymer (3).
// Volume fractions throughout. n1 = n2 = 1 (lattice sites per molecule), n3 = polymer.

function muFn(f1, f2, f3, n3, c12, c13, c23, lnf2, lnf3) {
  // Chemical potentials over RT, per molecule (Tompa / Altena-Smolders form, n1 = n2 = 1).
  // lnf2 and lnf3 may be passed directly: a polymer-lean phase can hold less polymer than
  // a double can represent (ln f3 < -745), and the logarithm must not be clamped.
  const L = (f) => Math.log(Math.max(f, 1e-300));
  const l2 = lnf2 !== undefined ? lnf2 : L(f2);
  const l3 = lnf3 !== undefined ? lnf3 : L(f3);
  const p1 = f1, p2 = f2, p3 = f3;
  const mu1 = L(p1) + 1 - p1 - p2 - p3 / n3
            + (c12 * p2 + c13 * p3) * (p2 + p3) - c23 * p2 * p3;
  const mu2 = l2 + 1 - p2 - p1 - p3 / n3
            + (c12 * p1 + c23 * p3) * (p1 + p3) - c13 * p1 * p3;
  const mu3 = l3 + 1 - p3 - n3 * p1 - n3 * p2
            + n3 * ((c13 * p1 + c23 * p2) * (p1 + p2) - c12 * p1 * p2);
  return [mu1, mu2, mu3];
}

// Spinodal: determinant of the Hessian of dG/RT in (f2, f3) with f1 = 1 - f2 - f3.
function spinodalDet(f1, f2, f3, n3, c12, c13, c23) {
  if (f1 <= 0 || f2 <= 0 || f3 <= 0) return NaN;
  const G22 = 1 / f1 + 1 / f2 - 2 * c12;
  const G33 = 1 / f1 + 1 / (n3 * f3) - 2 * c13;
  const G23 = 1 / f1 - c12 - c13 + c23;
  return G22 * G33 - G23 * G23;
}

function solve(A, b) {
  const n = b.length;
  const M = A.map((r, i) => r.concat([b[i]]));
  for (let i = 0; i < n; i++) {
    let piv = i;
    for (let k = i + 1; k < n; k++) if (Math.abs(M[k][i]) > Math.abs(M[piv][i])) piv = k;
    if (Math.abs(M[piv][i]) < 1e-300) return null;
    [M[i], M[piv]] = [M[piv], M[i]];
    for (let k = i + 1; k < n; k++) {
      const f = M[k][i] / M[i][i];
      for (let j = i; j <= n; j++) M[k][j] -= f * M[i][j];
    }
  }
  const x = new Array(n).fill(0);
  for (let i = n - 1; i >= 0; i--) {
    let s = M[i][n];
    for (let j = i + 1; j < n; j++) s -= M[i][j] * x[j];
    x[i] = s / M[i][i];
  }
  return x;
}

// Damped Newton with a numerical Jacobian. F returns null when x is outside the domain.
function newton(F, x0, tol, maxIter) {
  let x = x0.slice();
  let f = F(x);
  if (!f) return null;
  const norm = (v) => Math.max(...v.map(Math.abs));
  for (let it = 0; it < maxIter; it++) {
    const fn0 = norm(f);
    if (fn0 < tol) return x;
    const n = x.length;
    const J = [];
    for (let i = 0; i < n; i++) J.push(new Array(n).fill(0));
    for (let j = 0; j < n; j++) {
      const h = 1e-7 * Math.max(1, Math.abs(x[j]));
      const xp = x.slice(); xp[j] += h;
      const fp = F(xp);
      if (!fp) return null;
      for (let i = 0; i < n; i++) J[i][j] = (fp[i] - f[i]) / h;
    }
    const dx = solve(J, f.map((v) => -v));
    if (!dx) return null;
    let a = 1, ok = false;
    for (let k = 0; k < 25; k++) {
      const xn = x.map((v, i) => v + a * dx[i]);
      const fnew = F(xn);
      if (fnew && norm(fnew) < fn0) { x = xn; f = fnew; ok = true; break; }
      a *= 0.5;
    }
    if (!ok) return null;
  }
  return norm(f) < tol ? x : null;
}

// Binary demixing on the polymer / non-solvent edge (f2 = 0).
// Unknowns: [f3 rich, ln f3 lean]. The lean phase can hold less than 1e-300 polymer, so its
// logarithm is carried directly and exp() is only used where underflow to zero is harmless.
function binaryEdge(n3, c13) {
  const chiCrit = 0.5 * Math.pow(1 + 1 / Math.sqrt(n3), 2);
  if (c13 <= chiCrit + 1e-4) return null;
  const F = (x) => {
    const r = x[0], lnl = x[1];
    if (!(r > 0 && r < 1 && lnl < Math.log(r))) return null;
    const l = Math.exp(lnl);
    const m1r = Math.log(1 - r) + (1 - 1 / n3) * r + c13 * r * r;
    const m1l = Math.log1p(-l) + (1 - 1 / n3) * l + c13 * l * l;
    const m3r = Math.log(r) / n3 - (1 - 1 / n3) * (1 - r) + c13 * (1 - r) * (1 - r);
    const m3l = lnl / n3 - (1 - 1 / n3) * (1 - l) + c13 * (1 - l) * (1 - l);
    return [m1r - m1l, m3r - m3l];
  };
  for (const r0 of [0.95, 0.9, 0.8, 0.65, 0.5, 0.35, 0.2, 0.1]) {
    for (const l0 of [-3000, -1200, -500, -200, -60, -30, -15, -8, -4, -2, -1]) {
      const x = newton(F, [r0, l0], 1e-12, 80);
      if (x && x[0] - Math.exp(x[1]) > 1e-3) return { rich: x[0], lean: Math.exp(x[1]), lnLean: x[1] };
    }
  }
  return null;
}

// One tie-line at a fixed polymer fraction f3a in the polymer-rich phase.
// Unknowns: [ln f2a, ln f2b, ln f3b]. The logarithms are passed straight into muFn.
function tieAt(f3a, guess, n3, c12, c13, c23) {
  const F = (x) => {
    if (!(x[0] < 0 && x[1] < 0 && x[2] < Math.log(f3a))) return null;
    const f2a = Math.exp(x[0]), f2b = Math.exp(x[1]), f3b = Math.exp(x[2]);
    const f1a = 1 - f2a - f3a, f1b = 1 - f2b - f3b;
    if (!(f1a > 0 && f1b > 0)) return null;
    const A = muFn(f1a, f2a, f3a, n3, c12, c13, c23, x[0]);
    const B = muFn(f1b, f2b, f3b, n3, c12, c13, c23, x[1], x[2]);
    return [A[0] - B[0], A[1] - B[1], (A[2] - B[2]) / n3];
  };
  const x = newton(F, guess, 1e-10, 40);
  if (!x) return null;
  const f2a = Math.exp(x[0]), f2b = Math.exp(x[1]), f3b = Math.exp(x[2]);
  const a = [1 - f2a - f3a, f2a, f3a];
  const b = [1 - f2b - f3b, f2b, f3b];
  if (a[0] <= 0 || b[0] <= 0) return null;
  const len = Math.hypot(a[0] - b[0], a[1] - b[1], a[2] - b[2]);
  if (len < 2e-3) return null;          // collapsed onto a single phase
  return { a, b, x, len };
}

// Trace the binodal from the polymer / non-solvent edge to the critical point.
// Returns { ties, curve, critical } or null when no gap starts on that edge.
function traceBinodal(n3, c12, c13, c23) {
  const edge = binaryEdge(n3, c13);
  if (!edge) return null;
  const ties = [];
  let f3a = edge.rich;
  // start just inside the triangle
  let d = Math.min(0.01, 0.1 * (edge.rich - edge.lean));
  let guess = [Math.log(d), Math.log(d), edge.lnLean];
  let step = 0.01;
  let first = null;
  for (let k = 0; k < 6 && !first; k++) {
    first = tieAt(edge.rich - d, guess, n3, c12, c13, c23);
    if (!first) { d *= 0.3; guess = [Math.log(d), Math.log(d), guess[2]]; }
  }
  if (!first) return null;
  // edge tie (f2 = 0) for completeness
  ties.push({ a: [1 - edge.rich, 0, edge.rich], b: [1 - edge.lean, 0, edge.lean], len: edge.rich - edge.lean });
  ties.push(first);
  f3a = first.a[2];
  guess = first.x;
  let prevX = null;
  while (f3a > 1e-4 && ties.length < 400) {
    let next = null;
    while (step >= 1e-5) {
      // linear extrapolation of the unknowns improves the guess
      let g = guess;
      if (prevX) g = guess.map((v, i) => v + (v - prevX[i]) * Math.min(1, step / Math.max(1e-6, Math.abs(ties[ties.length - 1].a[2] - ties[ties.length - 2].a[2]))));
      next = tieAt(f3a - step, g, n3, c12, c13, c23) || tieAt(f3a - step, guess, n3, c12, c13, c23);
      if (next) break;
      step *= 0.5;
    }
    if (!next) break;
    prevX = guess;
    ties.push(next);
    f3a = next.a[2];
    guess = next.x;
    step = Math.min(step * 1.4, 0.03, 0.4 * next.len);  // finer steps as the tie-lines shorten
    if (next.len < 6e-3) break;  // essentially at the critical point
  }
  const last = ties[ties.length - 1];
  const critical = [0, 1, 2].map((i) => (last.a[i] + last.b[i]) / 2);
  // closed binodal: rich branch (edge -> critical), then lean branch (critical -> edge)
  const curve = ties.map((t) => t.a).concat([critical]).concat(ties.slice().reverse().map((t) => t.b));
  return { ties, curve, critical };
}

function pointInPolygon(pt, poly) {
  // pt, poly vertices in any consistent 2D projection [u, v]
  let inside = false;
  for (let i = 0, j = poly.length - 1; i < poly.length; j = i++) {
    const [xi, yi] = poly[i], [xj, yj] = poly[j];
    if ((yi > pt[1]) !== (yj > pt[1]) && pt[0] < ((xj - xi) * (pt[1] - yi)) / (yj - yi) + xi) inside = !inside;
  }
  return inside;
}

// Tie-line through a two-phase composition p, plus the lever rule (volume share of the
// polymer-rich phase). The bracketing tie-lines give a first estimate, which is then refined
// by solving for the tie-line that passes exactly through p. Falls back to interpolation
// if the refinement does not converge (very close to the critical point).
function tieThrough(p, ties, n3, c12, c13, c23) {
  const sideOf = (a, b) => (b[0] - a[0]) * (p[2] - a[2]) - (b[2] - a[2]) * (p[0] - a[0]);
  const lever = (a, b) => {
    const L = Math.hypot(a[0] - b[0], a[2] - b[2]);
    const fr = L > 0 ? Math.hypot(p[0] - b[0], p[2] - b[2]) / L : 0.5;
    return { a, b, fracRich: Math.max(0, Math.min(1, fr)) };
  };
  for (let i = 0; i < ties.length - 1; i++) {
    const s0 = sideOf(ties[i].a, ties[i].b), s1 = sideOf(ties[i + 1].a, ties[i + 1].b);
    if (!(s0 === 0 || s0 * s1 < 0)) continue;
    const w = s0 === 0 ? 0 : s0 / (s0 - s1);
    const interp = lever(
      [0, 1, 2].map((k) => ties[i].a[k] + w * (ties[i + 1].a[k] - ties[i].a[k])),
      [0, 1, 2].map((k) => ties[i].b[k] + w * (ties[i + 1].b[k] - ties[i].b[k])));
    if (n3 === undefined) return interp;
    // refine: bisection on the rich-phase polymer fraction between the two bracketing ties
    let hi = ties[i].a[2], lo = ties[i + 1].a[2], sHi = s0;
    let xHi = ties[i].x || ties[i + 1].x, xLo = ties[i + 1].x || ties[i].x;
    if (!xHi || !xLo) return interp;
    let best = null;
    for (let it = 0; it < 50; it++) {
      const m = 0.5 * (hi + lo);
      const guess = xHi.map((v, k) => 0.5 * (v + xLo[k]));
      const t = tieAt(m, guess, n3, c12, c13, c23) || tieAt(m, xLo, n3, c12, c13, c23) || tieAt(m, xHi, n3, c12, c13, c23);
      if (!t) return best || interp;
      best = lever(t.a, t.b);
      const sm = sideOf(t.a, t.b);
      if (Math.abs(sm) < 1e-13 || hi - lo < 1e-11) break;
      if ((sm < 0) === (sHi < 0)) { hi = m; sHi = sm; xHi = t.x; } else { lo = m; xLo = t.x; }
    }
    return best || interp;
  }
  return null;
}

// Minimum distance from composition p to the binodal curve (piecewise linear).
function distanceToCurve(p, curve) {
  let best = Infinity;
  for (let i = 0; i < curve.length - 1; i++) {
    const a = curve[i], b = curve[i + 1];
    const u = [b[0] - a[0], b[1] - a[1], b[2] - a[2]];
    const w = [p[0] - a[0], p[1] - a[1], p[2] - a[2]];
    const uu = u[0] * u[0] + u[1] * u[1] + u[2] * u[2];
    let t = uu > 0 ? (u[0] * w[0] + u[1] * w[1] + u[2] * w[2]) / uu : 0;
    t = Math.max(0, Math.min(1, t));
    const d = Math.hypot(w[0] - t * u[0], w[1] - t * u[1], w[2] - t * u[2]);
    if (d < best) best = d;
  }
  return best;
}

// Local stability: the Hessian must be positive definite (G22 > 0 and det > 0).
function locallyStable(f1, f2, f3, n3, c12, c13, c23) {
  if (f1 <= 0 || f2 <= 0 || f3 <= 0) return true;
  const G22 = 1 / f1 + 1 / f2 - 2 * c12;
  return G22 > 0 && spinodalDet(f1, f2, f3, n3, c12, c13, c23) > 0;
}

// Projection used for 2D geometry: [non-solvent, polymer].
function proj(p) { return [p[0], p[2]]; }

function inTwoPhase(p, curve) {
  return pointInPolygon(proj(p), curve.map(proj));
}

// Cloud point along a titration path: start from a binary solvent/polymer solution with
// polymer volume fraction r (on a non-solvent-free basis) and add non-solvent.
// Returns the non-solvent volume fraction at which the path first enters the two-phase region.
function cloudPoint(r, curve) {
  const at = (t) => [t, (1 - t) * (1 - r), (1 - t) * r];
  let prev = 0;
  for (let t = 0.0025; t < 0.999; t += 0.0025) {
    if (inTwoPhase(at(t), curve)) {
      let lo = prev, hi = t;
      for (let k = 0; k < 40; k++) {
        const m = 0.5 * (lo + hi);
        if (inTwoPhase(at(m), curve)) hi = m; else lo = m;
      }
      return 0.5 * (lo + hi);
    }
    prev = t;
  }
  return null;
}

// Volume fractions <-> weight fractions, given densities [rho_ns, rho_s, rho_p].
function volToWt(v, rho) {
  const m = v.map((x, i) => x * rho[i]); const s = m[0] + m[1] + m[2];
  return m.map((x) => x / s);
}
function wtToVol(w, rho) {
  const v = w.map((x, i) => x / rho[i]); const s = v[0] + v[1] + v[2];
  return v.map((x) => x / s);
}

// ---------- casting: dope, coagulation bath and the mixing line ----------
// Compositions are volume fractions [non-solvent, solvent, polymer]. A mixture of dope and bath
// lies on the straight line between them (mass balance with no volume change on mixing, as in
// the lattice model). The line is also straight in weight fractions, because the change from
// volume to weight fractions maps straight lines to straight lines.

// Fraction t of the way from p to q at which the line first enters the two-phase region
// (null if it never does). The entry is found by a scan and then by bisection.
function firstEntry(p, q, curve) {
  const at = (t) => [0, 1, 2].map((k) => p[k] + t * (q[k] - p[k]));
  if (inTwoPhase(p, curve)) return { t: 0, point: p.slice() };
  const N = 800;
  let prev = 0;
  for (let i = 1; i <= N; i++) {
    const t = i / N;
    if (inTwoPhase(at(t), curve)) {
      let lo = prev, hi = t;
      for (let k = 0; k < 50; k++) { const m = 0.5 * (lo + hi); if (inTwoPhase(at(m), curve)) hi = m; else lo = m; }
      const tt = 0.5 * (lo + hi);
      return { t: tt, point: at(tt) };
    }
    prev = t;
  }
  return null;
}

// Cloud point of the dope: add pure non-solvent until the solution turns turbid.
function dopeCloudPoint(dope, curve) {
  return firstEntry(dope, [1, 0, 0], curve);
}

// The most solvent a solvent / non-solvent bath can hold while the mixing line between dope and
// bath still reaches the two-phase region. Returns the solvent volume fraction of that bath,
// 0 if even pure non-solvent does not reach it, or 1 if every bath does.
function bathSolventLimit(dope, curve) {
  const reaches = (b) => firstEntry(dope, [1 - b, b, 0], curve) !== null;
  if (!reaches(0)) return 0;
  let lo = 0, hi = null;
  for (let b = 0.01; b < 1; b += 0.01) { if (!reaches(b)) { hi = b; break; } lo = b; }
  if (hi === null) return 1;
  for (let k = 0; k < 40; k++) { const m = 0.5 * (lo + hi); if (reaches(m)) lo = m; else hi = m; }
  return 0.5 * (lo + hi);
}

// ---------- export ----------

// Spinodal over the whole triangle, in volume fractions [non-solvent, solvent, polymer].
// The determinant is sampled on a triangular grid in (phi1, phi3); every grid edge where it changes
// sign is bisected on the exact determinant, so each point lies on the spinodal to rounding error.
// Points are joined into lines by the grid edge they lie on.
function spinodalLines(n3, c12, c13, c23, M) {
  M = M || 240;
  const eps = 1e-6;
  const node = (i, j) => {
    const v = [Math.max(i / M, eps), Math.max((M - i - j) / M, eps), Math.max(j / M, eps)];
    const t = v[0] + v[1] + v[2];
    return v.map((x) => x / t);
  };
  const det = (v) => spinodalDet(v[0], v[1], v[2], n3, c12, c13, c23);
  const val = new Map();
  const at = (i, j) => { const k = i * (M + 1) + j; if (!val.has(k)) val.set(k, det(node(i, j))); return val.get(k); };
  const crossing = new Map();
  const cross = (p, q) => {
    const kp = p[0] * (M + 1) + p[1], kq = q[0] * (M + 1) + q[1];
    const key = kp < kq ? kp + ':' + kq : kq + ':' + kp;
    if (crossing.has(key)) return key;
    let a = node(p[0], p[1]), b = node(q[0], q[1]);
    let da = at(p[0], p[1]);
    for (let it = 0; it < 60; it++) {
      const m = [0.5 * (a[0] + b[0]), 0.5 * (a[1] + b[1]), 0.5 * (a[2] + b[2])];
      const dm = det(m);
      if ((dm < 0) === (da < 0)) { a = m; da = dm; } else b = m;
    }
    crossing.set(key, [0, 1, 2].map((k) => 0.5 * (a[k] + b[k])));
    return key;
  };
  const segs = [];
  const tri = (A, B, C) => {
    const ks = [];
    for (const [p, q] of [[A, B], [B, C], [C, A]]) {
      const dp = at(p[0], p[1]), dq = at(q[0], q[1]);
      if (isFinite(dp) && isFinite(dq) && (dp < 0) !== (dq < 0)) ks.push(cross(p, q));
    }
    if (ks.length === 2) segs.push(ks);
  };
  for (let i = 0; i < M; i++) {
    for (let j = 0; j < M - i; j++) {
      tri([i, j], [i + 1, j], [i, j + 1]);
      if (i + j + 2 <= M) tri([i + 1, j], [i + 1, j + 1], [i, j + 1]);
    }
  }
  const ends = new Map();
  segs.forEach((sg, n) => { for (const e of [0, 1]) { if (!ends.has(sg[e])) ends.set(sg[e], []); ends.get(sg[e]).push([n, e]); } });
  const used = new Uint8Array(segs.length);
  const lines = [];
  for (let n = 0; n < segs.length; n++) {
    if (used[n]) continue;
    used[n] = 1;
    const fwd = [segs[n][1]], back = [segs[n][0]];
    for (const part of [fwd, back]) {
      for (;;) {
        const next = (ends.get(part[part.length - 1]) || []).find(([m]) => !used[m]);
        if (!next) break;
        used[next[0]] = 1;
        part.push(segs[next[0]][1 - next[1]]);
      }
    }
    lines.push(back.reverse().concat(fwd).map((k) => crossing.get(k)));
  }
  return lines;
}

// Number for a CSV cell: percent, ten significant digits, plain decimal where it is readable.
function csvNum(x) {
  const v = 100 * x;
  if (v === 0) return '0';
  return Math.abs(v) >= 1e-4 ? String(Number(v.toPrecision(10))) : v.toPrecision(10);
}

// Header lines that describe the system, the parameters and the units.
function csvHeader(what, info) {
  const P = info.params;
  return [
    `# ${what}`,
    `# System: ${info.polymer} (polymer) / ${info.solvent} (solvent) / ${info.nonsolvent} (non-solvent)`,
    `# Flory-Huggins parameters: chi12 = ${P.c12}, chi13 = ${P.c13}, chi23 = ${P.c23}, n3 = ${P.n3}`,
    `# Densities used for wt %, g/cm3: non-solvent ${info.rho[0]}, solvent ${info.rho[1]}, polymer ${info.rho[2]}`,
    '# All compositions in percent. vol % are the model volume fractions; wt % follow from them with the densities above.',
  ];
}

const comp = (v, rho) => { const w = volToWt(v, rho); return [w[0], w[1], w[2], v[0], v[1], v[2]].map(csvNum); };
const COLS = 'non-solvent wt %,solvent wt %,polymer wt %,non-solvent vol %,solvent vol %,polymer vol %';

// Binodal as one closed curve: polymer-rich branch from the polymer / non-solvent side to the
// critical point, then the polymer-lean branch back.
function binodalCSV(B, info) {
  const closed = B.ties[B.ties.length - 1].len < 0.03;
  const out = csvHeader('Binodal', info);
  if (!closed) out.push('# The calculation stopped before the critical point. The middle point is the centre of the last tie-line, not a critical point.');
  out.push('point,branch,' + COLS);
  const n = B.ties.length;
  B.curve.forEach((v, i) => {
    const branch = i < n ? 'polymer-rich' : i === n ? (closed ? 'critical point' : 'end') : 'polymer-lean';
    out.push([i + 1, branch].concat(comp(v, info.rho)).join(','));
  });
  return out.join('\n') + '\n';
}

// Tie-lines: the two coexisting phases on each line.
function tieLinesCSV(B, info) {
  const out = csvHeader('Tie-lines (all computed, polymer-rich phase first)', info);
  out.push('# The first tie-line lies on the polymer / non-solvent side, where there is no solvent.');
  const half = (p) => COLS.split(',').map((c) => `${p} ${c}`).join(',');
  out.push('tie-line,' + half('rich') + ',' + half('lean'));
  B.ties.forEach((t, i) => out.push([i + 1].concat(comp(t.a, info.rho), comp(t.b, info.rho)).join(',')));
  return out.join('\n') + '\n';
}

// Spinodal: one or more lines, each a sequence of points along the curve.
function spinodalCSV(lines, info) {
  const out = csvHeader('Spinodal', info);
  if (!lines.length) out.push('# No spinodal: the mixture is locally stable everywhere in the triangle.');
  out.push('line,point,' + COLS);
  lines.forEach((ln, k) => ln.forEach((v, i) => out.push([k + 1, i + 1].concat(comp(v, info.rho)).join(','))));
  return out.join('\n') + '\n';
}

if (typeof module !== 'undefined') {
  module.exports = { muFn, spinodalDet, solve, newton, binaryEdge, tieAt, traceBinodal, pointInPolygon, tieThrough, distanceToCurve, locallyStable, proj, inTwoPhase, cloudPoint, volToWt, wtToVol, spinodalLines, csvNum, firstEntry, dopeCloudPoint, bathSolventLimit, binodalCSV, tieLinesCSV, spinodalCSV };
}
