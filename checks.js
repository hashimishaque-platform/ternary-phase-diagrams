// Checks of the Flory-Huggins solver and of the fit to the published cloud points. Usage: node checks.js
'use strict';
const S = require('./solver.js');
const fs = require('fs');

// the nine measured cloud points (weight fractions), the only measured data in this project
const points = {};
for (const line of fs.readFileSync(__dirname + '/data/cloud_points.csv', 'utf8').split('\n')) {
  if (!line || line.startsWith('#') || line.startsWith('solvent,')) continue;
  const [s, ns, sv, p] = line.split(',');
  (points[s] = points[s] || []).push([+ns, +sv, +p]);
}

// fitted parameters, as on the page
const c13 = 1.53, n3 = 600;
const fits = { NMP: { c12: 1.10, c23: 0.35, rho: [0.781, 1.028, 1.40] }, DMF: { c12: 0.73, c23: 0.49, rho: [0.781, 0.944, 1.40] } };

let ok = true;
for (const s of ['NMP', 'DMF']) {
  const F = fits[s];
  const B = S.traceBinodal(n3, F.c12, c13, F.c23);
  if (!B) { console.log(s, ': no binodal'); ok = false; continue; }
  // 1. equal chemical potentials on every tie-line (per lattice site)
  let res = 0;
  for (const t of B.ties) {
    if (!t.x) continue;                       // the edge tie-line is solved separately
    const A = S.muFn(t.a[0], t.a[1], t.a[2], n3, F.c12, c13, F.c23, t.x[0]);
    const C = S.muFn(t.b[0], t.b[1], t.b[2], n3, F.c12, c13, F.c23, t.x[1], t.x[2]);
    res = Math.max(res, Math.abs(A[0] - C[0]), Math.abs(A[1] - C[1]), Math.abs(A[2] - C[2]) / n3);
  }
  // 2. the binodal ends on the spinodal (critical point): the spinodal determinant changes sign there
  const c = B.critical;
  const det = S.spinodalDet(c[0], c[1], c[2], n3, F.c12, c13, F.c23);
  const scale = (1 / c[0] + 1 / c[1]) * (1 / c[0] + 1 / (n3 * c[2]));
  // 3. cloud points against the measurements, in wt % of 2-propanol
  let sq = 0, worst = 0;
  const rows = [];
  for (const w of points[s]) {
    const v = S.wtToVol(w, F.rho);
    const rr = v[2] / (v[1] + v[2]);
    const t = S.cloudPoint(rr, B.curve);
    const m = S.volToWt([t, (1 - t) * (1 - rr), (1 - t) * rr], F.rho);
    const d = 100 * (m[0] - w[0]);
    sq += d * d; worst = Math.max(worst, Math.abs(d));
    rows.push(`    ${(100 * w[2]).toFixed(1).padStart(5)} wt % Torlon: measured ${(100 * w[0]).toFixed(1)}, model ${(100 * m[0]).toFixed(1)} wt % 2-propanol`);
  }
  const rms = Math.sqrt(sq / points[s].length);
  console.log(`Torlon / ${s} / 2-propanol: ${B.ties.length} tie-lines`);
  console.log(`  largest chemical potential difference on a tie-line: ${res.toExponential(1)} (should be below 1e-9)`);
  console.log(`  spinodal determinant at the critical point, relative: ${(det / scale).toExponential(1)} (should be near 0)`);
  console.log(`  critical point: ${(100 * c[2]).toFixed(2)} vol % polymer, ${(100 * c[0]).toFixed(2)} vol % non-solvent`);
  console.log(rows.join('\n'));
  console.log(`  RMS deviation ${rms.toFixed(2)} wt %, largest ${worst.toFixed(2)} wt %`);
  if (!(res < 1e-9 && Math.abs(det / scale) < 1e-2 && rms < 1.0)) ok = false;
}
// 4. CSV downloads: read the files back and check them against the solver
const readCSV = (text) => text.split('\n').filter((l) => l && !l.startsWith('#')).slice(1).map((l) => l.split(','));
const near = (a, b, tol) => Math.abs(a - b) <= tol;
for (const s of ['NMP', 'DMF']) {
  const F = fits[s];
  const B = S.traceBinodal(n3, F.c12, c13, F.c23);
  if (!B) continue;
  const info = { polymer: 'Torlon', solvent: s, nonsolvent: '2-propanol', params: { c12: F.c12, c13, c23: F.c23, n3 }, rho: F.rho };
  let bad = [];
  // a row of six percentages: both triples sum to 100, and the wt % follow from the vol % with the densities
  const rowOk = (r) => {
    const w = r.slice(0, 3).map(Number), v = r.slice(3, 6).map(Number);
    if (!(w.concat(v).every((x) => isFinite(x) && x >= 0))) return false;
    const ww = S.volToWt(v.map((x) => x / 100), F.rho).map((x) => 100 * x);
    return near(w[0] + w[1] + w[2], 100, 1e-7) && near(v[0] + v[1] + v[2], 100, 1e-7) && ww.every((x, i) => near(x, w[i], 1e-7));
  };
  // binodal: same points as the curve the page draws
  const bin = readCSV(S.binodalCSV(B, info));
  if (bin.length !== B.curve.length) bad.push('binodal row count');
  bin.forEach((r, i) => {
    if (!rowOk(r.slice(2))) bad.push(`binodal row ${i + 1}`);
    if (!B.curve[i].every((x, k) => near(100 * x, +r[5 + k], 1e-7))) bad.push(`binodal point ${i + 1}`);
  });
  if (bin[B.ties.length][1] !== 'critical point') bad.push('critical point row');
  // tie-lines: the two phases read back from the file still have equal chemical potentials
  const tie = readCSV(S.tieLinesCSV(B, info));
  if (tie.length !== B.ties.length) bad.push('tie-line row count');
  let dmu = 0;
  tie.forEach((r, i) => {
    if (!rowOk(r.slice(1, 7)) || !rowOk(r.slice(7, 13))) bad.push(`tie-line row ${i + 1}`);
    if (i === 0) return;                     // the edge tie-line has no solvent, ln 0 is not defined
    const a = r.slice(4, 7).map((x) => x / 100), b = r.slice(10, 13).map((x) => x / 100);
    const A = S.muFn(a[0], a[1], a[2], n3, F.c12, c13, F.c23), C = S.muFn(b[0], b[1], b[2], n3, F.c12, c13, F.c23);
    dmu = Math.max(dmu, Math.abs(A[0] - C[0]), Math.abs(A[1] - C[1]), Math.abs(A[2] - C[2]) / n3);
  });
  if (!(dmu < 1e-7)) bad.push(`tie-line chemical potentials differ by ${dmu.toExponential(1)}`);
  // spinodal: every point has a zero determinant, and the spinodal passes through the critical point
  const lines = S.spinodalLines(n3, F.c12, c13, F.c23);
  const sp = readCSV(S.spinodalCSV(lines, info));
  let worstDet = 0, dCrit = Infinity;
  sp.forEach((r, i) => {
    if (!rowOk(r.slice(2))) bad.push(`spinodal row ${i + 1}`);
    const v = r.slice(5, 8).map((x) => x / 100);
    const sc = (1 / v[0] + 1 / v[1]) * (1 / v[0] + 1 / (n3 * v[2]));
    worstDet = Math.max(worstDet, Math.abs(S.spinodalDet(v[0], v[1], v[2], n3, F.c12, c13, F.c23) / sc));
    dCrit = Math.min(dCrit, Math.hypot(v[0] - B.critical[0], v[1] - B.critical[1], v[2] - B.critical[2]));
  });
  if (!(worstDet < 1e-8)) bad.push(`spinodal determinant ${worstDet.toExponential(1)}`);
  if (!(dCrit < 5e-3)) bad.push(`spinodal misses the critical point by ${dCrit.toExponential(1)}`);
  console.log(`CSV files, Torlon / ${s}: ${bin.length} binodal points, ${tie.length} tie-lines, ${sp.length} spinodal points in ${lines.length} line(s)`);
  console.log(`  chemical potentials read back from the tie-line file agree to ${dmu.toExponential(1)} (should be below 1e-7)`);
  console.log(`  spinodal determinant at the exported points, relative: at most ${worstDet.toExponential(1)}; nearest point to the critical point ${dCrit.toExponential(1)}`);
  if (bad.length) { console.log('  problems: ' + bad.slice(0, 5).join('; ')); ok = false; }
}
console.log(ok ? 'All checks passed.' : 'Some checks failed.');
process.exit(ok ? 0 : 1);
