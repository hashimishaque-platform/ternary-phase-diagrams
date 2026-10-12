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
console.log(ok ? 'All checks passed.' : 'Some checks failed.');
process.exit(ok ? 0 : 1);
