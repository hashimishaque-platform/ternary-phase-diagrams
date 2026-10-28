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
// 5. Casting: the mixing line between a dope and a coagulation bath
for (const s of ['NMP', 'DMF']) {
  const F = fits[s];
  const B = S.traceBinodal(n3, F.c12, c13, F.c23);
  if (!B) continue;
  const bad = [];
  const ns = [1, 0, 0];
  // (a) for a dope with no non-solvent, the dope cloud point equals the titration cloud point used in the fit
  let dc = 0;
  for (const w of points[s]) {
    const v = S.wtToVol(w, F.rho), rr = v[2] / (v[1] + v[2]);
    const t = S.cloudPoint(rr, B.curve), c = S.dopeCloudPoint([0, 1 - rr, rr], B.curve);
    dc = Math.max(dc, Math.abs(c.point[0] - t));
  }
  if (!(dc < 1e-6)) bad.push(`dope cloud point differs from the titration by ${dc.toExponential(1)}`);
  // (b) the entry point of a mixing line lies on the binodal, and the line is straight in wt % as well
  const dope = S.wtToVol([0, 0.85, 0.15], F.rho);
  const bath = S.wtToVol([0.7, 0.3, 0], F.rho);
  const e = S.firstEntry(dope, bath, B.curve);
  if (!e) bad.push('mixing line does not enter'); else {
    const dist = S.distanceToCurve(e.point, B.curve);
    if (!(dist < 1e-9)) bad.push(`entry point is ${dist.toExponential(1)} off the binodal`);
    const [a, b, q] = [dope, bath, e.point].map((v) => S.volToWt(v, F.rho));
    const cross = (b[0] - a[0]) * (q[2] - a[2]) - (b[2] - a[2]) * (q[0] - a[0]);
    if (!(Math.abs(cross) < 1e-12)) bad.push(`mixing line is not straight in wt %: ${cross.toExponential(1)}`);
  }
  // (c) the bath limit: a little less solvent still reaches the two-phase region, a little more does not
  const L = S.bathSolventLimit(dope, B.curve);
  const reach = (b) => S.firstEntry(dope, [1 - b, b, 0], B.curve) !== null;
  if (!(L > 0.01 && L < 0.99 && reach(L - 0.005) && !reach(L + 0.005))) bad.push(`bath limit ${L}`);
  // (d) pure non-solvent always reaches the two-phase region from a clear dope
  if (!S.firstEntry(dope, ns, B.curve)) bad.push('pure non-solvent bath does not reach');
  const Lw = S.volToWt([1 - L, L, 0], F.rho)[1];
  console.log(`Casting, Torlon / ${s}, dope 15 wt % Torlon: bath limit ${(100 * Lw).toFixed(1)} wt % ${s}; dope cloud point matches the titration to ${dc.toExponential(1)}`);
  if (bad.length) { console.log('  problems: ' + bad.join('; ')); ok = false; }
}
// 6. Recipe bound and porosity formula
for (const s of ['NMP', 'DMF']) {
  const F = fits[s];
  const B = S.traceBinodal(n3, F.c12, c13, F.c23);
  if (!B) continue;
  const bad = [];
  for (const bs of [0, 0.3]) {
    const bath = S.wtToVol([1 - bs, bs, 0], F.rho);
    const pmin = S.minPolymerAboveCritical(0, bath, B.curve, B.critical);
    if (pmin === null) { bad.push(`no lower polymer bound for bath ${bs}`); continue; }
    // just above the bound the entry is above the critical point, just below it is not
    const eAt = (p) => S.firstEntry([0, 1 - p, p], bath, B.curve);
    const hi = eAt(pmin + 1e-6), lo = eAt(pmin - 1e-4);
    if (!(hi && hi.point[2] > B.critical[2])) bad.push('entry above the bound is not above the critical point');
    if (lo && lo.point[2] > B.critical[2]) bad.push('entry below the bound is already above the critical point');
    console.log(`Recipe, Torlon / ${s}, bath ${100 * bs} wt % ${s}: polymer at least ${(100 * S.volToWt([0, 1 - pmin, pmin], F.rho)[2]).toFixed(1)} wt % for entry above the critical point`);
  }
  if (bad.length) { console.log('  problems: ' + bad.join('; ')); ok = false; }
}
{
  // hand calculation: 1.0 g wet, 0.4 g dry, kerosene 0.82, PVDF 1.78 g/cm3 (densities from Ali et al. 2018)
  const e = S.gravimetricPorosity(1.0, 0.4, 0.82, 1.78);
  const hand = (0.6 / 0.82) / (0.6 / 0.82 + 0.4 / 1.78);
  const ok6 = Math.abs(e - hand) < 1e-12 && Math.abs(e - 0.76505) < 1e-4 && S.gravimetricPorosity(0.4, 0.5, 1, 1) === null;
  console.log(`Gravimetric porosity of the worked example: ${(100 * e).toFixed(2)} % (by hand 76.50 %)`);
  if (!ok6) ok = false;
}
// 7. Hansen distance, RED and chi (polymer films on wafers)
{
  const H = {};
  for (const line of fs.readFileSync(__dirname + '/data/hansen_parameters.csv', 'utf8').split('\n')) {
    if (!line || line.startsWith('#') || line.startsWith('name,')) continue;
    const [name, kind, dD, dP, dH, R0, M, rho] = line.split(',');
    H[name] = { kind, d: [+dD, +dP, +dH], R0: +R0, V: +M / +rho };
  }
  const bad = [];
  const P = H.PMMA, et = H.ethanol;
  // by hand: 4 (18.64 - 15.8)^2 + (10.52 - 8.8)^2 + (7.51 - 19.4)^2 = 32.2624 + 2.9584 + 141.3721 = 176.5929
  const ra = S.hansenRa(P.d, et.d), handRa = Math.sqrt(176.5929);
  if (Math.abs(ra - handRa) > 1e-12) bad.push(`Ra ${ra} against ${handRa} by hand`);
  const red = S.hansenRED(P.d, et.d, P.R0);
  if (Math.abs(red - handRa / 8.59) > 1e-12) bad.push('RED');
  // chi by hand: V = 46.07 / 0.789 cm3/mol, V Ra^2 / (4 R T)
  const handChi = (46.07 / 0.789) * 176.5929 / (4 * 8.314462618 * 298.15);
  const chi = S.hansenChi(P.d, et.d, et.V);
  if (Math.abs(chi - handChi) > 1e-12) bad.push('chi');
  // the 2-propanol row of chi12 on the page was estimated this way: reproduce it to within 0.01
  const page = { NMP: 1.10, DMF: 0.73, DMAc: 0.55, DMSO: 1.33, THF: 0.58, acetone: 0.83 };
  let worst12 = 0;
  for (const s in page) worst12 = Math.max(worst12, Math.abs(S.hansenChi(H['2-propanol'].d, H[s].d, H['2-propanol'].V) - page[s]));
  if (!(worst12 < 0.01)) bad.push(`2-propanol chi12 row off by ${worst12.toFixed(3)}`);
  // symmetry, and zero distance for identical points
  if (S.hansenRa(P.d, et.d) !== S.hansenRa(et.d, P.d) || S.hansenRa(P.d, P.d) !== 0) bad.push('distance not symmetric');
  // every liquid has a molar volume between 10 and 200 cm3/mol
  for (const n in H) if (H[n].kind === 'liquid' && !(H[n].V > 10 && H[n].V < 200)) bad.push(`molar volume of ${n}`);
  console.log(`Hansen: PMMA / ethanol Ra ${ra.toFixed(2)} MPa^1/2, RED ${red.toFixed(3)}, chi ${chi.toFixed(3)}; 2-propanol chi12 row reproduced to ${worst12.toFixed(3)}`);
  if (bad.length) { console.log('  problems: ' + bad.join('; ')); ok = false; }
}
console.log(ok ? 'All checks passed.' : 'Some checks failed.');
process.exit(ok ? 0 : 1);
