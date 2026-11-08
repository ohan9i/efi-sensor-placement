import { fromArray, invSPD } from '../js/linalg.js';
import { efiDefinition, efiProposed, efiFactor, factorSqrt, METHODS, buildSigma, sequentialDeletion, activeAt } from '../js/efi.js';
import { buildBeam, selectModes } from '../js/beam.js';
const f = (a) => Array.from(a, (v) => v.toFixed(3)).join(', ');
// Example 1
let S = fromArray([[5, 3], [1, 0], [0, 1]]), Sg = fromArray([[1, 0, 0], [0, 1, 0], [0, 0, 1]]);
console.log('ex1 def', f(efiDefinition(S, Sg)), '| perm', f(efiFactor(S, Sg, fromArray([[0, 1, 0], [1, 0, 0], [0, 0, 1]]))), '| prop', f(efiProposed(S, Sg)));
S = fromArray([[4, 5], [1, 0], [0, -1]]); Sg = fromArray([[1, 0, 0.8], [0, 1, 0], [0.8, 0, 1]]);
console.log('ex2 def', f(efiDefinition(S, Sg)), '| sqrt', f(efiFactor(S, Sg, factorSqrt(Sg))), '| prop', f(efiProposed(S, Sg)));
console.log('L', factorSqrt(Sg).map((r) => f(r)));
let t = performance.now();
const beam = buildBeam();
const { Phi, freqRatio } = selectModes(beam);
console.log('beam', (performance.now() - t).toFixed(0), 'ms; modes', Phi[0].length, freqRatio.map((v) => v.toFixed(2)).join(' '));
const Sigma = buildSigma(beam.x, Phi, { lambda: 4, model: 'eq29' });
for (const k of ['alt1', 'alt2', 'alt3', 'alt4']) {
  const r = sequentialDeletion(Phi, Sigma, k, 10);
  const a20 = activeAt(r, 20);
  console.log(k, r.time.toFixed(0) + 'ms', 'H10=', r.entropyByCount[10].toFixed(2), 'H20=', r.entropyByCount[20].toFixed(2), 'H40=', r.entropyByCount[40].toFixed(2), '\n  20:', a20.map((i) => beam.x[i]).join(' '));
}
