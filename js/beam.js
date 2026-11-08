// 연속보 유한요소 모델 (논문 4장, Eq. 27) — 휨 변형만 고려한 오일러-베르누이 보 요소
import { zeros, cholesky, invLower, matmul, transpose, symEig } from './linalg.js';

export function buildBeam(spans = [8, 12, 12, 8], le = 0.5) {
  const nEl = spans.reduce((a, s) => a + Math.round(s / le), 0);
  const nNode = nEl + 1;
  const x = Array.from({ length: nNode }, (_, i) => i * le);
  const supports = [0];
  let acc = 0;
  for (const s of spans) { acc += Math.round(s / le); supports.push(acc); }
  const l = le, l2 = l * l;
  // 상수 ρAL/420, EI/L³ 생략 (논문과 동일)
  const m = [[156, 22 * l, 54, -13 * l], [22 * l, 4 * l2, 13 * l, -3 * l2], [54, 13 * l, 156, -22 * l], [-13 * l, -3 * l2, -22 * l, 4 * l2]];
  const k = [[12, 6 * l, -12, 6 * l], [6 * l, 4 * l2, -6 * l, 2 * l2], [-12, -6 * l, 12, -6 * l], [6 * l, 2 * l2, -6 * l, 4 * l2]];
  const nd = 2 * nNode, M = zeros(nd, nd), K = zeros(nd, nd);
  for (let e = 0; e < nEl; e++) {
    const d = [2 * e, 2 * e + 1, 2 * e + 2, 2 * e + 3];
    for (let a = 0; a < 4; a++) for (let b = 0; b < 4; b++) { M[d[a]][d[b]] += m[a][b]; K[d[a]][d[b]] += k[a][b]; }
  }
  // 경계조건: 지점에서 수직변위 구속 (핀 1개 + 롤러 4개)
  const fixed = new Set(supports.map((n) => 2 * n));
  const free = []; for (let i = 0; i < nd; i++) if (!fixed.has(i)) free.push(i);
  const Mf = free.map((i) => Float64Array.from(free, (j) => M[i][j]));
  const Kf = free.map((i) => Float64Array.from(free, (j) => K[i][j]));
  // 일반 고유치 문제 Kφ = ω²Mφ → M = RRᵀ, C = R⁻¹KR⁻ᵀ
  const R = cholesky(Mf), Ri = invLower(R);
  const C = matmul(matmul(Ri, Kf), transpose(Ri));
  for (let i = 0; i < C.length; i++) for (let j = 0; j < i; j++) { const v = 0.5 * (C[i][j] + C[j][i]); C[i][j] = v; C[j][i] = v; }
  const { values, vectors } = symEig(C, 1e-13);
  const order = values.map((v, i) => i).sort((a, b) => values[a] - values[b]);
  const RiT = transpose(Ri);
  const modes = order.map((idx) => {
    const v = vectors.map((r) => r[idx]);
    const phiF = RiT.map((r) => r.reduce((s, a, q) => s + a * v[q], 0)); // φ = R⁻ᵀv (질량 정규화)
    const full = new Float64Array(nd); free.forEach((dof, q) => (full[dof] = phiF[q]));
    const w = Array.from({ length: nNode }, (_, n) => full[2 * n]);
    // 부호 정규화
    let mi = 0; w.forEach((a, i) => { if (Math.abs(a) > Math.abs(w[mi])) mi = i; });
    const sg = w[mi] < 0 ? -1 : 1;
    return { omega2: values[idx], w: w.map((a) => a * sg) };
  });
  return { spans, le, nNode, x, supports, modes, L: acc * le };
}

// 고유진동수가 1차의 ratio배 이하인 모드 선택 → Φ (nNode × m)
export function selectModes(beam, ratio = 10) {
  const w1 = Math.sqrt(beam.modes[0].omega2);
  const sel = beam.modes.filter((md) => Math.sqrt(md.omega2) <= ratio * w1 * (1 + 1e-9));
  const Phi = Array.from({ length: beam.nNode }, (_, i) => Float64Array.from(sel, (md) => md.w[i]));
  return { Phi, modes: sel, freqRatio: sel.map((md) => Math.sqrt(md.omega2) / w1) };
}
