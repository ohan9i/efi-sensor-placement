// 유효독립성(EfI) 계산 — 논문의 4가지 대안
import { matmul, transpose, invSPD, logdetSPD, cholesky, invSqrtSPD, sub, rowsOf, zeros } from './linalg.js?v=1794145072';

export const fim = (S, W) => matmul(transpose(S), matmul(W, S)); // F = Sᵀ Σ⁻¹ S

// 정보 엔트로피 지표 −ln det F
export const entropy = (S, Sigma) => -logdetSPD(fim(S, invSPD(Sigma)));

// Alt.1  원래 정의 e_i = 1 − det(F₋ᵢ)/det(F)  (센서마다 Σ₋ᵢ 역행렬 필요)
export function efiDefinition(S, Sigma) {
  const n = S.length, ldF = logdetSPD(fim(S, invSPD(Sigma)));
  const e = new Float64Array(n);
  for (let i = 0; i < n; i++) {
    const idx = []; for (let j = 0; j < n; j++) if (j !== i) idx.push(j);
    const Fi = fim(rowsOf(S, idx), invSPD(sub(Sigma, idx)));
    const ld = logdetSPD(Fi);
    e[i] = ld === -Infinity ? 1 : 1 - Math.exp(ld - ldF);
  }
  return e;
}

// Alt.2  제안식 e = diag(Σ⁻¹SF⁻¹SᵀΣ⁻¹) ⊘ diag(Σ⁻¹)  (역행렬 1회)
export function efiProposed(S, Sigma) {
  const W = invSPD(Sigma), G = matmul(W, S), Fi = invSPD(matmul(transpose(S), G));
  return quadDiag(G, Fi).map((v, i) => v / W[i][i]);
}

// 기존 근사식 e_eq = diag(LᵀSF⁻¹SᵀL),  Σ⁻¹ = LLᵀ
export function efiFactor(S, Sigma, L) {
  const W = invSPD(Sigma), Fi = invSPD(fim(S, W));
  return quadDiag(matmul(transpose(L), S), Fi);
}
export const factorSqrt = (Sigma) => invSqrtSPD(Sigma); // Alt.3  L = Σ^(-1/2)
export const factorChol = (Sigma) => cholesky(invSPD(Sigma)); // Alt.4  Σ⁻¹ = L_c L_cᵀ

// diag(G F⁻¹ Gᵀ)
function quadDiag(G, Fi) {
  const H = matmul(G, Fi);
  return Float64Array.from(G, (g, i) => g.reduce((s, a, k) => s + a * H[i][k], 0));
}

export const METHODS = {
  alt1: { name: '원래 정의', tag: 'Alt.1', note: '센서마다 역행렬, 정확', color: '#7d93cf', f: efiDefinition },
  alt2: { name: '제안식 Eq.(26)', tag: 'Alt.2', note: '역행렬 한 번, 정확', color: '#2747a8', f: efiProposed },
  alt3: { name: '근사 L = Σ^-1/2', tag: 'Alt.3', note: '행렬 제곱근 분해', color: '#9a9286', f: (S, Sg) => efiFactor(S, Sg, factorSqrt(Sg)) },
  alt4: { name: '근사 L = 촐레스키', tag: 'Alt.4', note: '하삼각 분해', color: '#a0782c', f: (S, Sg) => efiFactor(S, Sg, factorChol(Sg)) },
};

// 동점은 낮은 인덱스 우선 (대칭 구조에서 방법 간 비교를 공정하게)
export function argminTol(e, tol = 1e-9) {
  let m = Infinity; for (const v of e) if (v < m) m = v;
  for (let i = 0; i < e.length; i++) if (e[i] <= m + tol) return i;
  return 0;
}

// 공분산 행렬 Σ (Eq. 28 / 29)
export function buildSigma(x, Phi, { lambda = 4, model = 'eq29' } = {}) {
  const n = x.length, Sg = zeros(n, n), m = Phi[0].length;
  for (let i = 0; i < n; i++) for (let j = i; j < n; j++) {
    let rho;
    if (model === 'identity') rho = i === j ? 1 : 0;
    else {
      rho = Math.exp(-Math.abs(x[i] - x[j]) / lambda);
      if (model === 'eq29') {
        let s = 0;
        for (let k = 0; k < m; k++) {
          const a = Math.abs(Phi[i][k]), b = Math.abs(Phi[j][k]), mx = Math.max(a, b);
          s += mx < 1e-12 ? 1 : Math.min(a, b) / mx;
        }
        rho *= s / m;
      }
    }
    Sg[i][j] = rho; Sg[j][i] = rho;
  }
  return Sg;
}

// 순차 제거 (backward sequential deletion)
export function sequentialDeletion(Phi, SigmaFull, methodKey, target) {
  const f = METHODS[methodKey].f;
  let active = Phi.map((_, i) => i);
  const removed = [], entropyByCount = {}, t0 = performance.now();
  entropyByCount[active.length] = entropy(rowsOf(Phi, active), sub(SigmaFull, active));
  while (active.length > target) {
    const S = rowsOf(Phi, active), Sg = sub(SigmaFull, active);
    const e = f(S, Sg);
    const k = argminTol(e);
    removed.push(active[k]);
    active = active.filter((_, q) => q !== k);
    entropyByCount[active.length] = entropy(rowsOf(Phi, active), sub(SigmaFull, active));
  }
  return { removed, entropyByCount, time: performance.now() - t0, n0: Phi.length };
}

// 제거 순서로부터 센서 수 n일 때의 활성 센서 집합
export function activeAt(result, n) {
  const gone = new Set(result.removed.slice(0, result.n0 - n));
  const a = []; for (let i = 0; i < result.n0; i++) if (!gone.has(i)) a.push(i);
  return a;
}
