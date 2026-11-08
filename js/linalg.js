// 작은 밀집(dense) 선형대수 라이브러리 — 행렬은 number[][] (행 우선)

export const zeros = (r, c) => Array.from({ length: r }, () => new Float64Array(c));
export const eye = (n) => { const I = zeros(n, n); for (let i = 0; i < n; i++) I[i][i] = 1; return I; };
export const clone = (A) => A.map((r) => Float64Array.from(r));
export const fromArray = (a) => a.map((r) => Float64Array.from(r));

export function transpose(A) {
  const r = A.length, c = A[0].length, T = zeros(c, r);
  for (let i = 0; i < r; i++) for (let j = 0; j < c; j++) T[j][i] = A[i][j];
  return T;
}

export function matmul(A, B) {
  const n = A.length, m = B[0].length, k = B.length, C = zeros(n, m);
  for (let i = 0; i < n; i++) {
    const Ai = A[i], Ci = C[i];
    for (let p = 0; p < k; p++) {
      const a = Ai[p]; if (a === 0) continue;
      const Bp = B[p];
      for (let j = 0; j < m; j++) Ci[j] += a * Bp[j];
    }
  }
  return C;
}

// 부분행렬 A[idx, idx]
export function sub(A, rows, cols = rows) {
  return rows.map((i) => { const r = new Float64Array(cols.length); cols.forEach((j, q) => (r[q] = A[i][j])); return r; });
}
export const rowsOf = (A, rows) => rows.map((i) => Float64Array.from(A[i]));

// 촐레스키: A = L Lᵀ (L 하삼각). 실패 시 null
export function cholesky(A) {
  const n = A.length, L = zeros(n, n);
  for (let j = 0; j < n; j++) {
    let s = A[j][j];
    const Lj = L[j];
    for (let k = 0; k < j; k++) s -= Lj[k] * Lj[k];
    if (!(s > 0)) return null;
    const d = Math.sqrt(s); Lj[j] = d;
    for (let i = j + 1; i < n; i++) {
      const Li = L[i]; let t = A[i][j];
      for (let k = 0; k < j; k++) t -= Li[k] * Lj[k];
      Li[j] = t / d;
    }
  }
  return L;
}

export function logdetSPD(A) {
  const L = cholesky(A); if (!L) return -Infinity;
  let s = 0; for (let i = 0; i < A.length; i++) s += Math.log(L[i][i]);
  return 2 * s;
}

// 하삼각 역행렬
export function invLower(L) {
  const n = L.length, X = zeros(n, n);
  for (let j = 0; j < n; j++) {
    X[j][j] = 1 / L[j][j];
    for (let i = j + 1; i < n; i++) {
      let s = 0; const Li = L[i];
      for (let k = j; k < i; k++) s += Li[k] * X[k][j];
      X[i][j] = -s / Li[i];
    }
  }
  return X;
}

// SPD 역행렬 (촐레스키 기반)
export function invSPD(A) {
  const L = cholesky(A); if (!L) throw new Error('행렬이 양의 정부호가 아닙니다');
  const Li = invLower(L), n = A.length, X = zeros(n, n);
  // A⁻¹ = L⁻ᵀ L⁻¹
  for (let i = 0; i < n; i++) for (let j = 0; j <= i; j++) {
    let s = 0; for (let k = i; k < n; k++) s += Li[k][i] * Li[k][j];
    X[i][j] = s; X[j][i] = s;
  }
  return X;
}

// 대칭 행렬 고유분해 (순환 Jacobi). 반환 {values, vectors(열=고유벡터)}
export function symEig(Ain, tol = 1e-12, maxSweeps = 60) {
  const n = Ain.length, A = clone(Ain), V = eye(n);
  for (let sweep = 0; sweep < maxSweeps; sweep++) {
    let off = 0, diag = 0;
    for (let i = 0; i < n; i++) { diag += A[i][i] * A[i][i]; for (let j = i + 1; j < n; j++) off += A[i][j] * A[i][j]; }
    if (off <= tol * tol * diag) break;
    for (let p = 0; p < n - 1; p++) for (let q = p + 1; q < n; q++) {
      const apq = A[p][q];
      if (Math.abs(apq) < 1e-300) continue;
      const theta = (A[q][q] - A[p][p]) / (2 * apq);
      const t = Math.sign(theta || 1) / (Math.abs(theta) + Math.sqrt(theta * theta + 1));
      const c = 1 / Math.sqrt(t * t + 1), s = t * c;
      for (let k = 0; k < n; k++) { // 열 회전
        const akp = A[k][p], akq = A[k][q];
        A[k][p] = c * akp - s * akq; A[k][q] = s * akp + c * akq;
      }
      const Ap = A[p], Aq = A[q];
      for (let k = 0; k < n; k++) { // 행 회전
        const apk = Ap[k], aqk = Aq[k];
        Ap[k] = c * apk - s * aqk; Aq[k] = s * apk + c * aqk;
      }
      for (let k = 0; k < n; k++) {
        const vkp = V[k][p], vkq = V[k][q];
        V[k][p] = c * vkp - s * vkq; V[k][q] = s * vkp + c * vkq;
      }
    }
  }
  return { values: Array.from({ length: n }, (_, i) => A[i][i]), vectors: V };
}

// SPD 행렬의 역제곱근 A^(-1/2) — 고유분해 사용
export function invSqrtSPD(A) {
  const { values, vectors: V } = symEig(A);
  const n = A.length, X = zeros(n, n);
  for (let k = 0; k < n; k++) {
    const w = 1 / Math.sqrt(Math.max(values[k], 1e-300));
    for (let i = 0; i < n; i++) { const a = V[i][k] * w; if (a === 0) continue; for (let j = 0; j < n; j++) X[i][j] += a * V[j][k]; }
  }
  return X;
}
