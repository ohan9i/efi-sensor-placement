import { fromArray, invSPD, invSqrtSPD, cholesky, matmul, transpose, rowsOf, sub, zeros } from './linalg.js';
import { efiDefinition, efiProposed, efiFactor, METHODS, buildSigma, sequentialDeletion, activeAt, argminTol, entropy } from './efi.js';
import { buildBeam, selectModes } from './beam.js';

const $ = (id) => document.getElementById(id);
const NS = 'http://www.w3.org/2000/svg';
function el(tag, attrs = {}, parent, text) {
  const e = document.createElementNS(NS, tag);
  for (const [k, v] of Object.entries(attrs)) e.setAttribute(k, v);
  if (text != null) e.textContent = text;
  if (parent) parent.appendChild(e);
  return e;
}
const clear = (s) => { while (s.firstChild) s.removeChild(s.firstChild); };
const fmt = (v, d = 3) => (Math.abs(v) < 5e-12 ? 0 : v).toFixed(d);
const svgW = (s) => s.getBoundingClientRect().width || 800;
const KEYS = ['alt1', 'alt2', 'alt3', 'alt4'];

// ───────────────────────── ② 3센서 예제 ─────────────────────────
const toy = { S: [[5, 3], [1, 0], [0, 1]], v: [1, 1, 1], r: { r12: 0, r13: 0, r23: 0 }, factor: 'perm' };
const PRESETS = {
  ex1: { S: [[5, 3], [1, 0], [0, 1]], v: [1, 1, 1], r: { r12: 0, r13: 0, r23: 0 }, factor: 'perm' },
  ex2: { S: [[4, 5], [1, 0], [0, -1]], v: [1, 1, 1], r: { r12: 0, r13: 0.8, r23: 0 }, factor: 'sqrt' },
};
const SENSOR_COL = ['#8250df', '#bf3989', '#0a7ea4'];

function buildToyInputs() {
  const t = $('Sinput'); t.innerHTML = '<tr><td></td><td>θ₁</td><td>θ₂</td></tr>';
  toy.S.forEach((row, i) => {
    const tr = document.createElement('tr');
    tr.innerHTML = `<td style="color:${SENSOR_COL[i]}">센서 ${i + 1}</td>`;
    row.forEach((v, j) => {
      const td = document.createElement('td'), inp = document.createElement('input');
      inp.type = 'number'; inp.step = '0.5'; inp.value = v;
      inp.addEventListener('input', () => { const x = parseFloat(inp.value); if (!isNaN(x)) { toy.S[i][j] = x; renderToy(); } });
      td.appendChild(inp); tr.appendChild(td);
    });
    t.appendChild(tr);
  });
}
function syncToyControls() {
  ['v1', 'v2', 'v3'].forEach((id, i) => { $(id).value = toy.v[i]; });
  ['r12', 'r13', 'r23'].forEach((id) => { $(id).value = toy.r[id]; });
  $('factor').value = toy.factor;
  document.querySelectorAll('#toy .controls label').forEach((lb) => { const inp = lb.querySelector('input[type=range]'); if (inp) lb.querySelector('span').textContent = (+inp.value).toFixed(2); });
}
function toySigma() {
  const { v, r } = toy, R = [[1, r.r12, r.r13], [r.r12, 1, r.r23], [r.r13, r.r23, 1]];
  return fromArray(R.map((row, i) => row.map((c, j) => c * Math.sqrt(v[i] * v[j]))));
}
function toyFactor(Sg) {
  if (toy.factor === 'chol') return cholesky(invSPD(Sg));
  const R = invSqrtSPD(Sg);
  if (toy.factor === 'perm') return matmul(R, fromArray([[0, 1, 0], [1, 0, 0], [0, 0, 1]]));
  return R;
}
const texMat = (A, d = 2) => '\\begin{bmatrix}' + Array.from(A, (r) => Array.from(r, (v) => fmt(v, d)).join(' & ')).join('\\\\') + '\\end{bmatrix}';
function tex(str, display = false) { return window.katex ? katex.renderToString(str, { displayMode: display, throwOnError: false }) : str; }

function renderToy() {
  const S = fromArray(toy.S), Sg = toySigma();
  const res = $('toyResult'), msg = $('toyMsg');
  if (!cholesky(Sg)) { res.innerHTML = ''; msg.className = 'msg'; msg.textContent = '이 상관계수 조합은 유효한 공분산(양의 정부호)이 아닙니다. 값을 조정해 주세요.'; return; }
  let eDef;
  try { eDef = efiDefinition(S, Sg); } catch { res.innerHTML = ''; msg.className = 'msg'; msg.textContent = 'F가 특이합니다 (S의 열이 독립이 아님).'; return; }
  const L = toyFactor(Sg);
  const W = invSPD(Sg);
  const rows = [
    { name: '정의식 e = 1 − det F₋ᵢ / det F', e: eDef, ref: true },
    { name: '제안식 Eq.(26)', e: efiProposed(S, Sg) },
    { name: '기존 근사 Eq.(18) diag(LᵀSF⁻¹SᵀL)', e: efiFactor(S, Sg, L) },
    { name: '기본 EfI (Σ 무시) diag(S(SᵀS)⁻¹Sᵀ)', e: efiProposed(S, fromArray([[1, 0, 0], [0, 1, 0], [0, 0, 1]])) },
  ];
  const kRef = argminTol(eDef, 1e-9);
  let html = '<table class="restbl"><tr><th>계산법</th><th>e₁</th><th>e₂</th><th>e₃</th><th>제거</th></tr>';
  let approxWrong = false;
  rows.forEach((r, idx) => {
    const k = argminTol(r.e, 1e-9);
    const same = r.e.every((v, i) => Math.abs(v - eDef[i]) < 1e-6);
    if (idx === 2 && k !== kRef) approxWrong = true;
    html += `<tr class="${r.ref ? '' : same ? 'right' : 'wrong'}"><td style="text-align:left">${r.name}</td>` +
      Array.from(r.e, (v, i) => `<td class="${i === k ? 'min' : ''}">${fmt(v)}</td>`).join('') + `<td>센서 ${k + 1}</td></tr>`;
  });
  res.innerHTML = html + '</table>';
  msg.className = approxWrong ? 'msg' : 'msg ok';
  msg.innerHTML = approxWrong
    ? `⚠️ 기존 근사식은 <b>센서 ${argminTol(rows[2].e) + 1}</b>을 지우라고 하지만, 실제로 기여가 가장 작은 센서는 <b>센서 ${kRef + 1}</b>입니다. 근사식이 순위를 뒤바꿨어요.`
    : `근사식이 이 경우엔 같은 센서(센서 ${kRef + 1})를 고르지만, 값 자체는 ${rows[2].e.every((v, i) => Math.abs(v - eDef[i]) < 1e-6) ? '정확히 일치합니다 (Σ가 대각이고 L이 좌표를 섞지 않는 경우).' : '정의식과 다릅니다. 상관 ρ를 키우거나 분해 방식을 바꿔 보세요.'}`;
  $('toyMats').innerHTML = [
    ['\\Sigma=', texMat(Sg)], ['\\Sigma^{-1}=', texMat(W)], ['L=', texMat(L)],
    ['F=S^T\\Sigma^{-1}S=', texMat(matmul(transpose(S), matmul(W, S)), 2)],
  ].map(([a, b]) => `<div>${tex(a + b)}</div>`).join('');
  drawVectors(S, kRef);
}

function drawVectors(S, kRef) {
  const svg = $('vecplot'); clear(svg);
  const W = 300, c = W / 2;
  const mx = Math.max(1, ...S.flatMap((r) => [Math.abs(r[0]), Math.abs(r[1])])) * 1.15, sc = (c - 18) / mx;
  el('line', { x1: 10, y1: c, x2: W - 10, y2: c, stroke: '#d0d7de' }, svg);
  el('line', { x1: c, y1: 10, x2: c, y2: W - 10, stroke: '#d0d7de' }, svg);
  el('text', { x: W - 22, y: c - 4 }, svg, 'θ₁'); el('text', { x: c + 4, y: 18 }, svg, 'θ₂');
  const defs = el('defs', {}, svg);
  SENSOR_COL.forEach((col, i) => {
    const m = el('marker', { id: 'ah' + i, viewBox: '0 0 10 10', refX: 9, refY: 5, markerWidth: 7, markerHeight: 7, orient: 'auto-start-reverse' }, defs);
    el('path', { d: 'M0,0 L10,5 L0,10 z', fill: col }, m);
  });
  S.forEach((r, i) => {
    const x = c + r[0] * sc, y = c - r[1] * sc;
    el('line', { x1: c, y1: c, x2: x, y2: y, stroke: SENSOR_COL[i], 'stroke-width': i === kRef ? 2 : 3, 'stroke-dasharray': i === kRef ? '5 4' : '', 'marker-end': `url(#ah${i})` }, svg);
    el('text', { x: x + 6, y: y - 6, style: `fill:${SENSOR_COL[i]};font-weight:700;font-size:13px` }, svg, `s${i + 1}` + (i === kRef ? ' (제거)' : ''));
  });
  el('text', { x: 8, y: W - 8 }, svg, '점선 = 정의식 기준 제거 대상');
}

function initToy() {
  buildToyInputs(); syncToyControls();
  document.querySelectorAll('[data-preset]').forEach((b) => b.addEventListener('click', () => {
    const p = PRESETS[b.dataset.preset];
    toy.S = p.S.map((r) => [...r]); toy.v = [...p.v]; toy.r = { ...p.r }; toy.factor = p.factor;
    buildToyInputs(); syncToyControls(); renderToy();
  }));
  ['v1', 'v2', 'v3'].forEach((id, i) => $(id).addEventListener('input', (e) => { toy.v[i] = +e.target.value; syncToyControls(); renderToy(); }));
  ['r12', 'r13', 'r23'].forEach((id) => $(id).addEventListener('input', (e) => { toy.r[id] = +e.target.value; syncToyControls(); renderToy(); }));
  $('factor').addEventListener('change', (e) => { toy.factor = e.target.value; renderToy(); });
  renderToy();
}

// ───────────────────────── ③ 연속보 ─────────────────────────
const B = { beam: null, Phi: null, modes: null, freq: null, Sigma: null, res: {}, N: 20, mode: 0, playing: null };
const TARGET = () => B.Phi[0].length;

function setStatus(t) { $('status').textContent = t; }

function rebuildBeam() {
  const spans = $('spans').value.split(',').map((s) => parseFloat(s)).filter((v) => v > 0);
  if (!spans.length || spans.some((s) => Math.abs(s * 2 - Math.round(s * 2)) > 1e-9)) { setStatus('경간은 0.5의 배수로 입력해 주세요.'); return false; }
  const total = spans.reduce((a, b) => a + b, 0);
  if (total * 2 + 1 > 141) { setStatus('후보 위치가 너무 많아요 (총 길이 70 이하).'); return false; }
  B.beam = buildBeam(spans, 0.5);
  const sel = selectModes(B.beam, 10);
  B.Phi = sel.Phi; B.modes = sel.modes; B.freq = sel.freqRatio;
  const ms = $('modeSel'); ms.innerHTML = '';
  B.freq.forEach((f, i) => { const o = document.createElement('option'); o.value = i; o.textContent = `모드 ${i + 1} (ω/ω₁ = ${f.toFixed(2)})`; ms.appendChild(o); });
  B.mode = 0;
  const n = B.beam.nNode, ns = $('nSens');
  ns.min = TARGET(); ns.max = n; B.N = Math.min(Math.max(B.N, TARGET()), n); ns.value = B.N;
  return true;
}

let runToken = 0;
function compute() {
  const token = ++runToken;
  B.Sigma = buildSigma(B.beam.x, B.Phi, { lambda: +$('lambda').value, model: $('covModel').value });
  B.res = {};
  setStatus('Alt.2–4 계산 중…');
  setTimeout(() => {
    if (token !== runToken) return;
    try {
      for (const k of ['alt2', 'alt3', 'alt4']) B.res[k] = sequentialDeletion(B.Phi, B.Sigma, k, TARGET());
    } catch (e) { setStatus('계산 실패: ' + e.message); return; }
    renderBeamAll();
    if ($('runAlt1').checked) {
      setStatus('Alt.1(정의식) 계산 중… 잠시만요');
      setTimeout(() => {
        if (token !== runToken) return;
        B.res.alt1 = sequentialDeletion(B.Phi, B.Sigma, 'alt1', TARGET());
        const same = B.res.alt1.removed.every((v, i) => v === B.res.alt2.removed[i]);
        setStatus(same ? '완료 ✓ Alt.1과 Alt.2의 제거 순서가 완전히 동일합니다.' : '완료 (Alt.1/Alt.2 순서 일부 상이 — 수치 동점)');
        renderBeamAll();
      }, 40);
    } else setStatus('완료 (Alt.1 생략)');
  }, 30);
}

function renderBeamAll() { drawBeam(); drawEntropy(); drawTimes(); drawEfi(); drawSigma(); drawEntTable(); }

function drawBeam() {
  const svg = $('beamPlot'); clear(svg);
  const W = svgW(svg), pad = 40, L = B.beam.L, X = (x) => pad + (x / L) * (W - 2 * pad);
  const md = B.modes[B.mode].w, amp = Math.max(...md.map(Math.abs)) || 1, y0 = 75, Y = (v) => y0 - (v / amp) * 55;
  el('line', { x1: X(0), y1: y0, x2: X(L), y2: y0, stroke: '#8c959f', 'stroke-width': 2 }, svg);
  el('polyline', { points: B.beam.x.map((x, i) => `${X(x)},${Y(md[i])}`).join(' '), fill: 'none', stroke: '#0969da', 'stroke-width': 2 }, svg);
  B.beam.supports.forEach((n, q) => {
    const x = X(B.beam.x[n]);
    el('path', { d: `M${x},${y0} l-7,12 h14 z`, fill: q === 0 ? '#57606a' : '#fff', stroke: '#57606a' }, svg);
    if (q > 0) el('circle', { cx: x, cy: y0 + 15, r: 2.5, fill: 'none', stroke: '#57606a' }, svg);
  });
  el('text', { x: pad, y: 14 }, svg, `모드 ${B.mode + 1} 형상 (ω/ω₁ = ${B.freq[B.mode].toFixed(2)})`);
  const legend = $('legend'); legend.innerHTML = '';
  KEYS.forEach((k, r) => {
    const yy = 130 + r * 50, m = METHODS[k];
    el('text', { x: pad, y: yy - 12, style: `fill:${m.color};font-weight:700` }, svg, m.name);
    el('line', { x1: X(0), y1: yy, x2: X(L), y2: yy, stroke: '#d0d7de', 'stroke-width': 3 }, svg);
    B.beam.supports.forEach((n) => el('path', { d: `M${X(B.beam.x[n])},${yy} l-4,7 h8 z`, fill: '#afb8c1' }, svg));
    const res = B.res[k];
    if (!res) { el('text', { x: X(L) - 90, y: yy - 12 }, svg, k === 'alt1' ? '(계산 대기/생략)' : ''); return; }
    const act = activeAt(res, B.N);
    act.forEach((i) => el('circle', { cx: X(B.beam.x[i]), cy: yy, r: 5, fill: m.color, stroke: '#fff', 'stroke-width': 1 }, svg));
    el('text', { x: X(L) - 150, y: yy - 12 }, svg, `−ln det F = ${res.entropyByCount[B.N].toFixed(2)}`);
    const li = document.createElement('span'); li.innerHTML = `<i style="background:${m.color}"></i>${m.short}`; legend.appendChild(li);
  });
  for (let x = 0; x <= L + 1e-9; x += 4) el('text', { x: X(x) - 4, y: 325 }, svg, x);
  svg.setAttribute('height', 330);
}

function axes(svg, W, H, pad, xr, yr, X, Y, xlab, ylab) {
  el('line', { x1: pad.l, y1: H - pad.b, x2: W - pad.r, y2: H - pad.b, stroke: '#8c959f' }, svg);
  el('line', { x1: pad.l, y1: pad.t, x2: pad.l, y2: H - pad.b, stroke: '#8c959f' }, svg);
  const ny = 5;
  for (let i = 0; i <= ny; i++) {
    const v = yr[0] + (yr[1] - yr[0]) * i / ny;
    el('line', { x1: pad.l, y1: Y(v), x2: W - pad.r, y2: Y(v), stroke: '#eaeef2' }, svg);
    el('text', { x: pad.l - 6, y: Y(v) + 4, 'text-anchor': 'end' }, svg, Math.abs(yr[1] - yr[0]) < 2 ? v.toFixed(2) : v.toFixed(1));
  }
  if (xlab) el('text', { x: (W + pad.l) / 2, y: H - 4, 'text-anchor': 'middle' }, svg, xlab);
  if (ylab) el('text', { x: 12, y: pad.t - 6 }, svg, ylab);
}

function drawEntropy() {
  const svg = $('entPlot'); clear(svg);
  const W = svgW(svg), H = 300, pad = { l: 50, r: 20, t: 20, b: 34 };
  const n0 = B.beam.nNode, t = TARGET();
  const vals = KEYS.flatMap((k) => (B.res[k] ? Object.values(B.res[k].entropyByCount) : []));
  if (!vals.length) return;
  let lo = Math.min(...vals), hi = Math.max(...vals); const m = (hi - lo) * 0.05 || 1; lo -= m; hi += m;
  const X = (n) => pad.l + ((n - t) / (n0 - t)) * (W - pad.l - pad.r), Y = (v) => pad.t + (1 - (v - lo) / (hi - lo)) * (H - pad.t - pad.b);
  axes(svg, W, H, pad, [t, n0], [lo, hi], X, Y, '센서 수 (오른쪽 81개에서 시작해 왼쪽으로 제거)', '−ln det F');
  for (let n = t; n <= n0; n += 5) el('text', { x: X(n), y: H - pad.b + 14, 'text-anchor': 'middle' }, svg, n);
  el('line', { x1: X(B.N), y1: pad.t, x2: X(B.N), y2: H - pad.b, stroke: '#0969da', 'stroke-dasharray': '4 3' }, svg);
  KEYS.forEach((k) => {
    const r = B.res[k]; if (!r) return;
    const pts = []; for (let n = t; n <= n0; n++) pts.push(`${X(n)},${Y(r.entropyByCount[n])}`);
    el('polyline', { points: pts.join(' '), fill: 'none', stroke: METHODS[k].color, 'stroke-width': k === 'alt1' ? 5 : 2, 'stroke-opacity': k === 'alt1' ? 0.35 : 1, 'stroke-dasharray': k === 'alt2' ? '6 3' : '' }, svg);
  });
  KEYS.forEach((k, i) => { el('rect', { x: W - pad.r - 160, y: pad.t + i * 16, width: 12, height: 4, fill: METHODS[k].color }, svg); el('text', { x: W - pad.r - 142, y: pad.t + i * 16 + 6 }, svg, METHODS[k].short); });
  const hit = el('rect', { x: pad.l, y: pad.t, width: W - pad.l - pad.r, height: H - pad.t - pad.b, fill: 'transparent', style: 'cursor:crosshair' }, svg);
  hit.addEventListener('click', (ev) => {
    const bx = svg.getBoundingClientRect(); const n = Math.round(t + ((ev.clientX - bx.left - pad.l) / (W - pad.l - pad.r)) * (n0 - t));
    setN(Math.max(t, Math.min(n0, n)));
  });
}

function drawTimes() {
  const svg = $('timePlot'); clear(svg);
  const W = 300, rows = KEYS.filter((k) => B.res[k]); if (!rows.length) return;
  const mx = Math.max(...rows.map((k) => B.res[k].time));
  KEYS.forEach((k, i) => {
    const y = 10 + i * 46, r = B.res[k];
    el('text', { x: 0, y: y + 10, style: `fill:${METHODS[k].color};font-weight:700` }, svg, METHODS[k].short);
    if (!r) { el('text', { x: 0, y: y + 30 }, svg, '—'); return; }
    const w = Math.max(2, (r.time / mx) * (W - 80));
    el('rect', { x: 0, y: y + 16, width: w, height: 16, fill: METHODS[k].color, rx: 3 }, svg);
    el('text', { x: w + 6, y: y + 29 }, svg, r.time < 1000 ? `${r.time.toFixed(0)} ms` : `${(r.time / 1000).toFixed(2)} s`);
  });
}

function drawEntTable() {
  const best = Math.min(...KEYS.filter((k) => B.res[k]).map((k) => B.res[k].entropyByCount[B.N]));
  $('entTable').innerHTML = '<table class="restbl"><tr><th>방법</th><th>−ln det F</th><th>최선 대비</th></tr>' + KEYS.map((k) => {
    const r = B.res[k]; if (!r) return `<tr><td>${METHODS[k].short}</td><td>—</td><td></td></tr>`;
    const v = r.entropyByCount[B.N];
    return `<tr><td style="color:${METHODS[k].color};text-align:left">${METHODS[k].short}</td><td>${v.toFixed(3)}</td><td>${v - best < 1e-9 ? '최선' : '+' + (v - best).toFixed(3)}</td></tr>`;
  }).join('') + '</table>';
}

function drawEfi() {
  const svg = $('efiPlot'); clear(svg);
  if (!B.res.alt2) return;
  const act = activeAt(B.res.alt2, B.N);
  if (act.length <= TARGET()) { el('text', { x: 20, y: 40 }, svg, `N = ${TARGET()}: 더 제거하면 F가 특이해지므로 순차 제거가 끝난 상태입니다.`); return; }
  const S = rowsOf(B.Phi, act), Sg = sub(B.Sigma, act);
  const E = {}; KEYS.forEach((k) => { E[k] = METHODS[k].f(S, Sg); });
  const W = svgW(svg), H = 280, pad = { l: 46, r: 16, t: 24, b: 30 }, L = B.beam.L;
  const all = KEYS.flatMap((k) => Array.from(E[k]));
  const hi = Math.max(1e-6, ...all) * 1.08;
  const X = (x) => pad.l + (x / L) * (W - pad.l - pad.r), Y = (v) => pad.t + (1 - v / hi) * (H - pad.t - pad.b);
  axes(svg, W, H, pad, [0, L], [0, hi], X, Y, '보 위치 x', 'EfI');
  for (let x = 0; x <= L + 1e-9; x += 4) el('text', { x: X(x), y: H - pad.b + 14, 'text-anchor': 'middle' }, svg, x);
  KEYS.forEach((k, q) => {
    const col = METHODS[k].color, xs = act.map((i) => B.beam.x[i]);
    el('polyline', { points: xs.map((x, j) => `${X(x)},${Y(E[k][j])}`).join(' '), fill: 'none', stroke: col, 'stroke-width': k === 'alt1' ? 5 : 1.6, 'stroke-opacity': k === 'alt1' ? 0.35 : 0.9, 'stroke-dasharray': k === 'alt2' ? '6 3' : '' }, svg);
    xs.forEach((x, j) => el('circle', { cx: X(x), cy: Y(E[k][j]), r: k === 'alt1' ? 4 : 2.4, fill: col, 'fill-opacity': k === 'alt1' ? 0.35 : 1 }, svg));
    const j = argminTol(E[k]);
    el('text', { x: X(xs[j]), y: pad.t - 6 + q * 0, 'text-anchor': 'middle', style: `fill:${col};font-size:${16 - q * 2}px` }, svg, '▼');
    el('line', { x1: X(xs[j]), y1: pad.t, x2: X(xs[j]), y2: H - pad.b, stroke: col, 'stroke-opacity': 0.4, 'stroke-dasharray': '2 3' }, svg);
  });
  KEYS.forEach((k, i) => { el('rect', { x: W - pad.r - 150, y: pad.t + 4 + i * 16, width: 12, height: 4, fill: METHODS[k].color }, svg); el('text', { x: W - pad.r - 132, y: pad.t + 10 + i * 16 }, svg, METHODS[k].short); });
}

function drawSigma() {
  const cv = $('sigmaCanvas'), n = B.Sigma.length, sc = Math.max(1, Math.floor(243 / n));
  cv.width = n * sc; cv.height = n * sc;
  const ctx = cv.getContext('2d');
  for (let i = 0; i < n; i++) for (let j = 0; j < n; j++) {
    const v = Math.max(0, Math.min(1, Math.abs(B.Sigma[i][j])));
    const r = Math.round(255 - v * 246), g = Math.round(255 - v * 150), b = Math.round(255 - v * 37);
    ctx.fillStyle = `rgb(${r},${g},${b})`; ctx.fillRect(j * sc, i * sc, sc, sc);
  }
}

function setN(n) {
  B.N = n; $('nSens').value = n; $('nSens').nextElementSibling.textContent = n;
  drawBeam(); drawEntropy(); drawEfi(); drawEntTable();
}

function initBeam() {
  const lam = $('lambda'); const showLam = () => (lam.nextElementSibling.textContent = (+lam.value).toFixed(1)); showLam();
  let deb;
  lam.addEventListener('input', () => { showLam(); clearTimeout(deb); deb = setTimeout(compute, 250); });
  $('covModel').addEventListener('change', compute);
  $('runAlt1').addEventListener('change', compute);
  $('nSens').addEventListener('input', (e) => setN(+e.target.value));
  $('modeSel').addEventListener('change', (e) => { B.mode = +e.target.value; drawBeam(); });
  $('rebuild').addEventListener('click', () => { if (rebuildBeam()) compute(); });
  $('play').addEventListener('click', () => {
    if (B.playing) { clearInterval(B.playing); B.playing = null; $('play').textContent = '▶ 제거 과정 재생'; return; }
    let n = B.beam.nNode; setN(n); $('play').textContent = '■ 정지';
    B.playing = setInterval(() => {
      n -= 1; if (n < TARGET()) { clearInterval(B.playing); B.playing = null; $('play').textContent = '▶ 제거 과정 재생'; return; }
      setN(n);
    }, 160);
  });
  window.addEventListener('resize', () => { if (B.res.alt2) renderBeamAll(); });
  rebuildBeam(); $('nSens').nextElementSibling.textContent = B.N;
  compute();
}

window.addEventListener('DOMContentLoaded', () => {
  initToy(); initBeam();
  const go = () => window.renderMathInElement && renderMathInElement(document.body, { delimiters: [{ left: '$$', right: '$$', display: true }, { left: '$', right: '$', display: false }], ignoredTags: ['script', 'noscript', 'style', 'textarea', 'pre', 'code', 'svg'] });
  if (window.renderMathInElement) { go(); renderToy(); } else window.addEventListener('load', () => { go(); renderToy(); });
});
