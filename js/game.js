// 게임: 직접 센서를 배치하고 EfI 봇과 정보량(det F)으로 겨루기
import { invSPD, logdetSPD, matmul, transpose, rowsOf, sub } from './linalg.js';
import { efiProposed, buildSigma, sequentialDeletion, activeAt, argminTol } from './efi.js';
import { buildBeam, selectModes } from './beam.js';

const $ = (id) => document.getElementById(id);
const NS = 'http://www.w3.org/2000/svg';
const C = { t1: '#1f1b16', t2: '#5e554a', t3: '#7b7266', line: 'rgba(31,27,22,0.1)', strong: 'rgba(31,27,22,0.28)', ink: '#2747a8', pencil: '#b83a2b', ochre: '#a0782c', grey: '#9a9286', raised: '#efe9df', card: '#fdfbf7', good: '#2f7a4f' };
function el(tag, attrs = {}, parent, text) {
  const e = document.createElementNS(NS, tag);
  for (const [k, v] of Object.entries(attrs)) if (v != null) e.setAttribute(k, v);
  if (text != null) e.textContent = text;
  if (parent) parent.appendChild(e);
  return e;
}

export const LEVELS = [
  { id: 1, name: '워밍업', desc: '오차가 서로 독립(Σ = I)인 쉬운 판입니다. 센서 14개로 모드 10개를 구분해 보세요', spans: '8,12,12,8', model: 'identity', lambda: 4, K: 14 },
  { id: 2, name: '닮은 오차', desc: '가까운 센서끼리 오차가 상관됩니다(Eq. 29, λ = 4). 붙여 놓으면 정보가 겹쳐요', spans: '8,12,12,8', model: 'eq29', lambda: 4, K: 12 },
  { id: 3, name: '빠듯한 예산', desc: '모드 10개를 센서 10개로. 하나라도 어긋나면 F가 거의 특이해집니다', spans: '8,12,12,8', model: 'eq29', lambda: 4, K: 10 },
  { id: 4, name: '긴 상관', desc: '상관 길이를 λ = 8로 늘렸습니다. 센서를 더 넓게 퍼뜨려야 해요', spans: '8,12,12,8', model: 'eq28', lambda: 8, K: 12 },
  { id: 5, name: '자유 실험', desc: '구조, 공분산, 예산을 직접 바꿔 가며 실험하세요', spans: '8,12,12,8', model: 'eq29', lambda: 4, K: 12, free: true },
];
const SPAN_PRESETS = [
  { v: '8,12,12,8', name: '4경간 8–12–12–8 (논문)' },
  { v: '10,14,10', name: '3경간 10–14–10' },
  { v: '16,16', name: '2경간 16–16' },
  { v: '6,9,12,9,6', name: '5경간 6–9–12–9–6' },
  { v: '24', name: '단순보 24' },
];

const G = {
  level: LEVELS[0], spans: '8,12,12,8', model: 'identity', lambda: 4, K: 14,
  beam: null, Phi: null, m: 10, Sigma: null, placed: new Set(), hints: 3, hint: null,
  revealed: false, bots: null, showModes: 'all', modeIdx: 0, hover: null, flash: '',
};
const beamCache = {};

function logdetOf(set) {
  if (set.length < G.m) return -Infinity;
  const idx = [...set].sort((a, b) => a - b);
  const W = invSPD(sub(G.Sigma, idx)), S = rowsOf(G.Phi, idx);
  return logdetSPD(matmul(transpose(S), matmul(W, S)));
}
const effOf = (ld, ref) => (ld === -Infinity ? 0 : Math.exp((ld - ref) / G.m) * 100);

// 교환 탐색으로 "알려진 최고 기록" 찾기
function localSearch(start) {
  let cur = [...start], best = logdetOf(cur);
  const n = G.beam.nNode;
  for (let pass = 0; pass < 25; pass++) {
    let improved = false;
    for (let a = 0; a < cur.length; a++) {
      const inSet = new Set(cur);
      let bestB = -1, bestV = best;
      for (let b = 0; b < n; b++) {
        if (inSet.has(b)) continue;
        const trial = cur.slice(); trial[a] = b;
        const v = logdetOf(trial);
        if (v > bestV + 1e-10) { bestV = v; bestB = b; }
      }
      if (bestB >= 0) { cur[a] = bestB; best = bestV; improved = true; }
    }
    if (!improved) break;
  }
  return { set: cur.sort((a, b) => a - b), ld: best };
}

function setup() {
  if (!beamCache[G.spans]) {
    const beam = buildBeam(G.spans.split(',').map(Number), 0.5);
    beamCache[G.spans] = { beam, ...selectModes(beam, 10) };
  }
  const c = beamCache[G.spans];
  G.beam = c.beam; G.Phi = c.Phi; G.m = c.Phi[0].length; G.modes = c.modes; G.freq = c.freqRatio;
  G.Sigma = buildSigma(G.beam.x, G.Phi, { lambda: G.lambda, model: G.model });
  G.K = Math.max(G.m, Math.min(G.K, G.beam.nNode));
  G.placed = new Set(); G.hints = 3; G.hint = null; G.revealed = false; G.flash = ''; G.modeIdx = Math.min(G.modeIdx, G.m - 1);
  G.bots = null;
  $('gModeIdx').max = G.m;
  $('gStatus').textContent = '상대 계산 중';
  setTimeout(() => {
    const r2 = sequentialDeletion(G.Phi, G.Sigma, 'alt2', G.K), r4 = sequentialDeletion(G.Phi, G.Sigma, 'alt4', G.K);
    const efi = activeAt(r2, G.K), apx = activeAt(r4, G.K);
    const ls1 = localSearch(efi), ls2 = localSearch(apx);
    const best = ls1.ld >= ls2.ld ? ls1 : ls2;
    G.bots = { efi: { set: efi, ld: logdetOf(efi) }, apx: { set: apx, ld: logdetOf(apx) }, best };
    $('gStatus').textContent = '';
    render();
  }, 30);
  render();
}

function stars(eff, botEff) {
  if (eff <= 0) return 0;
  if (eff >= botEff - 1e-6) return 3;
  if (eff >= 90) return 2;
  if (eff >= 70) return 1;
  return 0;
}

function toggle(i) {
  if (G.placed.has(i)) G.placed.delete(i);
  else if (G.placed.size >= G.K) { G.flash = `예산은 ${G.K}개입니다. 먼저 센서 하나를 빼세요`; render(); return; }
  else G.placed.add(i);
  G.flash = ''; G.hint = null; render();
}

function useHint() {
  if (G.hints <= 0) return;
  const cur = [...G.placed], n = G.beam.nNode;
  let hint = null;
  if (cur.length < G.K) {
    // 더했을 때 det F를 가장 키우는 자리 (센서가 모드 수보다 적으면 랭크를 가장 키우는 자리)
    let bestV = -Infinity, bestB = -1;
    for (let b = 0; b < n; b++) {
      if (G.placed.has(b)) continue;
      const trial = cur.concat(b);
      const v = trial.length >= G.m ? logdetOf(trial) : proxyScore(trial);
      if (v > bestV) { bestV = v; bestB = b; }
    }
    hint = { add: bestB };
  } else {
    const base = logdetOf(cur);
    let bestV = base, swap = null;
    for (let a = 0; a < cur.length; a++) for (let b = 0; b < n; b++) {
      if (G.placed.has(b)) continue;
      const t = cur.slice(); t[a] = b;
      const v = logdetOf(t);
      if (v > bestV + 1e-9) { bestV = v; swap = { remove: cur[a], add: b }; }
    }
    hint = swap || { none: true };
  }
  G.hints -= 1; G.hint = hint; render();
}
// 센서가 모드 수보다 적을 때: Φ_A Φ_Aᵀ 의 행렬식 (작은 정칙화) 으로 근사
function proxyScore(set) {
  const S = rowsOf(G.Phi, set), A = matmul(S, transpose(S));
  for (let i = 0; i < A.length; i++) A[i][i] += 1e-9;
  return logdetSPD(A);
}

// ─────────── 그리기 ───────────
function drawBoard() {
  const svg = $('board');
  while (svg.firstChild) svg.removeChild(svg.firstChild);
  const W = svg.getBoundingClientRect().width || 1000;
  const rows = G.revealed && G.bots ? 3 : 0;
  const yB = 230, H = yB + 64 + rows * 40;
  svg.setAttribute('viewBox', `0 0 ${W} ${H}`); svg.setAttribute('height', H);
  const L = G.beam.L, padL = 24, padR = 24, X = (x) => padL + (x / L) * (W - padL - padR);
  const xs = G.beam.x;

  // 모드 형상 (배경)
  const modeY = 110, amp = 70;
  const drawMode = (w, op, sw) => {
    const mx = Math.max(...w.map(Math.abs)) || 1;
    el('polyline', { points: xs.map((x, i) => `${X(x)},${modeY - (w[i] / mx) * amp}`).join(' '), fill: 'none', stroke: C.ink, 'stroke-opacity': op, 'stroke-width': sw }, svg);
  };
  el('line', { x1: X(0), y1: modeY, x2: X(L), y2: modeY, stroke: C.line }, svg);
  if (G.showModes === 'all') G.modes.forEach((md) => drawMode(md.w, 0.16, 1.2));
  else if (G.showModes === 'one') drawMode(G.modes[G.modeIdx].w, 0.8, 2);

  // 센서별 EfI 막대 (내 배치)
  const cur = [...G.placed].sort((a, b) => a - b);
  let efi = null, weakest = -1;
  if (cur.length > G.m) {
    efi = efiProposed(rowsOf(G.Phi, cur), sub(G.Sigma, cur));
    weakest = cur[argminTol(efi)];
  }
  if (efi) cur.forEach((i, j) => {
    const h = efi[j] * 52;
    el('rect', { x: X(xs[i]) - 2.5, y: yB - 14 - h, width: 5, height: Math.max(1, h), fill: i === weakest ? C.pencil : C.ink, opacity: i === weakest ? 0.9 : 0.35 }, svg);
  });

  // 보와 지점
  el('line', { x1: X(0), y1: yB, x2: X(L), y2: yB, stroke: C.t1, 'stroke-width': 1.5 }, svg);
  G.beam.supports.forEach((n, q) => {
    const x = X(xs[n]);
    el('path', { d: `M${x},${yB} l-7,12 h14 z`, fill: 'none', stroke: C.t2 }, svg);
    el('line', { x1: x - 9, y1: yB + 15, x2: x + 9, y2: yB + 15, stroke: C.t3, 'stroke-width': 0.8 }, svg);
  });
  // 치수
  const yd = yB + 40, sx = G.beam.supports.map((n) => xs[n]);
  el('line', { x1: X(0), y1: yd, x2: X(L), y2: yd, stroke: C.t3, 'stroke-width': 0.8 }, svg);
  sx.forEach((x) => { el('line', { x1: X(x), y1: yd - 5, x2: X(x), y2: yd + 5, stroke: C.t3, 'stroke-width': 0.8 }, svg); el('line', { x1: X(x) - 4, y1: yd + 4, x2: X(x) + 4, y2: yd - 4, stroke: C.t3, 'stroke-width': 0.8 }, svg); });
  for (let i = 0; i < sx.length - 1; i++) el('text', { x: (X(sx[i]) + X(sx[i + 1])) / 2, y: yd - 6, 'text-anchor': 'middle' }, svg, (sx[i + 1] - sx[i]).toFixed(1));

  // 슬롯
  const slotW = Math.max(6, Math.min(10, (W - padL - padR) / xs.length - 2));
  xs.forEach((x, i) => {
    const on = G.placed.has(i), cx = X(x);
    const isHint = G.hint && (G.hint.add === i || G.hint.remove === i);
    if (isHint) el('circle', { cx, cy: yB - 6, r: 11, fill: 'none', stroke: G.hint.remove === i ? C.pencil : C.ink, 'stroke-width': 1.5, 'stroke-dasharray': '3 2' }, svg);
    el('rect', { x: cx - slotW / 2, y: yB - slotW - 2, width: slotW, height: slotW, rx: 1.5, fill: on ? C.t1 : C.card, stroke: on ? C.t1 : G.hover === i ? C.ink : C.strong, 'stroke-width': 1 }, svg);
    const hit = el('rect', { x: cx - 6, y: yB - 70, width: 12, height: 84, fill: 'transparent', style: 'cursor:pointer' }, svg);
    hit.addEventListener('click', () => toggle(i));
    hit.addEventListener('mouseenter', () => { G.hover = i; $('gHover').textContent = `x = ${x.toFixed(1)}${on && efi ? ` · EfI ${efi[cur.indexOf(i)].toFixed(3)}` : ''}`; });
  });
  if (weakest >= 0) el('text', { x: X(xs[weakest]), y: yB - 74, 'text-anchor': xs[weakest] > L * 0.85 ? 'end' : xs[weakest] < L * 0.15 ? 'start' : 'middle', class: 'lab', style: `fill:${C.pencil}` }, svg, '가장 약한 센서');

  // 공개: 상대 배치
  if (rows) {
    const bots = [['EfI 봇 (정확식)', G.bots.efi.set, C.ink], ['근사 봇 (촐레스키)', G.bots.apx.set, C.ochre], ['탐색 최고 기록', G.bots.best.set, C.good]];
    bots.forEach(([name, set, col], r) => {
      const y = yB + 80 + r * 40;
      el('text', { x: X(0), y: y - 10, class: 'lab' }, svg, name);
      el('line', { x1: X(0), y1: y, x2: X(L), y2: y, stroke: C.strong }, svg);
      set.forEach((i) => el('rect', { x: X(xs[i]) - 3.5, y: y - 8, width: 7, height: 7, fill: col, opacity: G.placed.has(i) ? 1 : 0.85 }, svg));
      set.forEach((i) => { if (G.placed.has(i)) el('rect', { x: X(xs[i]) - 5, y: y - 9.5, width: 10, height: 10, fill: 'none', stroke: C.t1, 'stroke-width': 0.8 }, svg); });
    });
  }
}

function render() {
  if (!G.beam) return;
  const cur = [...G.placed];
  const ld = logdetOf(cur);
  const ref = G.bots ? G.bots.best.ld : null;
  const eff = ref != null ? effOf(ld, ref) : null;
  const botEff = G.bots ? effOf(G.bots.efi.ld, ref) : null;
  $('gUsed').textContent = `${cur.length} / ${G.K}`;
  $('gUsedBar').style.width = `${(cur.length / G.K) * 100}%`;
  $('gScore').textContent = eff == null ? '…' : cur.length < G.m ? '0.0' : eff.toFixed(1);
  $('gEntropy').textContent = ld === -Infinity ? `센서 ${G.m - cur.length}개 더 필요` : `−ln det F  ${(-ld).toFixed(2)}`;
  const st = eff != null && cur.length === G.K ? stars(eff, botEff) : 0;
  if (G.revealed && !G.level.free) saveStars(G.level.id, st);
  $('gStars').innerHTML = [0, 1, 2].map((i) => `<svg viewBox="0 0 16 16" class="star ${i < st ? 'on' : ''}"><path d="M8 1.5l1.9 4.1 4.4.5-3.3 3 .9 4.4L8 11.3 4.1 13.5l.9-4.4-3.3-3 4.4-.5z"/></svg>`).join('');
  $('gHint').disabled = G.hints <= 0 || !G.bots;
  $('gHintN').textContent = G.hints;
  $('gSubmit').disabled = !G.bots || cur.length !== G.K;
  $('gSubmit').textContent = G.revealed ? '결과 다시 보기' : '제출하고 대결';
  let msg = G.flash;
  if (!msg && G.hint) msg = G.hint.none ? '지금 배치에서는 한 개만 바꿔서 나아질 곳이 없습니다' : G.hint.remove != null ? `x = ${G.beam.x[G.hint.remove].toFixed(1)}을 빼고 x = ${G.beam.x[G.hint.add].toFixed(1)}에 두면 좋아집니다` : `x = ${G.beam.x[G.hint.add].toFixed(1)}에 두면 정보가 가장 많이 늘어납니다`;
  if (!msg && cur.length < G.m) msg = `모드 ${G.m}개를 구분하려면 센서가 최소 ${G.m}개 필요합니다`;
  $('gMsg').textContent = msg || '';
  $('gCopy').hidden = !G.revealed;
  drawBoard();
  renderResult(eff, ld);
}

function renderResult(eff, ld) {
  const box = $('gResult');
  if (!G.revealed || !G.bots) { box.hidden = true; return; }
  box.hidden = false;
  const ref = G.bots.best.ld;
  const rows = [
    ['나', eff, C.t1],
    ['EfI 봇 · 정확식 Eq.(26)', effOf(G.bots.efi.ld, ref), C.ink],
    ['근사 봇 · 촐레스키 Eq.(18)', effOf(G.bots.apx.ld, ref), C.ochre],
    ['탐색 최고 기록', 100, C.good],
  ].sort((a, b) => b[1] - a[1]);
  const meE = eff, botE = effOf(G.bots.efi.ld, ref);
  const verdict = meE >= botE - 1e-6
    ? (meE > botE + 0.05 ? 'EfI 봇을 이겼습니다. 순차 제거는 탐욕적이라 전역 최적을 보장하지 않거든요' : 'EfI 봇과 같은 수준입니다')
    : meE >= 90 ? `EfI 봇보다 ${(botE - meE).toFixed(1)}점 낮습니다. 빨간 막대(가장 약한 센서)를 옮겨 보세요` : '정보가 겹치는 센서가 많습니다. 모드 형상이 크게 움직이는 곳에 고르게 퍼뜨려 보세요';
  box.innerHTML = `<p class="gr-verdict">${verdict}</p><div class="gr-list">${rows.map(([n, v, c]) => `
    <div class="gr-row${n === '나' ? ' me' : ''}"><span>${n}</span><span class="v">${v.toFixed(1)}</span><span class="t"><i style="width:${Math.max(0, Math.min(100, v))}%;background:${c}"></i></span></div>`).join('')}</div>
    <p class="gr-note">점수는 탐색으로 찾은 최고 배치 대비 정보 효율, 즉 (det F / det F<sub>best</sub>)<sup>1/${G.m}</sup> × 100 입니다. 아래 도면의 테두리 친 칸은 내 배치와 겹치는 위치입니다.</p>`;
}

const STORE = 'efi-game-stars';
const loadStars = () => { try { return JSON.parse(localStorage.getItem(STORE)) || {}; } catch { return {}; } };
function saveStars(id, st) {
  const all = loadStars();
  if ((all[id] || 0) >= st) return;
  all[id] = st; try { localStorage.setItem(STORE, JSON.stringify(all)); } catch {}
  paintLevelStars();
}
function paintLevelStars() {
  const all = loadStars();
  document.querySelectorAll('#gLevels button').forEach((b) => {
    const n = all[b.dataset.id] || 0, slot = b.querySelector('.lvs');
    if (slot) slot.textContent = n ? '★'.repeat(n) : '';
  });
}

function setLevel(lv) {
  G.level = lv;
  Object.assign(G, { spans: lv.spans, model: lv.model, lambda: lv.lambda, K: lv.K });
  document.querySelectorAll('#gLevels button').forEach((b) => b.setAttribute('aria-pressed', +b.dataset.id === lv.id));
  $('gDesc').textContent = lv.desc;
  $('gFree').hidden = !lv.free;
  syncFree();
  setup();
}
function syncFree() {
  $('fSpans').value = G.spans; $('fModel').value = G.model; $('fLambda').value = G.lambda; $('fLambdaOut').textContent = (+G.lambda).toFixed(1);
  $('fK').value = G.K; $('fKOut').textContent = G.K;
  $('fLambda').disabled = G.model === 'identity';
}

export function initGame() {
  $('gLevels').innerHTML = LEVELS.map((lv) => `<button data-id="${lv.id}" aria-pressed="false"><span class="lvn">${lv.free ? '∞' : lv.id}</span>${lv.name}<span class="lvs"></span></button>`).join('');
  paintLevelStars();
  $('gLevels').querySelectorAll('button').forEach((b) => b.addEventListener('click', () => setLevel(LEVELS.find((l) => l.id === +b.dataset.id))));
  $('fSpans').innerHTML = SPAN_PRESETS.map((p) => `<option value="${p.v}">${p.name}</option>`).join('');
  $('fSpans').addEventListener('change', (e) => { G.spans = e.target.value; setup(); syncFree(); });
  $('fModel').addEventListener('change', (e) => { G.model = e.target.value; syncFree(); setup(); });
  let deb;
  $('fLambda').addEventListener('input', (e) => { G.lambda = +e.target.value; $('fLambdaOut').textContent = G.lambda.toFixed(1); clearTimeout(deb); deb = setTimeout(setup, 250); });
  $('fK').addEventListener('input', (e) => { G.K = +e.target.value; $('fKOut').textContent = G.K; clearTimeout(deb); deb = setTimeout(setup, 250); });
  $('gHint').addEventListener('click', useHint);
  $('gReset').addEventListener('click', () => { G.placed = new Set(); G.hint = null; G.revealed = false; G.flash = ''; render(); });
  $('gSubmit').addEventListener('click', () => { G.revealed = true; render(); $('gResult').scrollIntoView({ behavior: 'smooth', block: 'nearest' }); });
  $('gCopy').addEventListener('click', () => { if (!G.bots) return; G.placed = new Set(G.bots.efi.set); G.hint = null; render(); });
  document.querySelectorAll('#gModes button').forEach((b) => b.addEventListener('click', () => {
    G.showModes = b.dataset.v; document.querySelectorAll('#gModes button').forEach((x) => x.setAttribute('aria-pressed', x === b)); render();
  }));
  $('gModeIdx').addEventListener('input', (e) => { G.modeIdx = +e.target.value - 1; G.showModes = 'one'; document.querySelectorAll('#gModes button').forEach((x) => x.setAttribute('aria-pressed', x.dataset.v === 'one')); $('gModeIdxOut').textContent = e.target.value; render(); });
  $('board').addEventListener('mouseleave', () => { G.hover = null; $('gHover').textContent = '칸을 눌러 센서를 놓거나 빼세요'; });
  window.addEventListener('resize', () => { clearTimeout(deb); deb = setTimeout(render, 120); });
  setLevel(LEVELS[0]);
}
