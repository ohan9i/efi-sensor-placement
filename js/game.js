// 게임: 진동하는 연속보(논문 4장) 위에 센서를 꽂고 EfI 봇과 정보량(det F)으로 대결
import { invSPD, logdetSPD, matmul, transpose, rowsOf, sub } from './linalg.js';
import { efiProposed, buildSigma, sequentialDeletion, activeAt, argminTol } from './efi.js';
import { buildBeam, selectModes, buildShear, firstModes } from './beam.js';

const $ = (id) => document.getElementById(id);
const COL = {
  bg: '#0b0d14', panel: '#11141c', raised: '#1a1e28', t1: '#eef0f5', t2: '#a3a9b8', t3: '#7a8296',
  line: 'rgba(238,240,245,0.08)', line2: 'rgba(238,240,245,0.16)',
  me: '#f5a524', bot: '#6c9bff', apx: '#d07a55', best: '#4cc38a', bad: '#e5484d',
};

export const LEVELS = [
  { id: 1, name: '9층 건물', desc: '논문 그림 1의 예시입니다. 9개 층 중 3개 층에만 센서를 달아 1–3차 모드를 구분하세요', spans: 'b9:3', model: 'identity', lambda: 3, K: 3 },
  { id: 2, name: '20층 건물', desc: '20층에서 4개 모드를 센서 6개로. 이웃 층끼리 오차가 닮았습니다(Eq. 28, λ = 3층)', spans: 'b20:4', model: 'eq28', lambda: 3, K: 6 },
  { id: 3, name: '30층 건물', desc: '30층에서 5개 모드를 센서 8개로. 오차 상관이 더 넓습니다(Eq. 28, λ = 4층)', spans: 'b30:5', model: 'eq28', lambda: 4, K: 8 },
];
const SPAN_PRESETS = [
  { v: '8,12,12,8', name: '4경간 8–12–12–8 (논문)' },
  { v: '10,14,10', name: '3경간 10–14–10' },
  { v: '16,16', name: '2경간 16–16' },
  { v: '6,9,12,9,6', name: '5경간 6–9–12–9–6' },
  { v: '24', name: '단순보 24' },
  { v: 'b9:3', name: '건물 9층 · 모드 3개' },
  { v: 'b20:4', name: '건물 20층 · 모드 4개' },
  { v: 'b30:5', name: '건물 30층 · 모드 5개' },
];

const G = {
  level: LEVELS[0], spans: '8,12,12,8', model: 'identity', lambda: 4, K: 14,
  beam: null, Phi: null, m: 10, Sigma: null, placed: new Set(), hints: 3, hint: null,
  revealed: false, bots: null, view: 'all', modeIdx: 0, hover: -1, flash: '',
  efi: null, weakest: -1, score: 0, shown: 0, lastScore: null, pulse: new Map(), cur: [],
};
const beamCache = {};

// ─────────── 계산 ───────────
function logdetOf(set) {
  if (set.length < G.m) return -Infinity;
  const idx = [...set].sort((a, b) => a - b);
  const W = invSPD(sub(G.Sigma, idx)), S = rowsOf(G.Phi, idx);
  return logdetSPD(matmul(transpose(S), matmul(W, S)));
}
const effOf = (ld, ref) => (ld === -Infinity ? 0 : Math.exp((ld - ref) / G.m) * 100);
function proxyScore(set) {
  const S = rowsOf(G.Phi, set), A = matmul(S, transpose(S));
  for (let i = 0; i < A.length; i++) A[i][i] += 1e-9;
  return logdetSPD(A);
}
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
        const t = cur.slice(); t[a] = b;
        const v = logdetOf(t);
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
    if (G.spans.startsWith('b')) {
      const [N, m] = G.spans.slice(1).split(':').map(Number);
      const st = buildShear(N);
      beamCache[G.spans] = { beam: st, ...firstModes(st, m) };
    } else {
      const beam = buildBeam(G.spans.split(',').map(Number), 0.5);
      beamCache[G.spans] = { beam, ...selectModes(beam, 10) };
    }
  }
  const c = beamCache[G.spans];
  G.beam = c.beam; G.Phi = c.Phi; G.m = c.Phi[0].length; G.modes = c.modes; G.freq = c.freqRatio;
  G.norm = G.modes.map((md) => Math.max(...md.w.map(Math.abs)) || 1);
  G.phase = G.modes.map((_, k) => (k * 2.399) % (2 * Math.PI));
  G.Sigma = buildSigma(G.beam.x, G.Phi, { lambda: G.lambda, model: G.model });
  G.K = Math.max(G.m, Math.min(G.K, G.beam.nNode));
  G.bld = G.beam.kind === 'building';
  $('gTip').textContent = G.bld ? '층 오른쪽 소켓을 눌러 센서를 꽂으세요 · 논문 그림 1의 층별 센서 예시' : '보 위의 소켓을 눌러 센서를 꽂으세요 · 논문 4장의 4경간 연속보';
  $('fK').min = G.m; $('fK').max = Math.min(30, G.beam.nNode);
  G.placed = new Set(); G.hints = 3; G.hint = null; G.revealed = false; G.flash = '';
  G.modeIdx = Math.min(G.modeIdx, G.m - 1); G.bots = null; G.lastScore = null; G.pulse.clear();
  $('gResult').hidden = true;
  setView(G.view);
  setTimeout(() => {
    const r2 = sequentialDeletion(G.Phi, G.Sigma, 'alt2', G.K), r4 = sequentialDeletion(G.Phi, G.Sigma, 'alt4', G.K);
    const efi = activeAt(r2, G.K), apx = activeAt(r4, G.K);
    const ls1 = localSearch(efi), ls2 = localSearch(apx);
    G.bots = { efi: { set: efi, ld: logdetOf(efi) }, apx: { set: apx, ld: logdetOf(apx) }, best: ls1.ld >= ls2.ld ? ls1 : ls2 };
    update();
  }, 30);
  update();
}

function toggle(i) {
  if (G.placed.has(i)) G.placed.delete(i);
  else if (G.placed.size >= G.K) { G.flash = `센서는 ${G.K}개가 전부입니다. 꽂힌 센서를 눌러 먼저 빼세요`; update(); return; }
  else { G.placed.add(i); G.pulse.set(i, performance.now()); }
  G.flash = ''; G.hint = null; update();
}

function useHint() {
  if (G.hints <= 0 || !G.bots) return;
  const cur = [...G.placed], n = G.beam.nNode;
  let hint;
  if (cur.length < G.K) {
    let bestV = -Infinity, bestB = -1;
    for (let b = 0; b < n; b++) {
      if (G.placed.has(b)) continue;
      const t = cur.concat(b), v = t.length >= G.m ? logdetOf(t) : proxyScore(t);
      if (v > bestV) { bestV = v; bestB = b; }
    }
    hint = { add: bestB };
  } else {
    let bestV = logdetOf(cur), swap = null;
    for (let a = 0; a < cur.length; a++) for (let b = 0; b < n; b++) {
      if (G.placed.has(b)) continue;
      const t = cur.slice(); t[a] = b;
      const v = logdetOf(t);
      if (v > bestV + 1e-9) { bestV = v; swap = { remove: cur[a], add: b }; }
    }
    hint = swap || { none: true };
  }
  G.hints -= 1; G.hint = hint; update();
}

const posLabel = (i) => (G.bld ? `${i + 1}층` : `x = ${G.beam.x[i].toFixed(1)}`);
const starsFor = (eff, botEff) => (eff <= 0 ? 0 : eff >= botEff - 1e-6 ? 3 : eff >= 90 ? 2 : eff >= 70 ? 1 : 0);

// ─────────── 상태 → HUD ───────────
function update() {
  if (!G.beam) return;
  const cur = [...G.placed].sort((a, b) => a - b);
  G.efi = null; G.weakest = -1;
  if (cur.length > G.m) { G.efi = efiProposed(rowsOf(G.Phi, cur), sub(G.Sigma, cur)); G.weakest = cur[argminTol(G.efi)]; }
  G.cur = cur;
  const ld = logdetOf(cur), ref = G.bots ? G.bots.best.ld : null;
  const eff = ref == null ? 0 : effOf(ld, ref);
  const botEff = G.bots ? effOf(G.bots.efi.ld, ref) : 100;
  if (G.lastScore != null && G.lastScore > 0 && Math.abs(eff - G.lastScore) > 0.05 && cur.length >= G.m) {
    const d = eff - G.lastScore, el = $('gDelta');
    el.textContent = (d > 0 ? '+' : '') + d.toFixed(1);
    el.className = 'g-delta ' + (d > 0 ? 'up' : 'down');
    el.style.animation = 'none'; void el.offsetWidth; el.style.animation = '';
  }
  G.lastScore = ref == null ? null : eff;
  G.score = eff;

  const circ = 2 * Math.PI * 84;
  $('gArc').style.strokeDasharray = `${(Math.min(100, eff) / 100) * circ} ${circ}`;
  $('gBotTick').style.strokeDasharray = G.bots ? `3 ${circ}` : `0 ${circ}`;
  $('gBotTick').style.strokeDashoffset = G.bots ? `${-(Math.min(100, botEff) / 100) * circ + 1.5}` : '0';
  const st = cur.length === G.K && G.bots ? starsFor(eff, botEff) : 0;
  $('gStars').innerHTML = [0, 1, 2].map((i) => `<svg viewBox="0 0 16 16" class="star ${i < st ? 'on' : ''}"><path d="M8 1.5l1.9 4.1 4.4.5-3.3 3 .9 4.4L8 11.3 4.1 13.5l.9-4.4-3.3-3 4.4-.5z"/></svg>`).join('');
  $('gEntropy').textContent = ld === -Infinity ? '—' : (-ld).toFixed(2);
  $('gWeak').textContent = G.weakest >= 0 ? posLabel(G.weakest) : '—';
  $('gModesN').textContent = G.m;
  $('gBotMark').textContent = G.bots ? botEff.toFixed(1) : '…';

  $('gPips').innerHTML = Array.from({ length: G.K }, (_, i) => `<i class="${i < G.K - cur.length ? 'on' : ''}"></i>`).join('') + `<span>${G.K - cur.length}</span>`;
  $('gHint').disabled = G.hints <= 0 || !G.bots; $('gHintN').textContent = G.hints;
  $('gSubmit').disabled = !G.bots || cur.length !== G.K;
  let msg = G.flash;
  if (!msg && G.hint) msg = G.hint.none ? '한 개만 바꿔서 나아질 곳이 없습니다. 거의 최적이에요' : G.hint.remove != null ? `${posLabel(G.hint.remove)} 센서를 ${posLabel(G.hint.add)}로 옮겨 보세요` : `${posLabel(G.hint.add)}에 꽂으면 정보가 가장 많이 늘어납니다`;
  if (!msg && !G.bots) msg = '상대 봇이 배치를 계산하는 중';
  else if (!msg && cur.length < G.m) msg = `모드 ${G.m}개를 구분하려면 센서가 최소 ${G.m}개 필요합니다`;
  else if (!msg && cur.length === G.K) msg = '센서를 다 꽂았습니다. 대결을 시작하거나 위치를 다듬으세요';
  $('gMsg').textContent = msg || '';
  if (G.revealed) renderDuel();
}

function tickScore() {
  const target = G.score;
  G.shown += (target - G.shown) * 0.18;
  if (Math.abs(target - G.shown) < 0.05) G.shown = target;
  $('gScore').textContent = G.shown.toFixed(1);
}

// ─────────── 장면 (캔버스) ───────────
let cv, ctx, W = 0, H = 0, dpr = 1;
const deckY = () => H * 0.44;
const amp = () => Math.min(56, H * 0.14);
const padX = () => Math.max(32, W * 0.04);
const X = (x) => padX() + (x / G.beam.L) * (W - 2 * padX());

function resize() {
  const r = cv.getBoundingClientRect();
  dpr = Math.min(2, window.devicePixelRatio || 1);
  W = r.width; H = r.height;
  cv.width = Math.round(W * dpr); cv.height = Math.round(H * dpr);
}

function disp(t) {
  const n = G.beam.nNode, out = new Float64Array(n);
  if (G.view === 'still') return out;
  const ks = G.view === 'one' ? [G.modeIdx] : G.modes.map((_, k) => k);
  for (const k of ks) {
    const w = G.modes[k].w, a = (G.view === 'one' ? 1 : 0.62 / Math.sqrt(k + 1)) / G.norm[k];
    const s = Math.sin(2 * Math.PI * 0.35 * Math.sqrt(G.freq[k]) * t + G.phase[k]);
    for (let i = 0; i < n; i++) out[i] += a * w[i] * s;
  }
  return out;
}

function roundRect(x, y, w, h, r) { ctx.beginPath(); ctx.moveTo(x + r, y); ctx.arcTo(x + w, y, x + w, y + h, r); ctx.arcTo(x + w, y + h, x, y + h, r); ctx.arcTo(x, y + h, x, y, r); ctx.arcTo(x, y, x + w, y, r); ctx.closePath(); }

function draw(now) {
  if (!G.beam || !W) return;
  if (G.bld) return drawBuilding(now);
  const t = now / 1000, d = disp(t), xs = G.beam.x, n = xs.length;
  const y0 = deckY(), A = amp(), ground = Math.min(H * 0.72, y0 + A + 48);
  ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
  ctx.clearRect(0, 0, W, H);

  // 바닥, 거리 눈금
  ctx.strokeStyle = COL.line2; ctx.lineWidth = 1;
  ctx.beginPath(); ctx.moveTo(0, ground); ctx.lineTo(W, ground); ctx.stroke();
  ctx.strokeStyle = COL.line;
  for (let x = 0; x <= G.beam.L + 1e-9; x += 1) { ctx.beginPath(); ctx.moveTo(X(x), ground + 4); ctx.lineTo(X(x), ground + (x % 4 === 0 ? 12 : 8)); ctx.stroke(); }
  ctx.fillStyle = COL.t3; ctx.font = '11px "IBM Plex Mono", monospace'; ctx.textAlign = 'center';
  for (let x = 0; x <= G.beam.L + 1e-9; x += 4) ctx.fillText(x, X(x), ground + 26);

  // 지점: 논문 그림 2처럼 핀 1개 + 롤러 4개
  G.beam.supports.forEach((s, q) => {
    const px = X(xs[s]), top = y0 + 6;
    ctx.strokeStyle = COL.t2; ctx.lineWidth = 1.2;
    ctx.beginPath(); ctx.moveTo(px, top); ctx.lineTo(px - 9, top + 14); ctx.lineTo(px + 9, top + 14); ctx.closePath(); ctx.stroke();
    if (q === 0) {
      ctx.beginPath(); ctx.moveTo(px - 13, top + 14); ctx.lineTo(px + 13, top + 14); ctx.stroke();
      ctx.strokeStyle = COL.t3; ctx.lineWidth = 1;
      for (let k = -12; k < 12; k += 4) { ctx.beginPath(); ctx.moveTo(px + k + 4, top + 15); ctx.lineTo(px + k, top + 20); ctx.stroke(); }
    } else {
      [-5, 5].forEach((o) => { ctx.beginPath(); ctx.arc(px + o, top + 17.5, 3, 0, 2 * Math.PI); ctx.stroke(); });
      ctx.beginPath(); ctx.moveTo(px - 13, top + 21); ctx.lineTo(px + 13, top + 21); ctx.stroke();
    }
  });
  ctx.fillStyle = COL.t3; ctx.font = '11px "IBM Plex Sans KR", sans-serif'; ctx.textAlign = 'center';
  ctx.fillText('핀', X(xs[G.beam.supports[0]]), y0 + 40);
  ctx.fillText('롤러', X(xs[G.beam.supports[1]]), y0 + 40);

  if (G.view === 'one') {
    const w = G.modes[G.modeIdx].w, nm = G.norm[G.modeIdx];
    ctx.setLineDash([3, 4]); ctx.strokeStyle = 'rgba(108,155,255,0.4)'; ctx.lineWidth = 1;
    [1, -1].forEach((sg) => { ctx.beginPath(); xs.forEach((x, i) => { const py = y0 - sg * (w[i] / nm) * A; i ? ctx.lineTo(X(x), py) : ctx.moveTo(X(x), py); }); ctx.stroke(); });
    ctx.setLineDash([]);
  }

  // 상판
  const yy = (i) => y0 - d[i] * A;
  ctx.beginPath();
  xs.forEach((x, i) => (i ? ctx.lineTo(X(x), yy(i) - 4) : ctx.moveTo(X(x), yy(i) - 4)));
  for (let i = n - 1; i >= 0; i--) ctx.lineTo(X(xs[i]), yy(i) + 6);
  ctx.closePath(); ctx.fillStyle = '#272d3b'; ctx.fill();
  ctx.beginPath(); xs.forEach((x, i) => (i ? ctx.lineTo(X(x), yy(i) - 4) : ctx.moveTo(X(x), yy(i) - 4)));
  ctx.strokeStyle = 'rgba(238,240,245,0.4)'; ctx.lineWidth = 1; ctx.stroke();

  // 소켓, 센서
  const hintAdd = G.hint && G.hint.add, hintRem = G.hint && G.hint.remove;
  let maxR = 0; for (let i = 0; i < n; i++) maxR = Math.max(maxR, Math.abs(d[i]));
  for (let i = 0; i < n; i++) {
    const px = X(xs[i]), py = yy(i) - 4;
    if (G.placed.has(i)) {
      const p0 = G.pulse.get(i), age = p0 ? (now - p0) / 450 : 1;
      if (age < 1) { ctx.strokeStyle = `rgba(245,165,36,${0.7 * (1 - age)})`; ctx.lineWidth = 1.5; ctx.beginPath(); ctx.arc(px, py - 6, 6 + age * 16, 0, 2 * Math.PI); ctx.stroke(); }
      const sig = Math.abs(d[i]) / Math.max(0.35, maxR);
      ctx.fillStyle = i === G.weakest ? 'rgba(229,72,77,0.45)' : 'rgba(245,165,36,0.45)';
      ctx.fillRect(px - 1, py - 15 - sig * 30, 2, sig * 30);
      ctx.fillStyle = i === G.weakest ? COL.bad : COL.me;
      roundRect(px - 5.5, py - 14, 11, 11, 3); ctx.fill();
      ctx.fillStyle = COL.bg; ctx.beginPath(); ctx.arc(px, py - 8.5, 1.8, 0, 2 * Math.PI); ctx.fill();
    } else {
      const hov = i === G.hover;
      ctx.strokeStyle = hov ? COL.me : 'rgba(238,240,245,0.25)'; ctx.lineWidth = hov ? 1.5 : 1;
      ctx.beginPath(); ctx.arc(px, py - 5, hov ? 5 : 2.5, 0, 2 * Math.PI); ctx.stroke();
    }
    if (i === hintAdd || i === hintRem) {
      ctx.strokeStyle = i === hintRem ? COL.bad : COL.best; ctx.setLineDash([3, 3]); ctx.lineWidth = 1.5;
      ctx.beginPath(); ctx.arc(px, py - 8, 14, 0, 2 * Math.PI); ctx.stroke(); ctx.setLineDash([]);
    }
  }

  // 센서별 EfI 띠
  const top = ground + 52, bh = H - top - 10;
  ctx.font = '12px "IBM Plex Sans KR", sans-serif'; ctx.textAlign = 'left'; ctx.fillStyle = COL.t3;
  ctx.fillText(G.efi ? '센서별 EfI · 빨간색이 가장 쓸모없는 센서' : '센서별 EfI · 센서가 모드 수보다 많아지면 나타납니다', padX() - 16, top - 4);
  ctx.strokeStyle = COL.line; ctx.beginPath(); ctx.moveTo(padX() - 16, H - 10); ctx.lineTo(W - padX() + 16, H - 10); ctx.stroke();
  if (G.efi) G.cur.forEach((i, j) => {
    const h = Math.max(1, G.efi[j] * (bh - 14));
    ctx.fillStyle = i === G.weakest ? COL.bad : 'rgba(245,165,36,0.7)';
    roundRect(X(xs[i]) - 3, H - 10 - h, 6, h, 1.5); ctx.fill();
  });

  // 호버 라벨
  if (G.hover >= 0) {
    const i = G.hover, px = X(xs[i]), py = yy(i) - 40;
    const j = G.cur.indexOf(i), txt = `x ${xs[i].toFixed(1)}${j >= 0 && G.efi ? `  EfI ${G.efi[j].toFixed(2)}` : ''}${G.placed.has(i) ? '  눌러서 빼기' : ''}`;
    ctx.font = '12px "IBM Plex Mono", monospace'; const tw = ctx.measureText(txt).width + 16;
    const bx = Math.min(W - tw - 4, Math.max(4, px - tw / 2));
    ctx.fillStyle = COL.raised; roundRect(bx, py - 24, tw, 24, 6); ctx.fill();
    ctx.strokeStyle = COL.line2; ctx.stroke();
    ctx.fillStyle = COL.t1; ctx.textAlign = 'left'; ctx.fillText(txt, bx + 8, py - 8);
  }
}

// ─────────── 건물 그리기 ───────────
function bgeo() {
  const N = G.beam.nNode, ground = H - 40, fh = Math.min(56, (ground - 84) / N);
  const bw = Math.min(200, W * 0.18), bx0 = Math.max(80, W * 0.2), A = Math.min(70, bw * 0.45);
  const ex0 = Math.max(bx0 + bw + A + 120, W * 0.56), ew = Math.min(360, W - ex0 - 40);
  return { N, ground, fh, bw, bx0, A, ex0, ew };
}
function drawBuilding(now) {
  const t = now / 1000, d = disp(t), g = bgeo(), { N, ground, fh, bw, bx0, A, ex0, ew } = g;
  const fy = (i) => ground - (i + 1) * fh;
  let mx = 0; for (let i = 0; i < N; i++) mx = Math.max(mx, Math.abs(d[i]));
  const sc = mx > 1 ? 1 / mx : 1, dx = (i) => (i < 0 ? 0 : d[i] * A * sc);
  ctx.setTransform(dpr, 0, 0, dpr, 0, 0); ctx.clearRect(0, 0, W, H);
  // 지반
  ctx.strokeStyle = COL.line2; ctx.lineWidth = 1; ctx.beginPath(); ctx.moveTo(bx0 - 60, ground); ctx.lineTo(bx0 + bw + 60, ground); ctx.stroke();
  ctx.strokeStyle = COL.t3;
  for (let x = bx0 - 56; x < bx0 + bw + 56; x += 8) { ctx.beginPath(); ctx.moveTo(x + 6, ground + 1); ctx.lineTo(x, ground + 7); ctx.stroke(); }
  // 모드 하나: 형상 윤곽
  if (G.view === 'one') {
    const w = G.modes[G.modeIdx].w, nm = G.norm[G.modeIdx];
    ctx.setLineDash([3, 4]); ctx.strokeStyle = 'rgba(108,155,255,0.45)';
    [1, -1].forEach((sgn) => { ctx.beginPath(); ctx.moveTo(bx0 + bw, ground); for (let i = 0; i < N; i++) ctx.lineTo(bx0 + bw + sgn * (w[i] / nm) * A, fy(i)); ctx.stroke(); });
    ctx.setLineDash([]);
  }
  // 기둥과 슬래브
  for (let i = 0; i < N; i++) {
    const y1 = i === 0 ? ground : fy(i - 1), y2 = fy(i), a = dx(i - 1), b = dx(i);
    ctx.strokeStyle = 'rgba(238,240,245,0.35)'; ctx.lineWidth = 1.5;
    [0, bw / 2, bw].forEach((o) => { ctx.beginPath(); ctx.moveTo(bx0 + o + a, y1); ctx.lineTo(bx0 + o + b, y2); ctx.stroke(); });
    ctx.fillStyle = '#272d3b'; roundRect(bx0 - 6 + b, y2 - 4, bw + 12, 7, 2); ctx.fill();
    ctx.fillStyle = COL.t3; ctx.font = '11px "IBM Plex Mono", monospace'; ctx.textAlign = 'right';
    if (N <= 12 || (i + 1) % 5 === 0 || i === 0) ctx.fillText(`${i + 1}F`, bx0 - 44, y2 + 4);
  }
  // 소켓, 센서
  const hintAdd = G.hint && G.hint.add, hintRem = G.hint && G.hint.remove;
  for (let i = 0; i < N; i++) {
    const px = bx0 + bw + 18 + dx(i), py = fy(i);
    if (G.placed.has(i)) {
      const p0 = G.pulse.get(i), age = p0 ? (now - p0) / 450 : 1;
      if (age < 1) { ctx.strokeStyle = `rgba(245,165,36,${0.7 * (1 - age)})`; ctx.lineWidth = 1.5; ctx.beginPath(); ctx.arc(px, py, 6 + age * 16, 0, 2 * Math.PI); ctx.stroke(); }
      const sig = Math.abs(d[i]) / Math.max(0.35, mx);
      ctx.fillStyle = i === G.weakest ? 'rgba(229,72,77,0.45)' : 'rgba(245,165,36,0.45)';
      ctx.fillRect(px + 9, py - 1, sig * 30, 2);
      ctx.fillStyle = i === G.weakest ? COL.bad : COL.me; roundRect(px - 5.5, py - 5.5, 11, 11, 3); ctx.fill();
      ctx.fillStyle = COL.bg; ctx.beginPath(); ctx.arc(px, py, 1.8, 0, 2 * Math.PI); ctx.fill();
    } else {
      const hov = i === G.hover;
      ctx.strokeStyle = hov ? COL.me : 'rgba(238,240,245,0.25)'; ctx.lineWidth = hov ? 1.5 : 1;
      ctx.beginPath(); ctx.arc(px, py, hov ? 5 : 3, 0, 2 * Math.PI); ctx.stroke();
    }
    if (i === hintAdd || i === hintRem) {
      ctx.strokeStyle = i === hintRem ? COL.bad : COL.best; ctx.setLineDash([3, 3]); ctx.lineWidth = 1.5;
      ctx.beginPath(); ctx.arc(px, py, 13, 0, 2 * Math.PI); ctx.stroke(); ctx.setLineDash([]);
    }
  }
  // 센서별 EfI (가로 막대)
  ctx.font = '12px "IBM Plex Sans KR", sans-serif'; ctx.textAlign = 'left'; ctx.fillStyle = COL.t3;
  ctx.fillText(G.efi ? '센서별 EfI · 빨간색이 가장 쓸모없는 센서' : '센서별 EfI · 센서가 모드 수보다 많아지면 나타납니다', ex0, ground + 4);
  for (let i = 0; i < N; i++) { ctx.fillStyle = 'rgba(238,240,245,0.05)'; ctx.fillRect(ex0, fy(i) - 1, ew, 2); }
  if (G.efi) G.cur.forEach((i, j) => {
    ctx.fillStyle = i === G.weakest ? COL.bad : 'rgba(245,165,36,0.75)';
    roundRect(ex0, fy(i) - Math.min(4, fh * 0.2), Math.max(2, G.efi[j] * ew), Math.min(8, fh * 0.4), 2); ctx.fill();
    ctx.fillStyle = COL.t2; ctx.font = '11px "IBM Plex Mono", monospace'; ctx.fillText(G.efi[j].toFixed(2), ex0 + ew + 8, fy(i) + 4);
  });
  // 호버
  if (G.hover >= 0) {
    const i = G.hover, j = G.cur.indexOf(i);
    const txt = `${i + 1}층${j >= 0 && G.efi ? `  EfI ${G.efi[j].toFixed(2)}` : ''}${G.placed.has(i) ? '  눌러서 빼기' : ''}`;
    ctx.font = '12px "IBM Plex Mono", monospace'; const tw = ctx.measureText(txt).width + 16;
    const bx = bx0 + bw + A + 40, by = fy(i) - 12;
    ctx.fillStyle = COL.raised; roundRect(bx, by, tw, 24, 6); ctx.fill(); ctx.strokeStyle = COL.line2; ctx.stroke();
    ctx.fillStyle = COL.t1; ctx.textAlign = 'left'; ctx.fillText(txt, bx + 8, by + 16);
  }
}

function slotAt(ev) {
  const r = cv.getBoundingClientRect(), mx = ev.clientX - r.left, my = ev.clientY - r.top;
  if (G.bld) {
    const g = bgeo(), i = Math.round((g.ground - my) / g.fh) - 1;
    if (i < 0 || i >= G.beam.nNode || mx < g.bx0 - 60 || mx > g.ex0 + g.ew) return -1;
    return i;
  }
  if (my < deckY() - amp() - 70 || my > deckY() + amp() + 30) return -1;
  const x = ((mx - padX()) / (W - 2 * padX())) * G.beam.L;
  const i = Math.round(x / G.beam.le);
  if (i < 0 || i >= G.beam.nNode) return -1;
  return Math.abs(X(G.beam.x[i]) - mx) <= Math.max(6, (W - 2 * padX()) / G.beam.nNode) ? i : -1;
}

let visible = true;
function loop(now) {
  requestAnimationFrame(loop);
  if (!visible || $('gameView').hidden) return;
  tickScore(); draw(now);
}

// ─────────── 대결 결과 ───────────
function renderDuel() {
  const box = $('gResult'); box.hidden = false;
  const ref = G.bots.best.ld, me = G.score, bot = effOf(G.bots.efi.ld, ref), apx = effOf(G.bots.apx.ld, ref);
  const entries = [
    { k: 'me', name: '나', sub: '직접 배치', v: me, set: G.cur, c: COL.me },
    { k: 'bot', name: 'EfI 봇', sub: '논문 제안식 Eq.(26)으로 순차 제거', v: bot, set: G.bots.efi.set, c: COL.bot },
    { k: 'apx', name: '근사 봇', sub: '기존 촐레스키 근사 Eq.(18)', v: apx, set: G.bots.apx.set, c: COL.apx },
    { k: 'best', name: '최고 기록', sub: '교환 탐색으로 찾은 배치', v: 100, set: G.bots.best.set, c: COL.best },
  ].sort((a, b) => b.v - a.v);
  $('gVerdict').textContent = me >= bot - 1e-6
    ? (me > bot + 0.05 ? 'EfI 봇을 이겼습니다' : 'EfI 봇과 비겼습니다')
    : me >= apx ? `근사 봇은 이겼지만 EfI 봇에 ${(bot - me).toFixed(1)}점 모자랍니다` : '근사 봇에게도 졌습니다. 빨간 센서부터 옮겨 보세요';
  $('gLesson').textContent = bot > apx + 0.5
    ? `같은 순차 제거라도 정확한 EfI를 쓴 봇이 근사식 봇보다 ${(bot - apx).toFixed(1)}점 높습니다. 논문이 말하는 차이가 바로 이것입니다.`
    : '이 조건에서는 근사식도 거의 같은 답을 냅니다. 오차 상관이 약하면 근사의 좌표 왜곡이 작기 때문입니다.';
  const L = G.beam.L;
  $('gDuel').innerHTML = entries.map((e, r) => `
    <div class="duel-row ${e.k}" style="--d:${r * 80}ms">
      <span class="rank">${r + 1}</span>
      <div class="who"><b>${e.name}</b><span>${e.sub}</span></div>
      <svg class="mini" viewBox="0 0 400 28" preserveAspectRatio="none" aria-hidden="true">
        <line x1="0" y1="19" x2="400" y2="19" stroke="rgba(238,240,245,0.2)"/>
        ${G.bld ? G.beam.x.map((_, i) => `<line x1="${(i / (G.beam.nNode - 1)) * 392 + 4}" y1="15" x2="${(i / (G.beam.nNode - 1)) * 392 + 4}" y2="23" stroke="rgba(238,240,245,0.25)"/>`).join('') : G.beam.supports.map((s) => `<path d="M${(G.beam.x[s] / L) * 392 + 4},19 l-4,7 h8z" fill="#7a8296"/>`).join('')}
        ${e.set.map((i) => `<rect x="${(G.bld ? i / (G.beam.nNode - 1) : G.beam.x[i] / L) * 392 + 1}" y="10" width="6" height="7" rx="1.5" fill="${e.c}"/>`).join('')}
      </svg>
      <div class="bar"><i style="width:${Math.max(0, Math.min(100, e.v))}%;background:${e.c}"></i></div>
      <span class="val">${e.v.toFixed(1)}</span>
    </div>`).join('');
}

// ─────────── 단계, 저장 ───────────
const STORE = 'efi-game-stars-building';
const loadStars = () => { try { return JSON.parse(localStorage.getItem(STORE)) || {}; } catch { return {}; } };
function saveStars(id, st) {
  const all = loadStars(); if ((all[id] || 0) >= st) return;
  all[id] = st; try { localStorage.setItem(STORE, JSON.stringify(all)); } catch {}
  paintLevels();
}
function paintLevels() {
  const all = loadStars();
  $('gLevels').querySelectorAll('button').forEach((b) => {
    const n = all[b.dataset.id] || 0;
    b.querySelector('.lvs').innerHTML = LEVELS.find((l) => l.id === +b.dataset.id).free ? '' : [0, 1, 2].map((i) => `<i class="${i < n ? 'on' : ''}"></i>`).join('');
    b.setAttribute('aria-pressed', +b.dataset.id === G.level.id);
  });
}
function setLevel(lv) {
  G.level = lv; Object.assign(G, { spans: lv.spans, model: lv.model, lambda: lv.lambda, K: lv.K });
  $('gStep').textContent = `스테이지 ${lv.id} / ${LEVELS.length}`;
  $('gTitle').textContent = lv.name; $('gDesc').textContent = lv.desc;
  $('gFree').hidden = !lv.free; syncFree(); paintLevels(); setup();
  $('gNext').hidden = lv.id === LEVELS.length;
}
function syncFree() {
  $('fSpans').value = G.spans; $('fModel').value = G.model; $('fLambda').value = G.lambda; $('fLambdaOut').textContent = (+G.lambda).toFixed(1);
  $('fK').value = G.K; $('fKOut').textContent = G.K; $('fLambda').disabled = G.model === 'identity';
}
function setView(v) {
  G.view = v;
  document.querySelectorAll('#gModes button').forEach((x) => x.setAttribute('aria-pressed', x.dataset.v === v));
  $('gModeStep').hidden = v !== 'one';
  if (G.freq) $('gModeLabel').textContent = `모드 ${G.modeIdx + 1}  ω/ω₁ ${G.freq[G.modeIdx].toFixed(2)}`;
}

export function initGame() {
  cv = $('scene'); ctx = cv.getContext('2d');
  $('gLevels').innerHTML = LEVELS.map((lv) => `<button data-id="${lv.id}" aria-pressed="false"><span class="lvn">${lv.free ? '∞' : lv.id}</span><span class="lvt">${lv.name}</span><span class="lvs"></span></button>`).join('');
  $('gLevels').querySelectorAll('button').forEach((b) => b.addEventListener('click', () => setLevel(LEVELS.find((l) => l.id === +b.dataset.id))));
  $('fSpans').innerHTML = SPAN_PRESETS.map((p) => `<option value="${p.v}">${p.name}</option>`).join('');
  $('fSpans').addEventListener('change', (e) => { G.spans = e.target.value; setup(); syncFree(); });
  $('fModel').addEventListener('change', (e) => { G.model = e.target.value; syncFree(); setup(); });
  let deb;
  $('fLambda').addEventListener('input', (e) => { G.lambda = +e.target.value; $('fLambdaOut').textContent = G.lambda.toFixed(1); clearTimeout(deb); deb = setTimeout(setup, 250); });
  $('fK').addEventListener('input', (e) => { G.K = +e.target.value; $('fKOut').textContent = G.K; clearTimeout(deb); deb = setTimeout(setup, 250); });
  $('gHint').addEventListener('click', useHint);
  $('gReset').addEventListener('click', () => { G.placed = new Set(); G.hint = null; G.revealed = false; G.flash = ''; $('gResult').hidden = true; update(); });
  $('gSubmit').addEventListener('click', () => {
    G.revealed = true; update();
    if (!G.level.free) saveStars(G.level.id, starsFor(G.score, effOf(G.bots.efi.ld, G.bots.best.ld)));
    setTimeout(() => $('gResult').scrollIntoView({ behavior: 'smooth', block: 'start' }), 60);
  });
  $('gCopy').addEventListener('click', () => { G.placed = new Set(G.bots.efi.set); G.hint = null; G.revealed = false; $('gResult').hidden = true; update(); window.scrollTo({ top: 0, behavior: 'smooth' }); });
  $('gNext').addEventListener('click', () => { const i = LEVELS.indexOf(G.level); setLevel(LEVELS[Math.min(LEVELS.length - 1, i + 1)]); window.scrollTo({ top: 0, behavior: 'smooth' }); });
  document.querySelectorAll('#gModes button').forEach((b) => b.addEventListener('click', () => setView(b.dataset.v)));
  $('gModePrev').addEventListener('click', () => { G.modeIdx = (G.modeIdx + G.m - 1) % G.m; setView('one'); });
  $('gModeNext').addEventListener('click', () => { G.modeIdx = (G.modeIdx + 1) % G.m; setView('one'); });

  cv.addEventListener('mousemove', (ev) => { const i = slotAt(ev); G.hover = i; cv.style.cursor = i >= 0 ? 'pointer' : 'default'; });
  cv.addEventListener('mouseleave', () => { G.hover = -1; });
  cv.addEventListener('click', (ev) => { const i = slotAt(ev); if (i >= 0) toggle(i); });
  new ResizeObserver(resize).observe(cv);
  new IntersectionObserver(([e]) => (visible = e.isIntersecting)).observe(cv);
  resize();
  setLevel(LEVELS[0]);
  requestAnimationFrame(loop);
}
