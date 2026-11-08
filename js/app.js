import { fromArray, invSPD, invSqrtSPD, cholesky, matmul, rowsOf, sub } from './linalg.js';
import { efiDefinition, efiProposed, efiFactor, METHODS, buildSigma, sequentialDeletion, activeAt, argminTol } from './efi.js';
import { buildBeam, selectModes } from './beam.js';

const $ = (id) => document.getElementById(id);
const NS = 'http://www.w3.org/2000/svg';
const KEYS = ['alt1', 'alt2', 'alt3', 'alt4'];
const C = { t1: '#1f1b16', t2: '#5e554a', t3: '#7b7266', line: 'rgba(31,27,22,0.1)', strong: 'rgba(31,27,22,0.28)', ink: '#2747a8', inkSoft: '#7d93cf', pencil: '#b83a2b', raised: '#efe9df', card: '#fdfbf7' };

function el(tag, attrs = {}, parent, text) {
  const e = document.createElementNS(NS, tag);
  for (const [k, v] of Object.entries(attrs)) if (v != null) e.setAttribute(k, v);
  if (text != null) e.textContent = text;
  if (parent) parent.appendChild(e);
  return e;
}
const clear = (s) => { while (s.firstChild) s.removeChild(s.firstChild); };
const fx = (v, d = 3) => (Math.abs(v) < 5e-12 ? 0 : v).toFixed(d);
const fmtX = (x) => x.toFixed(1);
const fmtMs = (t) => (t < 1000 ? `${Math.round(t)} ms` : `${(t / 1000).toFixed(2)} s`);
const tex = (s) => (window.katex ? katex.renderToString(s, { throwOnError: false }) : s);

// ─────────────── 상태 ───────────────
const S = {
  beam: null, Phi: null, modes: null, freq: null, Sigma: null, res: {},
  N: 20, mode: 0, ch: 1, model: 'eq29', lambda: 4, sel: null, hover: null, playing: null, busy: false,
};
const target = () => S.Phi[0].length;

// ─────────────── 계산 ───────────────
function setupBeam() {
  S.beam = buildBeam([8, 12, 12, 8], 0.5);
  const sel = selectModes(S.beam, 10);
  S.Phi = sel.Phi; S.modes = sel.modes; S.freq = sel.freqRatio;
  S.sel = S.beam.x.findIndex((x) => Math.abs(x - 14) < 1e-9);
  const ns = $('nSens'); ns.min = target(); ns.max = S.beam.nNode;
}

let token = 0;
function compute() {
  const my = ++token;
  S.Sigma = buildSigma(S.beam.x, S.Phi, { lambda: S.lambda, model: S.model });
  S.res = {}; S.busy = true; meta();
  setTimeout(() => {
    if (my !== token) return;
    for (const k of ['alt2', 'alt3', 'alt4']) S.res[k] = sequentialDeletion(S.Phi, S.Sigma, k, target());
    render();
    setTimeout(() => {
      if (my !== token) return;
      S.res.alt1 = sequentialDeletion(S.Phi, S.Sigma, 'alt1', target());
      S.busy = false; render();
    }, 30);
  }, 20);
}
const activeSet = () => (S.res.alt2 ? activeAt(S.res.alt2, S.N) : S.Phi.map((_, i) => i));

// ─────────────── 도면 공통 ───────────────
function frame(svg, h) {
  clear(svg);
  const W = svg.getBoundingClientRect().width || 900;
  svg.setAttribute('viewBox', `0 0 ${W} ${h}`); svg.setAttribute('height', h);
  const padL = 16, padR = 16, L = S.beam.L;
  const X = (x) => padL + (x / L) * (W - padL - padR);
  return { W, h, X, L };
}

// 보, 지점, 치수선 (도면 기호)
function drawBeamLine(svg, f, y, { dims = false, faint = false, dimOff = 32 } = {}) {
  const { X, L } = f;
  el('line', { x1: X(0), y1: y, x2: X(L), y2: y, stroke: faint ? C.strong : C.t1, 'stroke-width': faint ? 1 : 1.5 }, svg);
  S.beam.supports.forEach((n, q) => {
    const x = X(S.beam.x[n]);
    el('path', { d: `M${x},${y} l-6,10 h12 z`, fill: 'none', stroke: faint ? C.strong : C.t2, 'stroke-width': 1 }, svg);
    if (q === 0) for (let k = -6; k < 6; k += 3) el('line', { x1: x + k, y1: y + 13, x2: x + k + 3, y2: y + 10, stroke: C.t3, 'stroke-width': 0.8 }, svg);
    else el('line', { x1: x - 8, y1: y + 13, x2: x + 8, y2: y + 13, stroke: faint ? C.strong : C.t3, 'stroke-width': 0.8 }, svg);
  });
  if (!dims) return;
  const yd = y + dimOff;
  const sx = S.beam.supports.map((n) => S.beam.x[n]);
  el('line', { x1: X(0), y1: yd, x2: X(L), y2: yd, stroke: C.t3, 'stroke-width': 0.8 }, svg);
  sx.forEach((x) => {
    el('line', { x1: X(x), y1: yd - 6, x2: X(x), y2: yd + 6, stroke: C.t3, 'stroke-width': 0.8 }, svg);
    el('line', { x1: X(x) - 4, y1: yd + 4, x2: X(x) + 4, y2: yd - 4, stroke: C.t3, 'stroke-width': 0.8 }, svg);
  });
  for (let i = 0; i < sx.length - 1; i++) el('text', { x: (X(sx[i]) + X(sx[i + 1])) / 2, y: yd - 6, 'text-anchor': 'middle' }, svg, fx(sx[i + 1] - sx[i], 1));
}

// 센서 기호: 보 위의 작은 사각형
function sensor(svg, x, y, color = C.t1, size = 6) {
  return el('rect', { x: x - size / 2, y: y - size - 1, width: size, height: size, fill: color }, svg);
}
function pencil(svg, x, y, r = 10, label) {
  el('ellipse', { cx: x, cy: y, rx: r + 1.5, ry: r - 1, fill: 'none', stroke: C.pencil, 'stroke-width': 1.5, transform: `rotate(-12 ${x} ${y})` }, svg);
  if (label) el('text', { x: x + r + 6, y: y - r + 2, style: `fill:${C.pencil};font-family:var(--sans);font-size:12px` }, svg, label);
}

// ─────────────── 챕터별 도면 ───────────────
const FIG = {
  1() {
    const md = S.modes[S.mode].w, amp = Math.max(...md.map(Math.abs)) || 1;
    $('figTitle').innerHTML = `<b>모드 ${S.mode + 1}</b> · ω/ω₁ = ${S.freq[S.mode].toFixed(2)} · 후보 위치 ${S.beam.nNode}곳`;
    const svg = $('beam'), f = frame(svg, 236), y = 100, X = f.X;
    for (let i = 0; i < S.beam.nNode; i++) el('line', { x1: X(S.beam.x[i]), y1: y - 3, x2: X(S.beam.x[i]), y2: y + 3, stroke: C.strong, 'stroke-width': 1 }, svg);
    drawBeamLine(svg, f, y, { dims: true, dimOff: 124 });
    const pts = S.beam.x.map((x, i) => `${X(x)},${y - (md[i] / amp) * 80}`).join(' ');
    el('polyline', { points: pts, fill: 'none', stroke: C.ink, 'stroke-width': 2, 'stroke-linejoin': 'round' }, svg);
  },
  2() {
    const act = activeSet(), Sm = rowsOf(S.Phi, act), Sg = sub(S.Sigma, act);
    const e = efiProposed(Sm, Sg), k = argminTol(e);
    S._efi = { act, e, k };
    const svg = $('beam'), f = frame(svg, 260), y = 200, X = f.X, H = 150;
    $('figTitle').innerHTML = `센서 <b>${act.length}개</b>의 EfI · 막대 높이 = 빼면 잃는 det F 비율`;
    [0.25, 0.5, 0.75, 1].forEach((v) => {
      el('line', { x1: X(0), y1: y - v * H, x2: X(f.L), y2: y - v * H, stroke: C.line }, svg);
      el('text', { x: X(f.L), y: y - v * H - 4, 'text-anchor': 'end' }, svg, v.toFixed(2));
    });
    drawBeamLine(svg, f, y, { dims: true });
    const hv = S.hover;
    act.forEach((i, j) => {
      const x = X(S.beam.x[i]), isMin = j === k, isH = hv === i;
      el('rect', { x: x - 2.5, y: y - 2 - e[j] * H, width: 5, height: Math.max(1, e[j] * H), fill: isMin ? C.pencil : isH ? C.ink : C.t2, opacity: isMin || isH ? 1 : 0.55 }, svg);
      sensor(svg, x, y, C.t1, 5);
      const hit = el('rect', { x: x - 7, y: y - H - 10, width: 14, height: H + 14, fill: 'transparent', style: 'cursor:pointer' }, svg);
      hit.addEventListener('mouseenter', () => { if (S.hover === i) return; S.hover = i; FIG[2](); PANEL[2](); });
    });
    pencil(svg, X(S.beam.x[act[k]]), y - 6, 9, '다음 제거');
  },
  3() {
    const svg = $('beam'), f = frame(svg, 240), y = 190, X = f.X, H = 130;
    const row = S.Sigma[S.sel];
    $('figTitle').innerHTML = `x = <b>${fmtX(S.beam.x[S.sel])}</b> 센서와의 오차 상관 ρ · 보를 클릭해 바꾸기`;
    [0.5, 1].forEach((v) => { el('line', { x1: X(0), y1: y - v * H, x2: X(f.L), y2: y - v * H, stroke: C.line }, svg); el('text', { x: X(f.L), y: y - v * H - 4, 'text-anchor': 'end' }, svg, v.toFixed(1)); });
    const pts = S.beam.x.map((x, i) => `${X(x)},${y - Math.max(0, row[i]) * H}`);
    el('polygon', { points: `${X(0)},${y} ${pts.join(' ')} ${X(f.L)},${y}`, fill: C.ink, 'fill-opacity': 0.1 }, svg);
    el('polyline', { points: pts.join(' '), fill: 'none', stroke: C.ink, 'stroke-width': 1.5 }, svg);
    drawBeamLine(svg, f, y, { dims: true });
    sensor(svg, X(S.beam.x[S.sel]), y, C.ink, 8);
    const hit = el('rect', { x: 0, y: 0, width: f.W, height: f.h, fill: 'transparent', style: 'cursor:crosshair' }, svg);
    hit.addEventListener('click', (ev) => {
      const b = svg.getBoundingClientRect(), px = ((ev.clientX - b.left) / b.width) * f.W;
      const x = ((px - X(0)) / (X(f.L) - X(0))) * f.L;
      S.sel = Math.max(0, Math.min(S.beam.nNode - 1, Math.round(x / S.beam.le)));
      FIG[3](); PANEL[3]();
    });
  },
  4() {
    const act = activeSet(), Sm = rowsOf(S.Phi, act), Sg = sub(S.Sigma, act);
    const E = {}; KEYS.forEach((k) => (E[k] = METHODS[k].f(Sm, Sg)));
    const picks = {}; KEYS.forEach((k) => (picks[k] = act[argminTol(E[k])]));
    S._cmp = { act, E, picks };
    const svg = $('beam'), f = frame(svg, 176), y = 112, X = f.X;
    $('figTitle').innerHTML = `센서 <b>${act.length}개</b>에서 각 방법이 다음에 지우려는 센서`;
    drawBeamLine(svg, f, y, { dims: true });
    act.forEach((i) => sensor(svg, X(S.beam.x[i]), y, C.t2, 5));
    const truth = picks.alt2;
    el('line', { x1: X(S.beam.x[truth]), y1: y - 10, x2: X(S.beam.x[truth]), y2: 24, stroke: C.ink, 'stroke-width': 1 }, svg);
    el('text', { x: X(S.beam.x[truth]), y: 16, 'text-anchor': 'middle', class: 'labstrong', style: `fill:${C.ink}` }, svg, '정확한 답 (Alt.1, Alt.2)');
    sensor(svg, X(S.beam.x[truth]), y, C.ink, 8);
    const wrong = {};
    ['alt3', 'alt4'].forEach((k) => { if (picks[k] !== truth) (wrong[picks[k]] ||= []).push(METHODS[k].tag); });
    Object.entries(wrong).forEach(([i, tags], q) => {
      const x = X(S.beam.x[+i]), ly = 56 + q * 24;
      sensor(svg, x, y, C.pencil, 8);
      pencil(svg, x, y - 5, 10);
      el('line', { x1: x, y1: y - 16, x2: x, y2: ly + 4, stroke: C.pencil, 'stroke-width': 1 }, svg);
      el('text', { x, y: ly, 'text-anchor': 'middle', class: 'lab', style: `fill:${C.pencil}` }, svg, `${tags.join(', ')}의 선택`);
    });
  },
  5() {
    const svg = $('beam'), rows = 4, gap = 44, top = 32, f = frame(svg, top + rows * gap + 36), X = f.X;
    $('figTitle').innerHTML = `센서 <b>${S.N}개</b>가 남았을 때 네 방법의 배치`;
    KEYS.forEach((k, r) => {
      const y = top + r * gap, m = METHODS[k], res = S.res[k];
      el('text', { x: X(0), y: y - 16, class: k === 'alt2' || k === 'alt1' ? 'labstrong' : 'lab' }, svg, `${m.tag}  ${m.name}`);
      drawBeamLine(svg, f, y, { faint: true, dims: r === rows - 1 });
      if (!res) { el('text', { x: X(f.L), y: y - 16, 'text-anchor': 'end' }, svg, '계산 중'); return; }
      el('text', { x: X(f.L), y: y - 16, 'text-anchor': 'end' }, svg, `−ln det F  ${res.entropyByCount[S.N].toFixed(2)}`);
      activeAt(res, S.N).forEach((i) => sensor(svg, X(S.beam.x[i]), y, m.color, 7));
    });
  },
};
FIG[6] = FIG[5];

// ─────────────── 챕터별 패널 ───────────────
const PANEL = {
  1() {
    const mx = Math.max(...S.freq);
    $('panel').innerHTML = `<p class="ptitle"><b>목표 모드 ${S.freq.length}개</b> · 1차 진동수 대비 비율</p><div class="list" id="modeList">${S.freq.map((fq, i) => `
      <button aria-pressed="${i === S.mode}" data-i="${i}"><span>모드 ${i + 1}</span><span class="track"><span class="fill" style="display:block;width:${(fq / mx) * 100}%"></span></span><span class="v">${fq.toFixed(2)}</span></button>`).join('')}</div>
      <p class="note">10차 모드까지가 1차의 10배 이하입니다. 진동수가 비슷한 모드는 형상도 짝을 이룹니다 (3–4차, 7–8차).</p>`;
    $('modeList').querySelectorAll('button').forEach((b) => b.addEventListener('click', () => { S.mode = +b.dataset.i; FIG[1](); PANEL[1](); }));
  },
  2() {
    const { act, e, k } = S._efi;
    const j = S.hover != null && act.includes(S.hover) ? act.indexOf(S.hover) : k;
    const i = act[j];
    $('panel').innerHTML = `<p class="ptitle">${j === k ? '<b>가장 작은 EfI</b> · 다음에 제거될 센서' : '<b>선택한 센서</b>'}</p>
      <dl class="readout"><div><dt>위치 x</dt><dd>${fmtX(S.beam.x[i])}</dd></div><div><dt>EfI e<sub>i</sub></dt><dd>${fx(e[j])}</dd></div><div><dt>빼면 det F</dt><dd>×${fx(1 - e[j])}</dd></div></dl>
      <p class="note">센서가 많을 때는 각자의 몫이 작고, 센서가 줄어들수록 남은 센서의 EfI가 1에 가까워집니다. 위쪽 센서 슬라이더로 개수를 바꿔 보세요.</p>`;
  },
  3() {
    $('panel').innerHTML = `<div class="pgrid"><div>
      <p class="ptitle"><b>오차 공분산 Σ</b> · ${S.beam.nNode} × ${S.beam.nNode}</p>
      <p class="note" style="margin-top:0">밝은 대각 띠의 폭이 상관 길이 λ입니다. 띠 안의 센서들은 오차를 공유하므로, 둘을 함께 두어도 정보가 두 배가 되지 않습니다.</p>
      <p class="note">Eq. 29는 모드 형상이 비슷한 위치끼리만 강하게 상관시켜, 지점처럼 응답이 0인 곳에 센서가 남는 부작용을 막습니다. 지점 위치는 십자 모양의 빈 줄로 보입니다.</p>
      <p class="note">지금 설정: ${{ eq29: 'Eq. 29', eq28: 'Eq. 28', identity: 'Σ = I' }[S.model]}${S.model === 'identity' ? '' : `, λ = ${S.lambda.toFixed(1)}`}</p>
      </div><canvas class="heat" id="heat"></canvas></div>`;
    const cv = $('heat'), n = S.Sigma.length, sc = 3;
    cv.width = n * sc; cv.height = n * sc;
    const ctx = cv.getContext('2d');
    for (let i = 0; i < n; i++) for (let j = 0; j < n; j++) {
      const v = Math.min(1, Math.abs(S.Sigma[i][j]));
      ctx.fillStyle = `rgb(${Math.round(253 - v * 214)},${Math.round(251 - v * 180)},${Math.round(247 - v * 79)})`;
      ctx.fillRect(j * sc, i * sc, sc, sc);
    }
    ctx.fillStyle = 'rgba(184,58,43,0.9)'; ctx.fillRect(0, S.sel * sc, n * sc, 1); ctx.fillRect(S.sel * sc, 0, 1, n * sc);
  },
  4() {
    const { act, E, picks } = S._cmp;
    const p = $('panel');
    p.innerHTML = `<p class="ptitle"><b>같은 ${act.length}개 센서, 네 가지 EfI</b> · Alt.1과 Alt.2는 완전히 겹칩니다</p><svg id="cmp" class="beam"></svg><div class="legend">${KEYS.map((k) => `<span><i style="background:${METHODS[k].color};height:${k === 'alt1' ? 4 : 2}px"></i>${METHODS[k].tag} ${METHODS[k].name}</span>`).join('')}</div>`;
    const svg = $('cmp'), H = 220, f = frame(svg, H), X = f.X, pt = 12, pb = 24;
    const hi = Math.max(...KEYS.flatMap((k) => Array.from(E[k]))) * 1.05;
    const Y = (v) => pt + (1 - v / hi) * (H - pt - pb);
    [0, 0.25, 0.5, 0.75, 1].filter((v) => v <= hi).forEach((v) => { el('line', { x1: X(0), y1: Y(v), x2: X(f.L), y2: Y(v), stroke: C.line }, svg); el('text', { x: X(0), y: Y(v) - 4 }, svg, v.toFixed(2)); });
    for (let x = 0; x <= f.L; x += 4) el('text', { x: X(x), y: H - 6, 'text-anchor': 'middle' }, svg, x);
    const xs = act.map((i) => S.beam.x[i]);
    ['alt1', 'alt4', 'alt3', 'alt2'].forEach((k) => {
      const m = METHODS[k];
      el('polyline', { points: xs.map((x, j) => `${X(x)},${Y(E[k][j])}`).join(' '), fill: 'none', stroke: m.color, 'stroke-width': k === 'alt1' ? 5 : 1.5, 'stroke-opacity': k === 'alt1' ? 0.45 : 1, 'stroke-linejoin': 'round' }, svg);
    });
    KEYS.forEach((k) => {
      const j = act.indexOf(picks[k]), x = X(xs[j]), yv = Y(E[k][j]);
      if (k === 'alt2') el('circle', { cx: x, cy: yv, r: 4, fill: C.ink }, svg);
      if ((k === 'alt3' || k === 'alt4') && picks[k] !== picks.alt2) pencil(svg, x, yv, 7);
    });
  },
  5() {
    const p = $('panel');
    const done = KEYS.filter((k) => S.res[k]);
    const best = Math.min(...done.map((k) => S.res[k].entropyByCount[S.N]));
    const tmax = Math.max(...done.map((k) => S.res[k].time));
    p.innerHTML = `<div class="pgrid"><div><p class="ptitle"><b>정보 엔트로피 −ln det F</b> · 오른쪽 81개에서 왼쪽으로 제거</p><svg id="ent" class="beam"></svg><div class="legend">${KEYS.map((k) => `<span><i style="background:${METHODS[k].color};height:${k === 'alt1' ? 4 : 2}px"></i>${METHODS[k].tag}</span>`).join('')}</div></div>
      <div><p class="ptitle"><b>N = ${S.N}</b> · 최선 대비</p><div class="rank">${KEYS.map((k) => { const r = S.res[k]; const v = r ? r.entropyByCount[S.N] : null;
        return `<div><span>${METHODS[k].tag} ${METHODS[k].name}</span><span class="v">${v == null ? '계산 중' : v - best < 1e-9 ? '최선' : '+' + (v - best).toFixed(3)}</span></div>`; }).join('')}</div>
      <p class="ptitle" style="margin-top:24px"><b>81 → 10 계산 시간</b></p><div class="rank">${KEYS.map((k) => { const r = S.res[k];
        return `<div><span>${METHODS[k].tag}</span><span class="v">${r ? fmtMs(r.time) : '계산 중'}</span><span class="t"><i style="width:${r ? Math.max(1, (r.time / tmax) * 100) : 0}%;background:${k === 'alt1' ? C.t2 : METHODS[k].color}"></i></span></div>`; }).join('')}</div></div></div>`;
    const svg = $('ent'), H = 220, f = frame(svg, H), W = f.W;
    const n0 = S.beam.nNode, t = target(), pl = 36, pr = 8, pt = 12, pb = 24;
    const vals = done.flatMap((k) => Object.values(S.res[k].entropyByCount));
    let lo = Math.min(...vals), hi = Math.max(...vals); const pad = (hi - lo) * 0.06; lo -= pad; hi += pad;
    const X = (n) => pl + ((n - t) / (n0 - t)) * (W - pl - pr), Y = (v) => pt + (1 - (v - lo) / (hi - lo)) * (H - pt - pb);
    for (let q = 0; q <= 4; q++) { const v = lo + ((hi - lo) * q) / 4; el('line', { x1: pl, y1: Y(v), x2: W - pr, y2: Y(v), stroke: C.line }, svg); el('text', { x: pl - 6, y: Y(v) + 4, 'text-anchor': 'end' }, svg, v.toFixed(0)); }
    for (let n = t; n <= n0; n += 10) el('text', { x: X(n), y: H - 6, 'text-anchor': 'middle' }, svg, n);
    el('line', { x1: X(S.N), y1: pt, x2: X(S.N), y2: H - pb, stroke: C.t3, 'stroke-dasharray': '2 3' }, svg);
    ['alt1', 'alt4', 'alt3', 'alt2'].forEach((k) => {
      const r = S.res[k]; if (!r) return;
      const pts = []; for (let n = t; n <= n0; n++) pts.push(`${X(n)},${Y(r.entropyByCount[n])}`);
      el('polyline', { points: pts.join(' '), fill: 'none', stroke: METHODS[k].color, 'stroke-width': k === 'alt1' ? 5 : 1.5, 'stroke-opacity': k === 'alt1' ? 0.45 : 1, 'stroke-linejoin': 'round' }, svg);
    });
    const hit = el('rect', { x: pl, y: pt, width: W - pl - pr, height: H - pt - pb, fill: 'transparent', style: 'cursor:crosshair' }, svg);
    hit.addEventListener('click', (ev) => { const b = svg.getBoundingClientRect(), px = ((ev.clientX - b.left) / b.width) * W; setN(Math.round(t + ((px - pl) / (W - pl - pr)) * (n0 - t))); });
  },
};
PANEL[6] = PANEL[5];

// ─────────────── 3센서 예제 (Table 1) ───────────────
const TOY = {
  ex1: { S: [[5, 3], [1, 0], [0, 1]], rho: 0, factor: 'perm' },
  ex2: { S: [[4, 5], [1, 0], [0, -1]], rho: 0.8, factor: 'sqrt' },
};
const toy = { ...TOY.ex1, key: 'ex1' };
const texMat = (A) => '\\begin{bmatrix}' + Array.from(A, (r) => Array.from(r, (v) => fx(v, 2).replace(/\.00$/, '')).join('&')).join('\\\\') + '\\end{bmatrix}';
function renderToy() {
  const Sm = fromArray(toy.S), Sg = fromArray([[1, 0, toy.rho], [0, 1, 0], [toy.rho, 0, 1]]);
  let L;
  if (toy.factor === 'chol') L = cholesky(invSPD(Sg));
  else { L = invSqrtSPD(Sg); if (toy.factor === 'perm') L = matmul(L, fromArray([[0, 1, 0], [1, 0, 0], [0, 0, 1]])); }
  $('toyMat').innerHTML = tex(`S=${texMat(Sm)}\\quad \\Sigma=${texMat(Sg)}\\quad L=${texMat(L)}`);
  const rows = [
    ['정의 1 − det F₋ᵢ/det F', efiDefinition(Sm, Sg)],
    ['제안식 Eq.(26)', efiProposed(Sm, Sg)],
    ['근사 diag(LᵀSF⁻¹SᵀL)', efiFactor(Sm, Sg, L)],
  ];
  const kT = argminTol(rows[0][1]);
  $('toyTable').innerHTML = '<tr><th>방법</th><th>e₁</th><th>e₂</th><th>e₃</th></tr>' + rows.map(([n, e], r) => {
    const k = argminTol(e), off = r > 0 && k !== kT;
    return `<tr class="${off ? 'off' : ''}"><td>${n}</td>${Array.from(e, (v, i) => `<td class="${i === k ? 'pick' : ''}${i === k && off ? ' circ' : ''}"><span>${fx(v)}</span></td>`).join('')}</tr>`;
  }).join('');
  const kA = argminTol(rows[2][1]);
  const same = rows[2][1].every((v, i) => Math.abs(v - rows[0][1][i]) < 1e-6);
  $('toyVerdict').innerHTML = kA !== kT
    ? `근사식은 센서 ${kA + 1}을 지우라고 하지만 기여가 가장 작은 센서는 센서 ${kT + 1}입니다. <b>순위가 뒤바뀌었습니다</b>`
    : same ? `이 경우엔 근사식도 값까지 정확합니다. ρ₁₃을 0이 아니게 하거나 L을 바꿔 보세요`
      : `근사식이 같은 센서를 고르긴 하지만 값은 정의와 다릅니다. 진하게 표시된 것이 각 방법의 최솟값입니다`;
}
function initToy() {
  const seg = $('toyPreset');
  seg.querySelectorAll('button').forEach((b) => b.addEventListener('click', () => {
    Object.assign(toy, TOY[b.dataset.v], { key: b.dataset.v });
    seg.querySelectorAll('button').forEach((x) => x.setAttribute('aria-pressed', x === b));
    $('toyRho').value = toy.rho; $('toyRhoOut').textContent = toy.rho.toFixed(2); $('toyFactor').value = toy.factor;
    renderToy();
  }));
  $('toyRho').addEventListener('input', (e) => { toy.rho = +e.target.value; $('toyRhoOut').textContent = toy.rho.toFixed(2); renderToy(); });
  $('toyFactor').addEventListener('change', (e) => { toy.factor = e.target.value; renderToy(); });
  $('toyFactor').value = toy.factor;
  renderToy();
}

// ─────────────── 렌더 및 상호작용 ───────────────
function meta() {
  const m = { eq29: 'Σ Eq.29', eq28: 'Σ Eq.28', identity: 'Σ = I' }[S.model];
  const same = S.res.alt1 && S.res.alt2 && S.res.alt1.removed.every((v, i) => v === S.res.alt2.removed[i]);
  $('meta').textContent = `${m}${S.model === 'identity' ? '' : ` · λ ${S.lambda.toFixed(1)}`} · 모드 ${S.Phi[0].length}개 · ${S.busy ? '정의식(Alt.1) 계산 중' : same ? 'Alt.1과 Alt.2 제거 순서 81→10 전부 일치' : '계산 완료'}`;
}
function render() {
  if (!S.res.alt2) return;
  document.querySelectorAll('.crumbs a').forEach((a) => a.classList.toggle('on', +a.dataset.ch === S.ch));
  FIG[S.ch](); PANEL[S.ch](); meta();
}
function setN(n) {
  S.N = Math.max(target(), Math.min(S.beam.nNode, n));
  $('nSens').value = S.N; $('nOut').textContent = S.N; S.hover = null;
  render();
}
function setCh(ch) { if (ch === S.ch) return; S.ch = ch; S.hover = null; render(); }

function initControls() {
  $('crumbs').innerHTML = [1, 2, 3, 4, 5, 6].map((c) => `<a href="#" data-ch="${c}" aria-label="${c}장">${String(c).padStart(2, '0')}</a>`).join('');
  $('crumbs').querySelectorAll('a').forEach((a) => a.addEventListener('click', (e) => {
    e.preventDefault();
    const ch = +a.dataset.ch, sec = document.querySelector(`.ch[data-ch="${ch}"]`);
    setCh(ch);
    if (window.matchMedia('(min-width: 1081px)').matches) sec.scrollIntoView({ behavior: 'smooth', block: 'start' });
  }));
  $('nSens').addEventListener('input', (e) => setN(+e.target.value));
  const seg = $('covModel');
  seg.querySelectorAll('button').forEach((b) => b.addEventListener('click', () => {
    S.model = b.dataset.v; seg.querySelectorAll('button').forEach((x) => x.setAttribute('aria-pressed', x === b));
    $('lambda').disabled = S.model === 'identity'; compute();
  }));
  let deb;
  $('lambda').addEventListener('input', (e) => { S.lambda = +e.target.value; $('lambdaOut').textContent = S.lambda.toFixed(1); clearTimeout(deb); deb = setTimeout(compute, 200); });
  $('play').addEventListener('click', togglePlay);
  window.addEventListener('resize', () => { clearTimeout(deb); deb = setTimeout(render, 120); });

  // 스크롤 위치로 챕터 전환
  const secs = [...document.querySelectorAll('.ch')];
  const io = new IntersectionObserver(() => {
    if (!window.matchMedia('(min-width: 1081px)').matches) return;
    const mid = window.innerHeight * 0.4;
    let best = secs[0];
    for (const s of secs) if (s.getBoundingClientRect().top < mid) best = s;
    secs.forEach((s) => s.classList.toggle('on', s === best));
    setCh(+best.dataset.ch);
  }, { threshold: [0, 0.25, 0.5, 0.75, 1], rootMargin: '-30% 0px -50% 0px' });
  secs.forEach((s) => io.observe(s));
  secs[0].classList.add('on');
}
function togglePlay() {
  const stop = () => { clearInterval(S.playing); S.playing = null; $('playLabel').textContent = '재생'; $('playIcon').setAttribute('d', 'M5 3.5v9l7-4.5z'); };
  if (S.playing) return stop();
  if (S.ch !== 5 && S.ch !== 6) setCh(5);
  let n = S.beam.nNode; setN(n);
  $('playLabel').textContent = '정지'; $('playIcon').setAttribute('d', 'M4.5 3.5h2.5v9h-2.5zM9 3.5h2.5v9h-2.5z');
  S.playing = setInterval(() => { n -= 1; if (n < target()) return stop(); setN(n); }, 140);
}

window.addEventListener('DOMContentLoaded', () => {
  setupBeam(); initControls(); initToy();
  $('nOut').textContent = S.N;
  const math = () => { if (window.renderMathInElement) renderMathInElement(document.querySelector('.story'), { delimiters: [{ left: '$$', right: '$$', display: true }, { left: '$', right: '$', display: false }], throwOnError: false }); renderToy(); };
  if (window.renderMathInElement) math(); else window.addEventListener('load', math);
  compute();
});
