// GNU Affero General Public License v3.0 — Copyright (c) 2026 Fair
// SPDX-License-Identifier: AGPL-3.0

// ============================================================
//  greens.js — 现代克拉尼板（点激励）的精确物理模型
//  依据论文《Exploring the Origin of Maximum Entropy States
//  Relevant to Resonant Modes in Modern Chladni Plates》
//  (Shu, Tseng, Lai, Yu, Huang & Chen, Entropy 2022, 24, 215)
//
//  核心结论：现代克拉尼板由「点激励器」驱动（如音源/震源置于板心），
//  其共振节线图样对应「非齐次亥姆霍兹方程的格林函数」的
//  *最大熵态*（Eqs. 4–15），而非经典自由板本征函数本身。
//
//  —— 经典 Chladni（本仓库旧模型）：ψ = cos(mπx)cos(nπy) ∓ cos(nπx)cos(mπy)
//  —— 现代 Chladni（本文件）：        Ψ(r;k) = Σ_n ψ_n(r_s)·ψ_n(r)/(k²−k_n²)
//     其中 ψ_n 为齐次亥姆霍兹本征函数（诺依曼边界，自由边近似），
//     r_s 为点激励位置（论文取板心），k 取「最大熵态」对应的共振波数。
//
//  本模块为纯数学（无 DOM 依赖），浏览器与 Node 均可运行。
// ============================================================

const PI = Math.PI;

// 阻尼（洛伦兹展宽），以波数 k 的比例给出，避免 k=k_n 处分母发散。
// 取值 ~2% 即可在「单模共振」与「最大熵展宽」之间取得平衡。
const DEFAULT_GAMMA = 0.02;

// ------------------------------------------------------------
//  构建基底：齐次亥姆霍兹本征函数 {ψ_n(r)} 及其本征值 k_n、
//  在激励点 r_s 的耦合强度 coupling_n = norm²·ψ_n(r_s)。
//  每个基底项用 (fx, fy) 记录空间频率，使位移场可统一写成
//    Ψ(r;k) = Σ_j coupling_j/(k² − k_j² + i·2γk) · cos(fx_j π u) cos(fy_j π v)
//  （u=x/L, v=y/L ∈ [0,1]）。方板 fx=n, fy=m；三角板每支 (n,m) 展开为 3 项。
// ------------------------------------------------------------
export function buildBasis(
  shape,
  N = 100,
  opts = {},
) {
  const L = 1; // 归一化边长，实际尺度由 freqToMode 的物理参数承载
  const gamma = opts.gamma ?? DEFAULT_GAMMA;
  const terms = [];

  if (shape === "triangle") {
    // 等边三角形：顶点 (0,0) (L/2, √3L/2) (−L/2, √3L/2)
    // 激励点取重心 r_s = (0, √3L/3)
    const us = 0; // x_s/L
    const vs = Math.sqrt(3) / 3; // y_s/L = 1/√3
    const norm = Math.sqrt(16 / (Math.sqrt(3) * L * L)); // = √(16/(√3 L²))
    const norm2 = norm * norm;
    for (let n = 1; n <= N; n++) {
      for (let m = 2 * n; m <= N; m++) {
        // 论文 m ≥ 2n，避免三角对称重复计数
        // 本征值（Eq 17）：k̂ = (4π/(√3 L))√(n²+m²−nm)
        const kHat = (4 * PI / (Math.sqrt(3) * L)) * Math.sqrt(n * n + m * m - n * m);
        const kHat2 = kHat * kHat;
        // 三项（Eq 16）的空间频率（以 πu 为单位）：fx = 2(..)/3
        const comps = [
          [2 * (2 * n - m) / 3, 2 * m / 3],
          [2 * (2 * m - n) / 3, 2 * n / 3],
          [2 * (n + m) / 3, 2 * (n - m) / 3],
        ];
        // 三项在激励点的取值之和（论文 S_A(r_s)+S_B(r_s)+S_C(r_s)）
        let sumExc = 0;
        const sub = [];
        for (const [fx, fy] of comps) {
          const vExc = Math.cos(fx * PI * us) * Math.cos(fy * PI * vs);
          sumExc += vExc;
          sub.push([fx, fy, vExc]);
        }
        // coupling = norm²·Σ S_i(r_s)（论文 a_n ∝ ψ_n(r_s)）
        const coupling = norm2 * sumExc;
        // 三项共享同一本征值与本征函数耦合，各自空间频率不同
        for (const [fx, fy] of sub) {
          terms.push({
            fx,
            fy,
            kHat2,
            coupling, // 注意：三项耦合相同（同一 mode 的 ψ_n(r_s)）
          });
        }
      }
    }
  } else {
    // 正方形（默认）：区域 0≤x,y≤L，诺依曼边界 → 余弦本征函数
    //   ψ_{n,m}(r) = (2/L) cos(nπx/L) cos(mπy/L)         （Eq 12）
    //   k_{n,m}   = (π/L)√(n²+m²)                        （Eq 13）
    // 激励点取板心 r_s=(L/2,L/2)：cos(nπ/2)cos(mπ/2)，
    // 故中心激励只激发偶-偶模态（奇 n/m 耦合为 0）。
    const xus = 0.5; // x_s/L
    const yvs = 0.5; // y_s/L
    const norm2 = (2 / L) * (2 / L); // (2/L)²
    for (let n = 0; n <= N; n++) {
      for (let m = 0; m <= N; m++) {
        if (n === 0 && m === 0) continue; // 平凡模（k=0）跳过
        const kHat2 = (PI / L) * (PI / L) * (n * n + m * m);
        const vExc = Math.cos(n * PI * xus) * Math.cos(m * PI * yvs);
        const coupling = norm2 * vExc;
        terms.push({
          fx: n,
          fy: m,
          kHat2,
          coupling,
        });
      }
    }
  }

  // 按本征值升序排序，便于扫描共振与剪枝
  terms.sort((a, b) => a.kHat2 - b.kHat2);
  const kMax = Math.sqrt(terms[terms.length - 1].kHat2);
  return {
    shape,
    N,
    L,
    gamma,
    terms,
    kMax,
  };
}

// ------------------------------------------------------------
//  给定波数 k，计算每项实系数 coeff_j = Re[ coupling_j / (k² − k_j² + i·2γk) ]
//  并剪枝到幅值最大的 topK 项（近共振时主导模态只有少数几项）。
//  返回 { terms:[{fx,fy,coeff}], k, peak }
// ------------------------------------------------------------
export function buildSuperposition(
  basis,
  k,
  topK = 48,
  opts = {},
) {
  const gamma = opts.gamma ?? basis.gamma ?? DEFAULT_GAMMA;
  const k2 = k * k;
  const broaden = 2 * gamma * k; // 虚部系数
  const broaden2 = broaden * broaden;
  const all = [];
  for (const t of basis.terms) {
    const denomRe = k2 - t.kHat2;
    const denomMag2 = denomRe * denomRe + broaden2;
    if (denomMag2 < 1e-12) continue;
    // 复系数 coupling/(denomRe + i·broaden)，取实部
    const coeff = (t.coupling * denomRe) / denomMag2;
    if (!isFinite(coeff)) continue;
    all.push({
      fx: t.fx,
      fy: t.fy,
      coeff,
      mag: Math.abs(coeff),
    });
  }
  all.sort((a, b) => b.mag - a.mag);
  const trimmed = all.slice(0, topK);
  let peak = 1e-6;
  for (const t of trimmed) peak = Math.max(peak, Math.abs(t.coeff));
  return {
    k,
    terms: trimmed.map((t) => ({
      fx: t.fx,
      fy: t.fy,
      coeff: t.coeff / peak * 2, // 峰值归一化到 ≈2，沿用经典量纲
    })),
    peak: 2,
  };
}

// 由模式 (M,N) 反推目标波数 k（方板 k=π√(M²+N²)/L；三角板用 Eq 17 同形估算）
export function modeToK(
  shape,
  M,
  N,
  L = 1,
) {
  if (shape === "triangle") {
    return (4 * PI / (Math.sqrt(3) * L)) * Math.sqrt(M * M + N * N - M * N);
  }
  return (PI / L) * Math.sqrt(M * M + N * N);
}

// ------------------------------------------------------------
//  信息熵 N_eff(k) = exp(S)，S = −Σ p_n ln p_n
//  p_n = |a_n|² / Σ|a_n|²，a_n = coupling_n/(k²−k_n²)   （Eq 10/11/14）
//  共振频率 = N_eff 的局部极大值（最大熵态）。
// ------------------------------------------------------------
export function neff(
  basis,
  k,
  opts = {},
) {
  const gamma = opts.gamma ?? basis.gamma ?? DEFAULT_GAMMA;
  const k2 = k * k;
  const broaden = 2 * gamma * k;
  const broaden2 = broaden * broaden;
  let sumP = 0;
  const ps = [];
  for (const t of basis.terms) {
    const denomRe = k2 - t.kHat2;
    const denomMag2 = denomRe * denomRe + broaden2;
    const p = (t.coupling * t.coupling) / denomMag2;
    ps.push(p);
    sumP += p;
  }
  if (sumP < 1e-30) return 1;
  let S = 0;
  for (const p of ps) {
    const pn = p / sumP;
    if (pn > 1e-12) S -= pn * Math.log(pn);
  }
  return Math.exp(S);
}

// 扫描 k ∈ [kStart,kEnd]，返回 N_eff 的局部极大值（共振点）
export function findResonances(
  basis,
  kStart,
  kEnd,
  opts = {},
) {
  const step = opts.step ?? (kEnd - kStart) / 4000;
  const gamma = opts.gamma ?? basis.gamma ?? DEFAULT_GAMMA;
  const peaks = [];
  let prev = neff(basis, kStart, { gamma });
  let prevSlope = 0;
  let prevK = kStart;
  for (let k = kStart + step; k < kEnd; k += step) {
    const cur = neff(basis, k, { gamma });
    const slope = cur - prev;
    if (prevSlope > 0 && slope <= 0 && prev > 1.05) {
      peaks.push({ k: prevK, neff: prev });
    }
    prevSlope = slope;
    prev = cur;
    prevK = k;
  }
  return peaks;
}

// 在已排序基底中，找离 targetK 最近的「被激励」模态本征值（coupling≠0）
export function nearestExcitedK(
  basis,
  targetK,
) {
  let best = null;
  let bestD = Infinity;
  for (const t of basis.terms) {
    if (Math.abs(t.coupling) < 1e-9) continue;
    const d = Math.abs(Math.sqrt(t.kHat2) - targetK);
    if (d < bestD) {
      bestD = d;
      best = Math.sqrt(t.kHat2);
    }
  }
  return best ?? targetK;
}

// 找离 targetK 最近的 N_eff 局部极大值（最大熵共振）
export function nearestResonance(
  basis,
  targetK,
  opts = {},
) {
  const span = opts.span ?? Math.max(targetK * 0.5, 2);
  const kStart = Math.max(0.01, targetK - span);
  const kEnd = targetK + span;
  const peaks = findResonances(basis, kStart, kEnd, opts);
  if (peaks.length === 0) {
    // 区间内无峰：回退到最近被激励模态
    return nearestExcitedK(basis, targetK);
  }
  let best = peaks[0];
  let bestD = Infinity;
  for (const p of peaks) {
    const d = Math.abs(p.k - targetK);
    if (d < bestD) {
      bestD = d;
      best = p;
    }
  }
  return best.k;
}

// ------------------------------------------------------------
//  计算位移场网格 Ψ(u,v)（实场），用于渲染与粒子物理。
//  默认在 (M,N) 对应的目标波数附近取「最大熵共振」k，
//  使图样即论文所述现代克拉尼板节线。
//  返回 { field:Float32Array(gridN*gridN), gridN, k, peak }
// ------------------------------------------------------------
export function resonantField(
  basis,
  M,
  N,
  opts = {},
) {
  const gridN = opts.gridN ?? 256;
  const topK = opts.topK ?? 48;
  const gamma = opts.gamma ?? basis.gamma ?? DEFAULT_GAMMA;
  const targetK = modeToK(basis.shape, M, N, basis.L);
  const k = opts.snap === false
    ? targetK
    : nearestResonance(basis, targetK, { gamma, span: opts.span });

  const sup = buildSuperposition(basis, k, topK, { gamma });
  const field = new Float32Array(gridN * gridN);
  let idx = 0;
  for (let j = 0; j < gridN; j++) {
    const v = (j + 0.5) / gridN;
    const cym = [];
    for (const t of sup.terms) cym.push(Math.cos(t.fy * PI * v));
    for (let i = 0; i < gridN; i++) {
      const u = (i + 0.5) / gridN;
      let s = 0;
      for (let q = 0; q < sup.terms.length; q++) {
        s += sup.terms[q].coeff * Math.cos(sup.terms[q].fx * PI * u) * cym[q];
      }
      field[idx++] = s;
    }
  }
  return {
    field,
    gridN,
    k,
    terms: sup.terms,
    peak: sup.peak,
  };
}
