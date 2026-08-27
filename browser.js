/* browser.js — MVP locuszoom browser.
 * Loads data/manifest.json → populates gene dropdown → on select, fetches
 * data/<gene>.json and renders a 5-row Plotly figure with shared x-axis. */

const $ = (id) => document.getElementById(id);
const CTS = ["PT", "PT-VCAM1", "TAL", "DCT", "CNT", "PC", "ICA", "ICB", "ENDO"];

const COLOR_SCP_STRONG = "#0B6E3B";
const COLOR_SCP_WEAK   = "#B8B8B8";
const COLOR_DAP        = "#4154A1";
const COLOR_GWAS       = "#8A8A8A";
const COLOR_TSS        = "#B22222";

async function fetchJSON(path) {
  const res = await fetch(path, { cache: "no-store" });
  if (!res.ok) throw new Error(`fetch ${path} → ${res.status}`);
  return res.json();
}

function fmt(n, d = 3) {
  if (n === null || n === undefined || Number.isNaN(n)) return "NA";
  const a = Math.abs(n);
  if (a !== 0 && (a < 1e-3 || a >= 1e4)) return n.toExponential(2);
  return Number(n).toFixed(d);
}

function getFilters() {
  const nn = (v, fb) => {
    const n = Number(v);
    return Number.isFinite(n) ? n : fb;
  };
  return {
    wccPct:  nn($("wcc-pct").value,  100),
    wccPadj: nn($("wcc-padj").value, 1),
    darPct:  nn($("dar-pct").value,  100),
    darPadj: nn($("dar-padj").value, 1),
  };
}

function applyFilters(g, f) {
  const dar = {}, p2g = {};
  let nDarKept = 0, nDarTotal = 0, nWccKept = 0, nWccTotal = 0;
  for (const ct of CTS) {
    const dArr = g.dar[ct] || [];
    nDarTotal += dArr.length;
    dar[ct] = dArr.filter(p =>
      (p.percentile === null || p.percentile <= f.darPct) &&
      (p.p_val_adj === null  || p.p_val_adj  <= f.darPadj));
    nDarKept += dar[ct].length;

    const pArr = g.p2g[ct] || [];
    nWccTotal += pArr.length;
    p2g[ct] = pArr.filter(p =>
      (p.percentile === null || p.percentile <= f.wccPct) &&
      (p.padj === null       || p.padj       <= f.wccPadj));
    nWccKept += p2g[ct].length;
  }
  const info =
    `WCC ${nWccKept.toLocaleString()}/${nWccTotal.toLocaleString()}   ·   ` +
    `DAR ${nDarKept.toLocaleString()}/${nDarTotal.toLocaleString()}`;
  $("filter-info").textContent = info;
  return { ...g, dar, p2g };
}

function buildFigure(g) {
  const win = g.window;
  const traces = [];

  // -------- Row 1: GWAS -log10(P) --------
  const gwX = g.gwas.map(d => d.pos);
  const gwY = g.gwas.map(d => d.mlog10p);
  const gwT = g.gwas.map(d => `pos ${d.pos.toLocaleString()}<br>-log10P ${fmt(d.mlog10p, 2)}<br>${d.rsid || ""}`);
  traces.push({
    type: "scattergl", mode: "markers",
    x: gwX, y: gwY, text: gwT, hoverinfo: "text",
    marker: { size: 4, color: COLOR_GWAS, opacity: 0.75 },
    name: "GWAS", xaxis: "x", yaxis: "y1",
  });

  // -------- Row 2: DAP PIP (colored by coloc SCP) --------
  const scpByPos = new Map(g.coloc.map(d => [d.pos, d.scp]));
  const dapX = g.dap.map(d => d.pos);
  const dapY = g.dap.map(d => d.pip);
  const dapColor = g.dap.map(d => {
    const s = scpByPos.get(d.pos);
    if (s === undefined) return COLOR_SCP_WEAK;
    return s >= 0.5 ? COLOR_SCP_STRONG : COLOR_SCP_WEAK;
  });
  const dapSize = g.dap.map(d => {
    const s = scpByPos.get(d.pos);
    if (s === undefined) return 4;
    return 6 + 6 * Math.min(1, s);
  });
  const dapT = g.dap.map(d => {
    const s = scpByPos.get(d.pos);
    return `pos ${d.pos.toLocaleString()}<br>PIP ${fmt(d.pip)}<br>SCP ${s === undefined ? "—" : fmt(s)}`;
  });
  traces.push({
    type: "scattergl", mode: "markers",
    x: dapX, y: dapY, text: dapT, hoverinfo: "text",
    marker: { size: dapSize, color: dapColor,
              line: { color: "black", width: g.dap.map(d => (scpByPos.get(d.pos) ?? 0) >= 0.5 ? 0.6 : 0) } },
    name: "DAP", xaxis: "x", yaxis: "y2",
  });

  // -------- Row 3: DAR log2FC per CT (strip) --------
  //   y = CT (categorical, from top to bottom = CTS[0]..[N-1] via reversed)
  //   x = center, color = log2FC (RdBu diverging), size = fixed
  const darX = [], darY = [], darC = [], darT = [];
  for (const ct of CTS) {
    for (const p of (g.dar[ct] || [])) {
      darX.push(p.center);
      darY.push(ct);
      darC.push(p.logFC ?? 0);
      darT.push(
        `CT ${ct}<br>peak ${p.peak_id}<br>pos ${p.center.toLocaleString()}<br>` +
        `logFC ${fmt(p.logFC)}<br>p_val ${fmt(p.p_val)}<br>` +
        `p_val_adj ${fmt(p.p_val_adj)}<br>pct ${fmt(p.percentile, 2)}`
      );
    }
  }
  const darMax = Math.max(1e-6, ...darC.map(Math.abs), 1);
  traces.push({
    type: "scattergl", mode: "markers",
    x: darX, y: darY, text: darT, hoverinfo: "text",
    marker: {
      size: 8, symbol: "square",
      color: darC, cmin: -darMax, cmax: darMax, colorscale: "RdBu", reversescale: true,
      colorbar: { title: "DAR<br>log2FC", x: 1.01, y: 0.53, len: 0.20, thickness: 8, tickfont: { size: 9 } },
    },
    name: "DAR", xaxis: "x", yaxis: "y3",
  });

  // -------- Row 4: WCC cor per CT (strip) --------
  const wccX = [], wccY = [], wccC = [], wccT = [];
  for (const ct of CTS) {
    for (const p of (g.p2g[ct] || [])) {
      wccX.push(p.center);
      wccY.push(ct);
      wccC.push(p.cor ?? 0);
      wccT.push(
        `CT ${ct}<br>peak ${p.peak_id}<br>pos ${p.center.toLocaleString()}<br>` +
        `cor ${fmt(p.cor)}<br>pval ${fmt(p.pval)}<br>` +
        `padj ${fmt(p.padj)}<br>pct ${fmt(p.percentile, 2)}`
      );
    }
  }
  const wccMax = Math.max(1e-6, ...wccC.map(Math.abs), 0.5);
  traces.push({
    type: "scattergl", mode: "markers",
    x: wccX, y: wccY, text: wccT, hoverinfo: "text",
    marker: {
      size: 8, symbol: "circle",
      color: wccC, cmin: -wccMax, cmax: wccMax, colorscale: "RdBu", reversescale: true,
      colorbar: { title: "WCC<br>cor", x: 1.01, y: 0.28, len: 0.20, thickness: 8, tickfont: { size: 9 } },
    },
    name: "WCC", xaxis: "x", yaxis: "y4",
  });

  // -------- Row 5: gene model (exons as bars + gene span line + TSS marker) --------
  //   Use one line for gene span, exons as thick segments.
  const geneLineY = 0;
  const gmTraces = [
    { // gene span
      type: "scatter", mode: "lines",
      x: [g.gene_start, g.gene_end], y: [geneLineY, geneLineY],
      line: { color: "#333333", width: 2 },
      hoverinfo: "skip", showlegend: false,
      xaxis: "x", yaxis: "y5",
    },
    { // TSS marker
      type: "scatter", mode: "markers+text",
      x: [g.tss], y: [geneLineY],
      text: ["TSS"], textposition: "top center",
      marker: { color: COLOR_TSS, size: 10, symbol: g.strand === "+" ? "triangle-right" : "triangle-left" },
      textfont: { size: 10, color: COLOR_TSS },
      hovertext: [`TSS ${g.tss.toLocaleString()}<br>strand ${g.strand}`],
      hoverinfo: "text", showlegend: false,
      xaxis: "x", yaxis: "y5",
    }
  ];
  // Exons as thick horizontal shapes via shapes list (rendered in layout below)
  const shapes = (g.exons || []).map(([s, e]) => ({
    type: "rect", xref: "x", yref: "y5",
    x0: s, x1: e, y0: -0.35, y1: 0.35,
    fillcolor: "#333333", line: { width: 0 }
  }));
  traces.push(...gmTraces);

  return { traces, shapes, gwYmax: Math.max(1, ...gwY.filter(v => !Number.isNaN(v))) };
}

function render(g) {
  const { traces, shapes, gwYmax } = buildFigure(g);
  const win = g.window;

  const layout = {
    margin: { t: 30, r: 90, l: 70, b: 40 },
    hovermode: "closest",
    showlegend: false,
    dragmode: "pan",
    height: 820,
    grid: { rows: 5, columns: 1, pattern: "independent" },
    // Row heights: GWAS, DAP, DAR, WCC, gene
    // Manual yaxis domain (top → bottom)
    xaxis: {
      title: `${g.chrom}  position`,
      range: [win.start, win.end],
      showgrid: true, gridcolor: "#eef1f4",
      zeroline: false,
      tickformat: ",",
      hoverformat: ",",
      matches: null,   // all yaxes share this x
    },
    // domains (top y1 down to y5)
    yaxis:  { domain: [0.82, 0.99], title: "-log10 P<br>(GWAS)",
              rangemode: "tozero", showgrid: true, gridcolor: "#eef1f4" },
    yaxis2: { domain: [0.63, 0.79], title: "PIP<br>(DAP)",
              range: [0, 1.05], showgrid: true, gridcolor: "#eef1f4" },
    yaxis3: { domain: [0.42, 0.60], title: "DAR", type: "category",
              categoryorder: "array", categoryarray: [...CTS].reverse(),
              showgrid: true, gridcolor: "#eef1f4" },
    yaxis4: { domain: [0.20, 0.38], title: "WCC", type: "category",
              categoryorder: "array", categoryarray: [...CTS].reverse(),
              showgrid: true, gridcolor: "#eef1f4" },
    yaxis5: { domain: [0.02, 0.14], title: "gene",
              range: [-1, 1], showgrid: false,
              zeroline: false, showticklabels: false },
    shapes,
    annotations: [{
      xref: "paper", yref: "paper", x: 0, y: 1.03,
      text: `<b>${g.gene}</b>  ${g.chrom}:${g.gene_start.toLocaleString()}-${g.gene_end.toLocaleString()}  strand ${g.strand}  ·  window ±250 kb  ·  ${g.gwas.length.toLocaleString()} GWAS SNPs · ${g.coloc.length} coloc SNP (SCP≥0.1) · ${g.dap.length.toLocaleString()} DAP SNPs`,
      showarrow: false, font: { size: 11 }
    }]
  };

  // Ensure x-axis matches across all 5 subplots
  ["y2", "y3", "y4", "y5"].forEach(y => {
    layout[`xaxis${y.slice(1)}`] = { matches: "x", showticklabels: false };
  });
  // xaxis (base) already covers row 1; give row 5 the label
  layout.xaxis5 = { matches: "x", showticklabels: true, title: `${g.chrom}  position` };
  layout.xaxis.showticklabels = false;
  layout.xaxis.title = "";

  Plotly.react("plot", traces, layout, {
    responsive: true,
    scrollZoom: true,
    displaylogo: false,
    modeBarButtonsToRemove: ["lasso2d", "select2d"],
  });
}

async function init() {
  const sel = $("gene-select");
  const info = $("gene-info");
  let manifest;
  try {
    manifest = await fetchJSON("data/manifest.json");
  } catch (e) {
    document.getElementById("plot").innerHTML =
      "<p style='color:#c00;padding:20px'>Failed to load data/manifest.json — build_browser_data.py 를 먼저 실행하세요.</p>";
    return;
  }
  manifest.sort((a, b) => a.gene.localeCompare(b.gene));
  for (const m of manifest) {
    const opt = document.createElement("option");
    opt.value = m.file;
    opt.textContent = `${m.gene}  (${m.chrom})`;
    sel.appendChild(opt);
  }

  let currentG = null;

  function rerender() {
    if (!currentG) return;
    render(applyFilters(currentG, getFilters()));
  }

  async function load(fp) {
    info.textContent = `loading ${fp} ...`;
    try {
      const g = await fetchJSON(fp);
      currentG = g;
      rerender();
      info.textContent = `${g.gene}  ·  ${g.chrom}:${g.gene_start.toLocaleString()}–${g.gene_end.toLocaleString()}  ·  TSS ${g.tss.toLocaleString()}  (${g.strand})`;
    } catch (e) {
      info.textContent = "load failed";
      console.error(e);
    }
  }

  sel.addEventListener("change", () => load(sel.value));
  for (const id of ["wcc-pct", "wcc-padj", "dar-pct", "dar-padj"]) {
    $(id).addEventListener("input", rerender);
  }
  $("reset-filters").addEventListener("click", () => {
    $("wcc-pct").value = 100;
    $("wcc-padj").value = 1;
    $("dar-pct").value = 100;
    $("dar-padj").value = 1;
    rerender();
  });

  load(manifest[0].file);
}

init();
