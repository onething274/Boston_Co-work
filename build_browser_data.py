#!/usr/bin/env python3
"""build_browser_data.py — locuszoom browser MVP data extractor.

For 19 flip genes, dump per-gene JSON containing:
  - meta: chrom, tss, strand, gene_start/end, gene_id
  - window: TSS ± 250 kb
  - gwas       : [{pos, mlog10p, rsid}]
  - coloc      : [{pos, scp}]                (fastENLOC egfr_tube, SCP≥0.1)
  - dap        : [{pos, pip}]                (int_prior Tube, PIP any)
  - dar[CT]    : [{peak_id, start, end, center, logFC, p_val, p_val_adj, percentile}]
  - p2g[CT]    : [{peak_id, start, end, center, cor, pval, padj, percentile}]
  - exons      : [[start, end], ...]         (canonical transcript from gencode v28)

Output: /home/cklee/projects/Boston_Co-work/github/data/<GENE>.json
"""
from __future__ import annotations
import importlib.util
import json
from pathlib import Path

import numpy as np
import pandas as pd


F10_PATH = "/home/cklee/projects/Boston_Co-work/scripts/replication_scripts/260521/figure_scripts/F10_panelC_locuszoom.py"
_spec = importlib.util.spec_from_file_location("f10", F10_PATH)
f10 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(f10)

DATA_DIR = f10.DATA_DIR
OUT_DIR  = Path("/home/cklee/projects/Boston_Co-work/github/data")
WINDOW   = 250_000
SCP_THR  = 0.1

CTS = ["PT", "PT-VCAM1", "TAL", "DCT", "CNT", "PC", "ICA", "ICB", "ENDO"]
FLIP_GENES = f10.FLIP_GENES  # 19 genes


def load_p2g_perc(ct):
    """P2G — non-perc bed (has gene-peak mapping) merged with perc for percentile."""
    df = pd.read_csv(DATA_DIR / f"p2g/p2g_{ct}/twotail/p2g_{ct}_filt_overlap.hg38.bed",
                     sep="\t", usecols=["chrom", "start", "end", "peak_ID", "gene",
                                        "pval", "padj", "cor"])
    perc = pd.read_csv(DATA_DIR / f"p2g/p2g_{ct}/twotail/p2g_{ct}_filt_overlap_perc.hg38.bed",
                       sep="\t", usecols=["peak_ID", "percentile"])
    df = df.merge(perc, on="peak_ID", how="left")
    df["center"] = ((df["start"] + df["end"]) // 2).astype(int)
    return df


def load_dar_perc(ct):
    """DAR — perc bed (dedup peak universe with percentile)."""
    df = pd.read_csv(DATA_DIR / f"dar/dar_{ct}/twotail/dar_{ct}_filt_overlap_perc.hg38.bed",
                     sep="\t", usecols=["chrom", "start", "end", "peak_ID",
                                        "avg_logFC", "p_val", "p_val_adj", "percentile"])
    df["center"] = ((df["start"] + df["end"]) // 2).astype(int)
    return df


def clean_num(x):
    """JSON-safe conversion."""
    if x is None:
        return None
    if isinstance(x, (float, np.floating)):
        if np.isnan(x) or np.isinf(x):
            return None
        return float(x)
    if isinstance(x, (int, np.integer)):
        return int(x)
    return str(x)


def dar_rows(df, chrom, ws, we):
    sub = df[(df["chrom"] == chrom) & (df["center"] >= ws) & (df["center"] <= we)]
    out = []
    for r in sub.itertuples(index=False):
        out.append({
            "peak_id":    r.peak_ID,
            "start":      int(r.start),
            "end":        int(r.end),
            "center":     int(r.center),
            "logFC":      clean_num(r.avg_logFC),
            "p_val":      clean_num(r.p_val),
            "p_val_adj":  clean_num(r.p_val_adj),
            "percentile": clean_num(r.percentile),
        })
    return out


def p2g_rows(df, gene, ws, we):
    sub = df[(df["gene"] == gene) & (df["center"] >= ws) & (df["center"] <= we)]
    out = []
    for r in sub.itertuples(index=False):
        out.append({
            "peak_id":    r.peak_ID,
            "start":      int(r.start),
            "end":        int(r.end),
            "center":     int(r.center),
            "cor":        clean_num(r.cor),
            "pval":       clean_num(r.pval),
            "padj":       clean_num(r.padj),
            "percentile": clean_num(r.percentile),
        })
    return out


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    print(f"OUT_DIR = {OUT_DIR}")

    print("[1/5] Loading GENCODE ...")
    tss_df = f10.load_gencode_tss()

    print("[2/5] Loading fastENLOC + DAP + GWAS ...")
    enloc = f10.load_enloc()          # cols: chr, pos, gene, SCP, RCP
    gwas_by_c = f10.load_gwas_by_chrom()  # {chr: DataFrame[pos, nlp, RSID]}
    dap_tube = f10.load_dap_by_gene(f10.DAP_TUBE)  # {ensg: DataFrame[chr, pos, pip]}
    sym2ensg = f10.load_sym2ensg()

    print("[3/5] Loading peak data (9 CTs, WCC + DAR perc-adjusted) ...")
    p2g_by_ct = {}
    dar_by_ct = {}
    for ct in CTS:
        p2g_by_ct[ct] = load_p2g_perc(ct)
        dar_by_ct[ct] = load_dar_perc(ct)
        print(f"  {ct:9s}  p2g={len(p2g_by_ct[ct]):,}  dar={len(dar_by_ct[ct]):,}")

    print("[4/5] Loading canonical exons (gencode v28) ...")
    exons_by_gene = f10.load_canonical_exons(f10.GTF_FULL, FLIP_GENES)

    print("[5/5] Building per-gene JSON ...")
    manifest = []
    for gene in FLIP_GENES:
        if gene not in tss_df.index:
            print(f"  skip {gene} (no TSS)")
            continue
        chrom = tss_df.at[gene, "chrom"]
        tss   = int(tss_df.at[gene, "tss"])
        ws, we = tss - WINDOW, tss + WINDOW

        # GWAS
        gw = gwas_by_c.get(chrom)
        gwas_out = []
        if gw is not None:
            sub = gw[(gw["pos"] >= ws) & (gw["pos"] <= we)]
            for r in sub.itertuples(index=False):
                gwas_out.append({
                    "pos":     int(r.pos),
                    "mlog10p": clean_num(r.nlp),
                    "rsid":    str(r.RSID),
                })

        # coloc (fastENLOC SNPs for this gene ONLY)
        en_g = enloc[(enloc["gene"] == gene) & (enloc["chr"] == chrom)]
        en_g = en_g[(en_g["pos"] >= ws) & (en_g["pos"] <= we) & (en_g["SCP"] >= SCP_THR)]
        coloc_out = [{"pos": int(r.pos), "scp": clean_num(r.SCP)}
                     for r in en_g.itertuples(index=False)]

        # DAP (Tube int_prior — for tube-associated flip genes)
        ensg = sym2ensg.get(gene)
        dap_out = []
        if ensg is not None and ensg in dap_tube:
            d = dap_tube[ensg]
            d = d[(d["chr"] == chrom) & (d["pos"] >= ws) & (d["pos"] <= we)]
            for r in d.itertuples(index=False):
                dap_out.append({"pos": int(r.pos), "pip": clean_num(r.pip)})

        # DAR per CT (all peaks in window)
        dar_per_ct = {ct: dar_rows(dar_by_ct[ct], chrom, ws, we) for ct in CTS}
        # P2G per CT (only peaks assoc with this gene)
        p2g_per_ct = {ct: p2g_rows(p2g_by_ct[ct], gene, ws, we) for ct in CTS}

        gene_row = tss_df.loc[gene]
        exons = exons_by_gene.get(gene, {}).get("exons", [])

        rec = {
            "gene":       gene,
            "chrom":      str(chrom),
            "tss":        tss,
            "strand":     str(gene_row["strand"]),
            "gene_start": int(gene_row["start"]),
            "gene_end":   int(gene_row["end"]),
            "gene_id":    str(gene_row["gene_id"]),
            "window":     {"start": int(ws), "end": int(we)},
            "cts":        CTS,
            "gwas":       gwas_out,
            "coloc":      coloc_out,
            "dap":        dap_out,
            "dar":        dar_per_ct,
            "p2g":        p2g_per_ct,
            "exons":      [[int(a), int(b)] for a, b in exons],
        }

        out_fp = OUT_DIR / f"{gene}.json"
        with out_fp.open("w") as fh:
            json.dump(rec, fh, separators=(",", ":"))
        n_dar = sum(len(v) for v in dar_per_ct.values())
        n_p2g = sum(len(v) for v in p2g_per_ct.values())
        print(f"  {gene:10s} {chrom}:{ws:,}-{we:,}  gwas={len(gwas_out):,}  "
              f"coloc={len(coloc_out):,}  dap={len(dap_out):,}  "
              f"dar_all_ct={n_dar:,}  p2g_all_ct={n_p2g:,}  → {out_fp.name}")

        manifest.append({
            "gene":    gene,
            "chrom":   str(chrom),
            "tss":     tss,
            "file":    f"data/{gene}.json",
            "n_gwas":  len(gwas_out),
            "n_coloc": len(coloc_out),
        })

    with (OUT_DIR / "manifest.json").open("w") as fh:
        json.dump(manifest, fh, indent=2)
    print(f"saved: {OUT_DIR}/manifest.json ({len(manifest)} genes)")


if __name__ == "__main__":
    main()
