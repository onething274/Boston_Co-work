# Locuszoom browser (MVP)

19 flip gene 각각에 대해 TSS ±250 kb 창의 locuszoom 을 브라우저에서 슬라이드/hover 확인.

라이브 데모: https://onething274.github.io/Boston_Co-work/  (GitHub Pages 활성화 필요)

## 구성

```
.
├── index.html            # gene dropdown + 5-track Plotly figure
├── browser.js            # 렌더 로직
├── style.css
├── build_browser_data.py # 데이터 추출 스크립트
└── data/
    ├── manifest.json     # gene 목록
    └── <GENE>.json       # gene 별 데이터 (19 개)
```

## 5 track (top → bottom)
1. **GWAS eGFR** — Wuttke 2019 raw p, y = −log10 P
2. **DAP PIP + coloc SCP** — Tube int_prior. SCP ≥ 0.5 = 짙은 초록 + 검은 테두리, SCP < 0.5 는 회색, DAP-only 는 작음 (build 시 coloc SCP ≥ 0.1 만 포함)
3. **DAR log2FC (per CT strip)** — 9 tube CT (`PT, PT-VCAM1, TAL, DCT, CNT, PC, ICA, ICB, ENDO`), 색 = log2FC (RdBu diverging). Hover → peak_id, log2FC, p_val, p_val_adj, percentile
4. **WCC cor (per CT strip)** — 9 CT, 색 = cor. Hover → peak_id, cor, pval, padj, percentile
5. **Gene model** — canonical transcript (gencode v28, APPRIS principal) exon + TSS 삼각형

## 데이터 소스
- GWAS: `data/260521/gene_coloc/twas_raw_data/20171016_MW_eGFR_overall_ALL_nstud61.hg38.txt.gz`
- DAP: `data/260521/gene_coloc/ptwas_raw_data/int_prior/DAP_snp_PIPS_Tube_hg38.txt`
- Coloc: `webpeace/.../egfr_tube.enloc.snp.out` (SCP ≥ 0.1)
- DAR: `data/260521/dar/dar_<CT>/twotail/dar_<CT>_filt_overlap_perc.hg38.bed`
- P2G: `data/260521/p2g/p2g_<CT>/twotail/p2g_<CT>_filt_overlap.hg38.bed` (+ perc merge)
- Exon: `data/reference/gencode.v28.annotation.gtf.gz`

## 데이터 재생성
```bash
conda activate cklee_py
python build_browser_data.py
```
19 gene 처리 ~1-2 분, 총 산출물 ~21 MB.

## 로컬 실행
`file://` 는 fetch 를 못 쓰므로 static 서버 필요:

```bash
python -m http.server 8000
```
그리고 브라우저에서 `http://localhost:8000/` (또는 SSH tunneling 시 `-L 8000:localhost:8000`).

## GitHub Pages 활성화
1. Repo → Settings → Pages
2. Source: **Deploy from a branch**
3. Branch: `main` / Folder: `/ (root)` → Save
4. 몇 분 뒤 위 라이브 데모 URL 로 접근

## 인터랙션
- 상단 dropdown 으로 gene 선택
- 마우스 드래그 = 팬, 스크롤 = 줌 (Plotly toolbar 옵션)
- Hover 로 각 point 상세값
- x축은 5 track 이 shared → 한 곳에서 zoom 하면 전부 같이 이동

## MVP 제한
- LD r² 색 미포함 (GWAS 점은 회색 단일). 필요하면 후속 확장
- Track 하드코딩, CT toggle / p-value slider 없음
- 26.08.26 미팅 확정된 F17-style rank dot + top20 fill + `*` 스펙 아직 반영 안 함
