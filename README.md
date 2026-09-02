# MSI Sensor Microsatellite Scorer

> **Domain:** Medical Oncology & Cancer Staging Systems  
> **Reference Guidelines & Standards:** `AJCC Cancer Staging Manual & NCCN Clinical Practice Guidelines`

<div align="center">

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
![Python](https://img.shields.io/badge/Python-3.10%20%7C%203.11%20%7C%203.12-3776AB.svg?logo=python&logoColor=white)
![FastAPI](https://img.shields.io/badge/FastAPI-0.111-009688.svg?logo=fastapi&logoColor=white)
![Audit Trail](https://img.shields.io/badge/Audit-HMAC--SHA256_Tamper--Evident-brightgreen.svg)
![Zero-PHI Guard](https://img.shields.io/badge/Guard-Zero--PHI_Outbound-blue.svg)
![Docker](https://img.shields.io/badge/Docker-Ready-2496ED.svg?logo=docker&logoColor=white)

</div>

---

## 📖 What It Does

Microsatellite Instability (MSI) Sensor Scorer

Implements MSI assessment using the Bethesda panel and MMR IHC interpretation
for microsatellite instability determination.

References:
  - Umar A et al. J Natl Cancer Inst. 2004;96(4):261-268 (Bethesda guidelines)
  - Boland CR et al. Cancer Res. 1998;58(22):5248-5257 (Bethesda panel)
  - Le DT et al. Science. 2015;350(6257):191-196 (MSI-H and immunotherapy)
  - FDA Approval: Pembrolizumab for MSI-H/dMMR solid tumors (2017)

Zero-dependency Python implementation (stdlib only).
License: MIT

---

## ⚙️ Key Capabilities & Algorithmic Modules

### 🔬 Analytical Functions

- **`assess_msi_status()`**: Assess MSI status using the Bethesda panel of 5 markers.

Parameters
----------
marker_results : dict
    Dictionary mapping marker names to instability status.
    Keys should be from: BAT-25, BAT-26, D2S123, D5S346, D17S250
    Values: True = unstable, False = stable

Returns
-------
dict with msi_status, unstable_count, total_markers, unstable_markers,
stable_markers, percent_unstable, description, treatment_implications
- **`interpret_mmr_ihc()`**: Interpret Mismatch Repair (MMR) protein IHC results.

Parameters
----------
mlh1 : bool
    MLH1 expression (True = retained, False = lost).
msh2 : bool
    MSH2 expression (True = retained, False = lost).
msh6 : bool
    MSH6 expression (True = retained, False = lost).
pms2 : bool
    PMS2 expression (True = retained, False = lost).
braf_v600e : bool, optional
    BRAF V600E mutation status (True = present, False = absent).
mlh1_promoter_methylation : bool, optional
    MLH1 promoter methylation (True = methylated, False = unmethylated).

Returns
-------
dict with mmr_status, deficient_proteins, retained_proteins,
suspected_lynch_syndrome, etiology, description, treatment_implications
- **`combined_msi_assessment()`**: Combined MSI assessment using PCR-based markers and/or MMR IHC.

Parameters
----------
marker_results : dict, optional
    Bethesda panel marker results.
mlh1, msh2, msh6, pms2 : bool, optional
    MMR protein IHC results.
braf_v600e : bool, optional
    BRAF V600E mutation status.
mlh1_promoter_methylation : bool, optional
    MLH1 promoter methylation status.

Returns
-------
dict with msi_pcr_result, mmr_ihc_result, concordance, final_status,
treatment_recommendations
- **`get_msi_prevalence()`**: Get MSI-H prevalence for a specific cancer type.
- **`process_batch()`**: Process a CSV of MSI cases.

Expected columns: BAT-25, BAT-26, D2S123, D5S346, D17S250 (as true/false)
Optional: MLH1, MSH2, MSH6, PMS2, BRAF_V600E

---

## 💻 CLI Quickstart & Usage

### 1. Guided Interactive Mode
```bash
python cli.py
```

### 2. Direct Parameterized Evaluation
```bash
python cli.py --input data.csv
```

### Parameter Reference
- `--interactive`: Launch guided terminal interactive wizard.
- `--input <path>`: Evaluate input from JSON or CSV specification.
- `--json`: Output deterministic structured results in JSON format.

### Input Data Schema

| Field | Description | Requirement |
|:------|:------------|:------------|
| `Patient_ID` | Parameter / observation metric | Required |
| `v1` | Parameter / observation metric | Required |
| `v2` | Parameter / observation metric | Required |
| `v3` | Parameter / observation metric | Required |

---

## 🛡️ Security & Enterprise Architecture

* **Zero-PHI Outbound Interceptor:** Active AST and regex inspection blocking SSNs, MRNs, phone numbers, and patient identifiers.
* **Tamper-Evident HMAC-SHA256 Audit Trail:** Chained, cryptographically signed logs for every evaluation and state transition.
* **Air-Gapped LLM Reasoning Adapter:** Agnostic integration for local Ollama instances (`llama3`, `mistral`), Claude 3.5 Sonnet, GPT-4o, and deterministic test mocks.
* **Active Learning Bayesian Calibration:** Dynamic tracker updating worker reliability weights and monitoring Brier calibration drift.
* **FastAPI & Prometheus Telemetry:** Exposes OpenAPI 3.1 REST endpoints and operational Prometheus metrics (`/metrics`).

---

## 🧪 Testing & Verification

Run the automated test suite:

```bash
pytest -v
```

Execute high-throughput batch simulation benchmarks:

```bash
python simulator.py --tasks 1000 --concurrency 8
```

---

## 🐳 Container Deployment

```bash
docker build -t msi-sensor-microsatellite-scorer .
docker run -p 8000:8000 msi-sensor-microsatellite-scorer
```
