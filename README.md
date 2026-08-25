# Microsatellite Instability (MSI) Sensor Scorer

> **MSI/MMR Assessment for Immunotherapy Eligibility**
> Reference: Umar A et al. J Natl Cancer Inst. 2004;96(4):261-268 (Bethesda guidelines)

[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)
![Python](https://img.shields.io/badge/Python-3.10+-3776AB.svg)

## Overview

Real implementation of Microsatellite Instability (MSI) assessment for immunotherapy eligibility:

- **Bethesda panel** (5 markers): BAT-25, BAT-26, D2S123, D5S346, D17S250
- **MSI classification**: MSI-H (≥2 unstable), MSI-L (1 unstable), MSS (0 unstable)
- **MMR IHC interpretation**: MLH1, MSH2, MSH6, PMS2 with Lynch syndrome assessment
- **Combined MSI + MMR** assessment with concordance checking
- **MSI-H prevalence** by cancer type
- **Treatment implications**: Pembrolizumab, nivolumab, dostarlimab eligibility

## Quick Start

```bash
# MSI from Bethesda panel
python msi_scorer.py msi --bat25 unstable --bat26 unstable \
    --d2s123 stable --d5s346 stable --d17s250 stable

# MMR IHC interpretation
python msi_scorer.py mmr --mlh1 lost --msh2 retained --msh6 retained --pms2 lost \
    --braf present

# Combined assessment
python msi_scorer.py assess --bat25 unstable --bat26 unstable \
    --d2s123 stable --d5s346 stable --d17s250 stable \
    --mlh1 lost --msh2 retained --msh6 retained --pms2 lost --braf present

# MSI prevalence lookup
python msi_scorer.py prevalence --cancer-type colorectal

# Batch processing
python msi_scorer.py batch -i cases.csv -o results.csv
```

## Bethesda Panel Markers

| Marker | Type | Description |
|--------|------|-------------|
| BAT-25 | Mononucleotide | Poly(A) tract in c-KIT |
| BAT-26 | Mononucleotide | Poly(A) tract in MSH2 |
| D2S123 | Dinucleotide | Chromosome 2q |
| D5S346 | Dinucleotide | Chromosome 5q (APC region) |
| D17S250 | Dinucleotide | Chromosome 17q |

## MSI Classification

| Status | Criteria | Clinical Action |
|--------|----------|----------------|
| **MSI-H** | ≥2/5 unstable markers | Checkpoint inhibitor eligible |
| **MSI-L** | 1/5 unstable markers | Checkpoint inhibitors NOT indicated |
| **MSS** | 0/5 unstable markers | Checkpoint inhibitors NOT indicated |

## MMR IHC Interpretation

| MMR Loss Pattern | Likely Etiology | Lynch Syndrome? |
|------------------|-----------------|-----------------|
| MLH1 + PMS2 loss + BRAF V600E | Sporadic | No |
| MLH1 + PMS2 loss + methylation | Sporadic | No |
| MLH1 + PMS2 loss, no BRAF/methylation | Lynch syndrome | **Yes** |
| MSH2 + MSH6 loss | Lynch syndrome | **Yes** |
| Isolated MSH6 loss | Possible Lynch | **Yes** |
| Isolated PMS2 loss | Possible Lynch | **Yes** |

## MSI-H Prevalence by Cancer Type

| Cancer Type | MSI-H Prevalence |
|-------------|-----------------|
| Endometrial | 20-30% |
| Colorectal | 15% |
| Gastric | 15-20% |
| Ovarian | 10-12% |
| Hepatobiliary | 10-15% |
| Breast | 1-3% |
| Lung | 0.5-1% |

## Treatment Implications

- **MSI-H/dMMR (any solid tumor)**: Pembrolizumab (FDA-approved 2017)
- **MSI-H/dMMR colorectal cancer**: Nivolumab (FDA-approved)
- **dMMR solid tumors**: Dostarlimab (FDA-approved)

## Python API

```python
from msi_scorer import assess_msi_status, interpret_mmr_ihc

# MSI from Bethesda panel
markers = {"BAT-25": True, "BAT-26": True, "D2S123": False, "D5S346": False, "D17S250": False}
result = assess_msi_status(markers)
print(result["msi_status"])  # "MSI-H"

# MMR IHC
result = interpret_mmr_ihc(mlh1=False, msh2=True, msh6=True, pms2=False, braf_v600e=True)
print(result["mmr_status"])  # "dMMR"
print(result["etiology"])  # "Sporadic (MLH1 loss + BRAF V600E = sporadic hypermethylation)"
```

## Running Tests

```bash
python -m pytest test_msi_scorer.py -v
```

## License

MIT License. See [LICENSE](LICENSE) for details.
