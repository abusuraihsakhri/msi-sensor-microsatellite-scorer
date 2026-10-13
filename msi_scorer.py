#!/usr/bin/env python3
"""
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
"""

import argparse
import csv
import json
import os
import sys
from pathlib import Path
from typing import Dict, Any, List, Optional


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

# Bethesda panel markers
BETHESDA_MARKERS = {
    "mononucleotide": ["BAT-25", "BAT-26"],
    "dinucleotide": ["D2S123", "D5S346", "D17S250"],
}

ALL_MARKERS = BETHESDA_MARKERS["mononucleotide"] + BETHESDA_MARKERS["dinucleotide"]

# MMR proteins
MMR_PROTEINS = ["MLH1", "MSH2", "MSH6", "PMS2"]

# MSI prevalence by cancer type
MSI_PREVALENCE = {
    "endometrial": {"range": "20-30%", "approximate": 25.0},
    "colorectal": {"range": "15%", "approximate": 15.0},
    "gastric": {"range": "15-20%", "approximate": 17.5},
    "ovarian": {"range": "10-12%", "approximate": 11.0},
    "hepatobiliary": {"range": "10-15%", "approximate": 12.5},
    "small_bowel": {"range": "5-10%", "approximate": 7.5},
    "urinary_tract": {"range": "5-10%", "approximate": 7.5},
    "glioblastoma": {"range": "5-10%", "approximate": 7.5},
    "breast": {"range": "1-3%", "approximate": 2.0},
    "prostate": {"range": "1-3%", "approximate": 2.0},
    "pancreatic": {"range": "1-2%", "approximate": 1.5},
    "lung": {"range": "0.5-1%", "approximate": 0.75},
}


# ---------------------------------------------------------------------------
# MSI Status Determination (Bethesda Panel)
# ---------------------------------------------------------------------------

def assess_msi_status(marker_results: Dict[str, bool]) -> Dict[str, Any]:
    """
    Assess MSI status using the Bethesda panel of 5 markers.

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
    """
    # Validate markers
    provided = set(marker_results.keys())
    expected = set(ALL_MARKERS)
    missing = expected - provided
    if missing:
        raise ValueError(f"Missing markers: {missing}. Required: {ALL_MARKERS}")

    unknown = provided - expected
    if unknown:
        raise ValueError(f"Unknown markers: {sorted(unknown)}. Allowed: {ALL_MARKERS}")

    for marker, status in marker_results.items():
        if not isinstance(status, bool):
            raise ValueError(f"Marker '{marker}' value must be bool; got {type(status)}")

    unstable_markers = [m for m in ALL_MARKERS if marker_results[m]]
    stable_markers = [m for m in ALL_MARKERS if not marker_results[m]]
    unstable_count = len(unstable_markers)
    total = len(ALL_MARKERS)
    percent = round((unstable_count / total) * 100, 1)

    # Classify by Bethesda criteria
    if unstable_count >= 2:
        msi_status = "MSI-H"
        description = (
            f"MSI-High: {unstable_count}/{total} markers unstable "
            f"({percent}%). >=2 unstable markers meets MSI-H criteria."
        )
    elif unstable_count == 1:
        msi_status = "MSI-L"
        description = (
            f"MSI-Low: {unstable_count}/{total} markers unstable "
            f"({percent}%). Only 1 unstable marker."
        )
    else:
        msi_status = "MSS"
        description = (
            f"Microsatellite Stable: {unstable_count}/{total} markers unstable "
            f"({percent}%). No instability detected."
        )

    # Categorize by marker type
    mono_unstable = [m for m in BETHESDA_MARKERS["mononucleotide"] if marker_results[m]]
    di_unstable = [m for m in BETHESDA_MARKERS["dinucleotide"] if marker_results[m]]

    result = {
        "msi_status": msi_status,
        "unstable_count": unstable_count,
        "total_markers": total,
        "unstable_markers": unstable_markers,
        "stable_markers": stable_markers,
        "percent_unstable": percent,
        "mononucleotide_unstable": mono_unstable,
        "dinucleotide_unstable": di_unstable,
        "description": description,
    }

    # Treatment implications
    if msi_status == "MSI-H":
        result["treatment_implications"] = {
            "pembrolizumab": "FDA-approved for MSI-H/dMMR solid tumors (any cancer type)",
            "nivolumab": "FDA-approved for MSI-H/dMMR colorectal cancer",
            "dostarlimab": "FDA-approved for dMMR solid tumors",
            "note": "MSI-H status is a strong predictor of response to immune checkpoint inhibitors.",
        }
    elif msi_status == "MSI-L":
        result["treatment_implications"] = {
            "pembrolizumab": "NOT indicated on the basis of MSI-L alone",
            "nivolumab": "NOT indicated for MSI-L tumors",
            "note": "MSI-L alone is not a validated positive predictive marker for checkpoint inhibitor response; assess other indications.",
        }
    else:
        result["treatment_implications"] = {
            "pembrolizumab": "NOT indicated on the basis of MSS alone",
            "nivolumab": "NOT indicated for MSS tumors",
            "note": "MSS alone does not determine overall immunotherapy eligibility; consider other biomarkers and indications.",
        }

    return result


# ---------------------------------------------------------------------------
# MMR IHC Interpretation
# ---------------------------------------------------------------------------

def interpret_mmr_ihc(mlh1: bool, msh2: bool, msh6: bool, pms2: bool,
                      braf_v600e: Optional[bool] = None,
                      mlh1_promoter_methylation: Optional[bool] = None) -> Dict[str, Any]:
    """
    Interpret Mismatch Repair (MMR) protein IHC results.

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
    """
    proteins = {"MLH1": mlh1, "MSH2": msh2, "MSH6": msh6, "PMS2": pms2}
    for protein, value in proteins.items():
        if not isinstance(value, bool):
            raise ValueError(f"{protein} expression must be bool")
    for name, value in {
        "braf_v600e": braf_v600e,
        "mlh1_promoter_methylation": mlh1_promoter_methylation,
    }.items():
        if value is not None and not isinstance(value, bool):
            raise ValueError(f"{name} must be bool or None")
    lost = [name for name, retained in proteins.items() if not retained]
    retained = [name for name, retained in proteins.items() if retained]

    if not lost:
        mmr_status = "pMMR"
        description = (
            "Proficient MMR (pMMR): All four MMR proteins are retained. "
            "Mismatch repair is intact."
        )
        etiology = "N/A - proficient repair"
        lynch = False
    else:
        mmr_status = "dMMR"
        description = (
            f"Deficient MMR (dMMR): Loss of expression of {', '.join(lost)}. "
            f"Retained: {', '.join(retained) if retained else 'none'}."
        )

        # Determine etiology
        if "MLH1" in lost and "PMS2" in lost:
            # MLH1 loss causes secondary PMS2 loss
            if braf_v600e is True:
                etiology = "Sporadic origin favored by BRAF V600E in colorectal tumors; correlate with methylation and clinical context"
                lynch = False
            elif mlh1_promoter_methylation is True:
                etiology = "Sporadic origin favored (MLH1 promoter methylation detected); rare constitutional causes require context"
                lynch = False
            elif braf_v600e is False and mlh1_promoter_methylation is False:
                etiology = "Lynch syndrome suspected (MLH1 loss without BRAF V600E or methylation)"
                lynch = True
            else:
                etiology = "Unknown etiology - further testing needed (BRAF/methylation)"
                lynch = None  # uncertain

        elif "MSH2" in lost and "MSH6" in lost:
            # MSH2 loss causes secondary MSH6 loss
            etiology = "Lynch syndrome suspected (MSH2 loss pattern)"
            lynch = True

        elif "MSH6" in lost and "MSH2" in retained:
            etiology = "Possible Lynch syndrome (isolated MSH6 loss)"
            lynch = True

        elif "PMS2" in lost and "MLH1" in retained:
            etiology = "Possible Lynch syndrome (isolated PMS2 loss)"
            lynch = True

        else:
            etiology = "Unusual MMR loss pattern - further evaluation recommended"
            lynch = None

    result = {
        "mmr_status": mmr_status,
        "deficient_proteins": lost,
        "retained_proteins": retained,
        "suspected_lynch_syndrome": lynch,
        "etiology": etiology,
        "description": description,
    }

    # Treatment implications
    if mmr_status == "dMMR":
        result["treatment_implications"] = {
            "pembrolizumab": "FDA-approved for dMMR/MSI-H solid tumors (any cancer type)",
            "nivolumab": "FDA-approved for dMMR/MSI-H colorectal cancer",
            "dostarlimab": "FDA-approved for dMMR solid tumors",
            "genetic_counseling": "Recommended if Lynch syndrome suspected",
        }
    else:
        result["treatment_implications"] = {
            "pembrolizumab": "NOT indicated on the basis of pMMR alone",
            "note": "pMMR alone does not determine overall immunotherapy eligibility; consider other biomarkers and indications.",
        }

    return result


# ---------------------------------------------------------------------------
# Combined MSI Assessment
# ---------------------------------------------------------------------------

def combined_msi_assessment(marker_results: Optional[Dict[str, bool]] = None,
                             mlh1: Optional[bool] = None,
                             msh2: Optional[bool] = None,
                             msh6: Optional[bool] = None,
                             pms2: Optional[bool] = None,
                             braf_v600e: Optional[bool] = None,
                             mlh1_promoter_methylation: Optional[bool] = None) -> Dict[str, Any]:
    """
    Combined MSI assessment using PCR-based markers and/or MMR IHC.

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
    """
    result = {
        "msi_pcr_result": None,
        "mmr_ihc_result": None,
        "concordance": None,
        "final_status": None,
    }

    # PCR-based MSI assessment
    if marker_results is not None:
        result["msi_pcr_result"] = assess_msi_status(marker_results)

    # IHC-based MMR assessment: do not silently discard partial panels.
    ihc_values = [mlh1, msh2, msh6, pms2]
    if any(v is not None for v in ihc_values) and not all(v is not None for v in ihc_values):
        raise ValueError("MMR IHC assessment requires all four proteins")
    if all(v is not None for v in ihc_values):
        result["mmr_ihc_result"] = interpret_mmr_ihc(
            mlh1, msh2, msh6, pms2,
            braf_v600e=braf_v600e,
            mlh1_promoter_methylation=mlh1_promoter_methylation,
        )

    # Determine final status
    pcr = result["msi_pcr_result"]
    ihc = result["mmr_ihc_result"]

    if pcr and ihc:
        pcr_msih = pcr["msi_status"] == "MSI-H"
        ihc_dmmr = ihc["mmr_status"] == "dMMR"
        if pcr_msih == ihc_dmmr:
            result["concordance"] = "Concordant"
        else:
            result["concordance"] = "Discordant"

        if result["concordance"] == "Discordant":
            result["final_status"] = "Discordant (PCR/IHC)"
        elif pcr_msih:
            result["final_status"] = "MSI-H/dMMR"
        elif pcr["msi_status"] == "MSI-L":
            result["final_status"] = "MSI-L/pMMR"
        else:
            result["final_status"] = "MSS/pMMR"

    elif pcr:
        result["final_status"] = pcr["msi_status"]
    elif ihc:
        result["final_status"] = ihc["mmr_status"]
    else:
        raise ValueError("At least one of marker_results or MMR IHC data must be provided")

    # Treatment recommendations
    if result["concordance"] == "Discordant":
        result["treatment_recommendations"] = {
            "note": "PCR and MMR IHC are discordant; review specimen, assay quality and clinical context before interpretation."
        }
    elif result["final_status"] in ("MSI-H", "MSI-H/dMMR", "dMMR"):
        result["treatment_recommendations"] = {
            "immune_checkpoint_inhibitors": "Indicated",
            "pembrolizumab": "FDA-approved (any solid tumor, MSI-H/dMMR)",
            "nivolumab": "FDA-approved (colorectal cancer, MSI-H/dMMR)",
            "genetic_counseling": "Recommended to evaluate for Lynch syndrome",
        }
    else:
        result["treatment_recommendations"] = {
            "immune_checkpoint_inhibitors": "NOT indicated solely from these MSI/MMR findings",
            "note": "Consider other biomarkers for treatment selection.",
        }

    return result


# ---------------------------------------------------------------------------
# MSI Prevalence Lookup
# ---------------------------------------------------------------------------

def get_msi_prevalence(cancer_type: str) -> Dict[str, Any]:
    """Get MSI-H prevalence for a specific cancer type."""
    key = cancer_type.lower().strip().replace(" ", "_").replace("-", "_")
    if key not in MSI_PREVALENCE:
        valid = ", ".join(sorted(MSI_PREVALENCE.keys()))
        raise ValueError(f"Unknown cancer type '{cancer_type}'. Valid: {valid}")

    info = MSI_PREVALENCE[key]
    return {
        "cancer_type": cancer_type,
        "msi_h_prevalence_range": info["range"],
        "msi_h_prevalence_approximate": info["approximate"],
    }


# ---------------------------------------------------------------------------
# Batch Processing
# ---------------------------------------------------------------------------

def _validate_safe_path(path_str: str, must_exist: bool = False) -> Path:
    """Normalize a caller-supplied local path; this trusted CLI is not a filesystem sandbox."""
    path = Path(path_str).resolve()
    # The local CLI caller explicitly controls the input and output paths.
    if must_exist and not path.exists():
        raise FileNotFoundError(f"File not found: {path_str}")
    if must_exist and not path.is_file():
        raise ValueError(f"Not a regular file: {path_str}")
    return path


def _parse_csv_bool(value: str, label: str, row_number: int, *, protein: bool = False) -> bool:
    """Accept only explicit stable/unstable or retained/lost values."""
    if value is None:
        raise ValueError(f"Row {row_number}: missing {label}")
    token = value.strip().lower()
    positives = {"true", "1", "yes", "retained", "positive", "present", "detected", "methylated"} if protein else {"true", "1", "yes", "unstable"}
    negatives = {"false", "0", "no", "lost", "negative", "absent", "not_detected", "unmethylated"} if protein else {"false", "0", "no", "stable"}
    if token in positives:
        return True
    if token in negatives:
        return False
    raise ValueError(f"Row {row_number}: invalid {label} value {value!r}")


def _csv_safe(value: str) -> str:
    """Prevent spreadsheet programs evaluating untrusted free-text cells as formulas."""
    if isinstance(value, str) and value.lstrip().startswith(("=", "+", "-", "@")):
        return "'" + value
    return value


def process_batch(input_csv: str, output_csv: str) -> int:
    """Classify a CSV of five-marker Bethesda panel observations.

    Requires all five marker headers. If any MMR IHC header is supplied,
    all four are required. Invalid and missing cell values fail closed.
    """
    in_path = _validate_safe_path(input_csv, must_exist=True)
    out_path = _validate_safe_path(output_csv, must_exist=False)
    if in_path == out_path:
        raise ValueError("Input and output CSV paths must differ")

    with open(in_path, mode="r", encoding="utf-8-sig", newline="") as f:
        reader = csv.DictReader(f)
        fieldnames = list(reader.fieldnames or [])
        if not fieldnames or len(fieldnames) != len(set(fieldnames)):
            raise ValueError("CSV must contain unique column names")
        missing = set(ALL_MARKERS) - set(fieldnames)
        if missing:
            raise ValueError(f"Missing required CSV columns: {sorted(missing)}")
        ihc_cols = {"MLH1", "MSH2", "MSH6", "PMS2"}
        if (ihc_cols & set(fieldnames)) and not ihc_cols.issubset(fieldnames):
            raise ValueError("MMR IHC CSV requires all four protein columns")
        result_columns = ["msi_status", "unstable_count", "mmr_status", "final_status", "treatment_note"]
        if set(result_columns) & set(fieldnames):
            raise ValueError("CSV contains reserved output column names")
        out_rows = []
        for row_number, row in enumerate(reader, start=2):
            if None in row:
                raise ValueError(f"Row {row_number}: more cells than CSV headers")
            markers = {marker: _parse_csv_bool(row.get(marker), marker, row_number)
                       for marker in ALL_MARKERS}
            kwargs = {}
            if ihc_cols.issubset(fieldnames):
                kwargs = {protein.lower(): _parse_csv_bool(row.get(protein), protein, row_number, protein=True)
                          for protein in MMR_PROTEINS}
                for field, arg in (("BRAF_V600E", "braf_v600e"),
                                   ("MLH1_PROMOTER_METHYLATION", "mlh1_promoter_methylation")):
                    value = row.get(field)
                    if value is not None and value.strip():
                        kwargs[arg] = _parse_csv_bool(value, field, row_number, protein=True)
            assessment = combined_msi_assessment(marker_results=markers, **kwargs)
            pcr = assessment["msi_pcr_result"]
            ihc = assessment["mmr_ihc_result"]
            item = {key: _csv_safe(value) for key, value in row.items()}
            item.update({
                "msi_status": pcr["msi_status"],
                "unstable_count": pcr["unstable_count"],
                "mmr_status": ihc["mmr_status"] if ihc else "N/A",
                "final_status": assessment["final_status"],
                "treatment_note": "Interpret with clinical context; not a treatment decision",
            })
            out_rows.append(item)

    with open(out_path, mode="w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames + result_columns)
        writer.writeheader()
        writer.writerows(out_rows)
    return len(out_rows)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main(argv=None):
    parser = argparse.ArgumentParser(
        description="Microsatellite Instability (MSI) Sensor Scorer"
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    # MSI status from Bethesda panel
    p_msi = subparsers.add_parser("msi", help="Assess MSI from Bethesda panel")
    p_msi.add_argument("--bat25", choices=["stable", "unstable"], required=True)
    p_msi.add_argument("--bat26", choices=["stable", "unstable"], required=True)
    p_msi.add_argument("--d2s123", choices=["stable", "unstable"], required=True)
    p_msi.add_argument("--d5s346", choices=["stable", "unstable"], required=True)
    p_msi.add_argument("--d17s250", choices=["stable", "unstable"], required=True)

    # MMR IHC
    p_mmr = subparsers.add_parser("mmr", help="Interpret MMR IHC results")
    p_mmr.add_argument("--mlh1", choices=["retained", "lost"], required=True)
    p_mmr.add_argument("--msh2", choices=["retained", "lost"], required=True)
    p_mmr.add_argument("--msh6", choices=["retained", "lost"], required=True)
    p_mmr.add_argument("--pms2", choices=["retained", "lost"], required=True)
    p_mmr.add_argument("--braf", choices=["present", "absent"], default=None)
    p_mmr.add_argument("--methylation", choices=["methylated", "unmethylated"], default=None)

    # Combined
    p_combined = subparsers.add_parser("assess", help="Combined MSI + MMR assessment")
    p_combined.add_argument("--bat25", choices=["stable", "unstable"], default=None)
    p_combined.add_argument("--bat26", choices=["stable", "unstable"], default=None)
    p_combined.add_argument("--d2s123", choices=["stable", "unstable"], default=None)
    p_combined.add_argument("--d5s346", choices=["stable", "unstable"], default=None)
    p_combined.add_argument("--d17s250", choices=["stable", "unstable"], default=None)
    p_combined.add_argument("--mlh1", choices=["retained", "lost"], default=None)
    p_combined.add_argument("--msh2", choices=["retained", "lost"], default=None)
    p_combined.add_argument("--msh6", choices=["retained", "lost"], default=None)
    p_combined.add_argument("--pms2", choices=["retained", "lost"], default=None)
    p_combined.add_argument("--braf", choices=["present", "absent"], default=None)
    p_combined.add_argument("--methylation", choices=["methylated", "unmethylated"], default=None)

    # Prevalence
    p_prev = subparsers.add_parser("prevalence", help="MSI-H prevalence by cancer type")
    p_prev.add_argument("--cancer-type", required=True)

    # Batch
    p_batch = subparsers.add_parser("batch", help="Batch process CSV")
    p_batch.add_argument("-i", "--input", required=True)
    p_batch.add_argument("-o", "--output", default="results.csv")

    # Serve (FastAPI)
    p_serve = subparsers.add_parser("serve", help="Start FastAPI REST server")
    p_serve.add_argument("--host", default="127.0.0.1")
    p_serve.add_argument("--port", type=int, default=8000)

    args = parser.parse_args(argv)

    if args.command == "msi":
        markers = {
            "BAT-25": args.bat25 == "unstable",
            "BAT-26": args.bat26 == "unstable",
            "D2S123": args.d2s123 == "unstable",
            "D5S346": args.d5s346 == "unstable",
            "D17S250": args.d17s250 == "unstable",
        }
        result = assess_msi_status(markers)
        print(json.dumps(result, indent=2))

    elif args.command == "mmr":
        braf = args.braf == "present" if args.braf else None
        meth = args.methylation == "methylated" if args.methylation else None
        result = interpret_mmr_ihc(
            args.mlh1 == "retained", args.msh2 == "retained",
            args.msh6 == "retained", args.pms2 == "retained",
            braf_v600e=braf, mlh1_promoter_methylation=meth,
        )
        print(json.dumps(result, indent=2))

    elif args.command == "assess":
        markers = None
        marker_args = ["bat25", "bat26", "d2s123", "d5s346", "d17s250"]
        if any(getattr(args, m) is not None for m in marker_args) and not all(
            getattr(args, m) is not None for m in marker_args
        ):
            parser.error("Combined assessment requires all five PCR marker arguments")
        if all(getattr(args, m) is not None for m in marker_args):
            markers = {
                "BAT-25": args.bat25 == "unstable",
                "BAT-26": args.bat26 == "unstable",
                "D2S123": args.d2s123 == "unstable",
                "D5S346": args.d5s346 == "unstable",
                "D17S250": args.d17s250 == "unstable",
            }
        mlh1 = args.mlh1 == "retained" if args.mlh1 else None
        msh2 = args.msh2 == "retained" if args.msh2 else None
        msh6 = args.msh6 == "retained" if args.msh6 else None
        pms2 = args.pms2 == "retained" if args.pms2 else None
        braf = args.braf == "present" if args.braf else None
        meth = args.methylation == "methylated" if args.methylation else None

        result = combined_msi_assessment(
            marker_results=markers,
            mlh1=mlh1, msh2=msh2, msh6=msh6, pms2=pms2,
            braf_v600e=braf, mlh1_promoter_methylation=meth,
        )
        print(json.dumps(result, indent=2))

    elif args.command == "prevalence":
        result = get_msi_prevalence(args.cancer_type)
        print(json.dumps(result, indent=2))

    elif args.command == "batch":
        count = process_batch(args.input, args.output)
        print(f"Processed {count} records -> {args.output}")

    elif args.command == "serve":
        import uvicorn
        from agents.api import app
        uvicorn.run(app, host=args.host, port=args.port)


if __name__ == "__main__":
    main()
