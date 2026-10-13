"""Regression tests for MSI data validation, discordance and audit integrity."""
import csv

import pytest

from msi_scorer import assess_msi_status, combined_msi_assessment, interpret_mmr_ihc, process_batch
from agents.base import AuditTrail

HEADERS = "BAT-25,BAT-26,D2S123,D5S346,D17S250\n"
VALUES = "true,true,false,false,false\n"
MARKERS = {"BAT-25": True, "BAT-26": True, "D2S123": False, "D5S346": False, "D17S250": False}


def test_unknown_marker_is_rejected():
    with pytest.raises(ValueError, match="Unknown markers"):
        assess_msi_status({**MARKERS, "EXTRA": True})


def test_ihc_rejects_non_boolean():
    with pytest.raises(ValueError, match="MLH1"):
        interpret_mmr_ihc("lost", True, True, False)


def test_partial_ihc_is_rejected():
    with pytest.raises(ValueError, match="all four"):
        combined_msi_assessment(marker_results=MARKERS, mlh1=True)


def test_discordant_pcr_ihc_not_collapsed_into_positive():
    result = combined_msi_assessment(marker_results=MARKERS, mlh1=True, msh2=True, msh6=True, pms2=True)
    assert result["concordance"] == "Discordant"
    assert result["final_status"] == "Discordant (PCR/IHC)"
    assert "review" in result["treatment_recommendations"]["note"].lower()


def test_invalid_batch_value_fails_closed(tmp_path):
    source = tmp_path / "input.csv"
    output = tmp_path / "out.csv"
    source.write_text(HEADERS + "yes,wrong,no,no,no\n")
    with pytest.raises(ValueError, match="Row 2"):
        process_batch(str(source), str(output))
    assert not output.exists()


def test_missing_batch_header_fails(tmp_path):
    source = tmp_path / "input.csv"
    source.write_text("BAT-25,BAT-26\ntrue,true\n")
    with pytest.raises(ValueError, match="Missing required"):
        process_batch(str(source), str(tmp_path / "out.csv"))


def test_input_output_cannot_be_same(tmp_path):
    source = tmp_path / "input.csv"
    source.write_text(HEADERS + VALUES)
    with pytest.raises(ValueError, match="differ"):
        process_batch(str(source), str(source))
    assert source.read_text() == HEADERS + VALUES


def test_discordant_batch_and_formula_protection(tmp_path):
    source = tmp_path / "input.csv"
    output = tmp_path / "out.csv"
    source.write_text(
        "BAT-25,BAT-26,D2S123,D5S346,D17S250,MLH1,MSH2,MSH6,PMS2,comment\n"
        "true,true,false,false,false,retained,retained,retained,retained,=2+2\n"
    )
    assert process_batch(str(source), str(output)) == 1
    with output.open(newline="") as fh:
        row = next(csv.DictReader(fh))
    assert row["final_status"] == "Discordant (PCR/IHC)"
    assert row["comment"] == "'=2+2"


def test_audit_hmac_detects_payload_tampering():
    audit = AuditTrail(secret_key="unit-test-key")
    audit.log("test", "worker", "PROCESS", {"value": "ok"})
    assert audit.verify_integrity()
    audit.logs[0]["payload_hash"] = "0" * 64
    assert not audit.verify_integrity()


def test_audit_hmac_detects_signature_tampering():
    audit = AuditTrail(secret_key="unit-test-key")
    audit.log("test", "worker", "PROCESS", {"value": "ok"})
    audit.logs[0]["current_hash"] = "0" * 64
    assert not audit.verify_integrity()
