#!/usr/bin/env python3
"""Tests for MSI Sensor Scorer - 20 real test cases."""
import pytest
from msi_scorer import (
    assess_msi_status, interpret_mmr_ihc, combined_msi_assessment,
    get_msi_prevalence, process_batch
)


# ---------------------------------------------------------------------------
# MSI Status Tests (Bethesda Panel)
# ---------------------------------------------------------------------------

class TestMSIStatus:

    def test_msih_two_unstable(self):
        markers = {
            "BAT-25": True, "BAT-26": True,
            "D2S123": False, "D5S346": False, "D17S250": False,
        }
        result = assess_msi_status(markers)
        assert result["msi_status"] == "MSI-H"
        assert result["unstable_count"] == 2

    def test_msih_all_unstable(self):
        markers = {
            "BAT-25": True, "BAT-26": True,
            "D2S123": True, "D5S346": True, "D17S250": True,
        }
        result = assess_msi_status(markers)
        assert result["msi_status"] == "MSI-H"
        assert result["unstable_count"] == 5

    def test_msih_three_unstable(self):
        markers = {
            "BAT-25": True, "BAT-26": True,
            "D2S123": True, "D5S346": False, "D17S250": False,
        }
        result = assess_msi_status(markers)
        assert result["msi_status"] == "MSI-H"

    def test_msil_one_unstable(self):
        markers = {
            "BAT-25": True, "BAT-26": False,
            "D2S123": False, "D5S346": False, "D17S250": False,
        }
        result = assess_msi_status(markers)
        assert result["msi_status"] == "MSI-L"
        assert result["unstable_count"] == 1

    def test_mss_zero_unstable(self):
        markers = {
            "BAT-25": False, "BAT-26": False,
            "D2S123": False, "D5S346": False, "D17S250": False,
        }
        result = assess_msi_status(markers)
        assert result["msi_status"] == "MSS"
        assert result["unstable_count"] == 0

    def test_msih_treatment_implications(self):
        markers = {
            "BAT-25": True, "BAT-26": True,
            "D2S123": False, "D5S346": False, "D17S250": False,
        }
        result = assess_msi_status(markers)
        assert "pembrolizumab" in result["treatment_implications"]
        assert "FDA-approved" in result["treatment_implications"]["pembrolizumab"]

    def test_mss_no_checkpoint(self):
        markers = {
            "BAT-25": False, "BAT-26": False,
            "D2S123": False, "D5S346": False, "D17S250": False,
        }
        result = assess_msi_status(markers)
        assert "NOT indicated" in result["treatment_implications"]["pembrolizumab"]

    def test_missing_markers(self):
        with pytest.raises(ValueError, match="Missing markers"):
            assess_msi_status({"BAT-25": True})

    def test_mono_vs_di_tracking(self):
        markers = {
            "BAT-25": True, "BAT-26": False,
            "D2S123": True, "D5S346": False, "D17S250": False,
        }
        result = assess_msi_status(markers)
        assert "BAT-25" in result["mononucleotide_unstable"]
        assert "D2S123" in result["dinucleotide_unstable"]


# ---------------------------------------------------------------------------
# MMR IHC Tests
# ---------------------------------------------------------------------------

class TestMMRIHC:

    def test_pmmr_all_retained(self):
        result = interpret_mmr_ihc(True, True, True, True)
        assert result["mmr_status"] == "pMMR"
        assert result["deficient_proteins"] == []
        assert result["suspected_lynch_syndrome"] is False

    def test_dmmr_mlh1_pms2_loss_sporadic(self):
        result = interpret_mmr_ihc(False, True, True, False, braf_v600e=True)
        assert result["mmr_status"] == "dMMR"
        assert "MLH1" in result["deficient_proteins"]
        assert "PMS2" in result["deficient_proteins"]
        assert "Sporadic" in result["etiology"]
        assert result["suspected_lynch_syndrome"] is False

    def test_dmmr_mlh1_pms2_loss_lynch(self):
        result = interpret_mmr_ihc(False, True, True, False,
                                    braf_v600e=False, mlh1_promoter_methylation=False)
        assert result["mmr_status"] == "dMMR"
        assert result["suspected_lynch_syndrome"] is True
        assert "Lynch" in result["etiology"]

    def test_dmmr_msh2_msh6_loss_lynch(self):
        result = interpret_mmr_ihc(True, False, False, True)
        assert result["mmr_status"] == "dMMR"
        assert result["suspected_lynch_syndrome"] is True
        assert "MSH2" in result["etiology"]

    def test_dmmr_isolated_msh6(self):
        result = interpret_mmr_ihc(True, True, False, True)
        assert result["mmr_status"] == "dMMR"
        assert result["suspected_lynch_syndrome"] is True

    def test_dmmr_isolated_pms2(self):
        result = interpret_mmr_ihc(True, True, True, False)
        assert result["mmr_status"] == "dMMR"
        assert result["suspected_lynch_syndrome"] is True

    def test_dmmr_methylation_confirmed(self):
        result = interpret_mmr_ihc(False, True, True, False,
                                    mlh1_promoter_methylation=True)
        assert "Sporadic" in result["etiology"]
        assert result["suspected_lynch_syndrome"] is False

    def test_pmmr_treatment(self):
        result = interpret_mmr_ihc(True, True, True, True)
        assert "NOT indicated" in result["treatment_implications"]["pembrolizumab"]


# ---------------------------------------------------------------------------
# Combined Assessment Tests
# ---------------------------------------------------------------------------

class TestCombinedAssessment:

    def test_combined_concordant_positive(self):
        markers = {
            "BAT-25": True, "BAT-26": True,
            "D2S123": False, "D5S346": False, "D17S250": False,
        }
        result = combined_msi_assessment(
            marker_results=markers,
            mlh1=False, msh2=True, msh6=True, pms2=False,
            braf_v600e=True,
        )
        assert result["concordance"] == "Concordant"
        assert result["final_status"] == "MSI-H/dMMR"

    def test_combined_pcr_only(self):
        markers = {
            "BAT-25": False, "BAT-26": False,
            "D2S123": False, "D5S346": False, "D17S250": False,
        }
        result = combined_msi_assessment(marker_results=markers)
        assert result["final_status"] == "MSS"
        assert result["mmr_ihc_result"] is None

    def test_combined_ihc_only(self):
        result = combined_msi_assessment(mlh1=True, msh2=True, msh6=True, pms2=True)
        assert result["final_status"] == "pMMR"
        assert result["msi_pcr_result"] is None

    def test_combined_no_data_raises(self):
        with pytest.raises(ValueError, match="At least one of"):
            combined_msi_assessment()


# ---------------------------------------------------------------------------
# Prevalence Tests
# ---------------------------------------------------------------------------

class TestPrevalence:

    def test_endometrial_prevalence(self):
        result = get_msi_prevalence("endometrial")
        assert result["msi_h_prevalence_approximate"] == 25.0

    def test_colorectal_prevalence(self):
        result = get_msi_prevalence("colorectal")
        assert result["msi_h_prevalence_approximate"] == 15.0

    def test_invalid_cancer_type(self):
        with pytest.raises(ValueError, match="Unknown cancer type"):
            get_msi_prevalence("unknown_cancer")


# ---------------------------------------------------------------------------
# Batch Processing Tests
# ---------------------------------------------------------------------------

class TestBatchProcessing:

    def test_batch_basic(self, tmp_path):
        csv_in = tmp_path / "in.csv"
        csv_out = tmp_path / "out.csv"
        csv_in.write_text(
            "BAT-25,BAT-26,D2S123,D5S346,D17S250\n"
            "true,true,false,false,false\n"
            "false,false,false,false,false\n",
            encoding="utf-8"
        )
        count = process_batch(str(csv_in), str(csv_out))
        assert count == 2
        content = csv_out.read_text(encoding="utf-8")
        assert "MSI-H" in content
        assert "MSS" in content
