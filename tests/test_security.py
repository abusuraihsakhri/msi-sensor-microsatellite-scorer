"""
Security-focused tests for MSI Sensor Scorer.
Tests path traversal protection, HMAC key generation, and PHI guard enhancements.
"""
import os
import sys
import warnings
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).parent.parent))

from msi_scorer import process_batch, _validate_safe_path
from agents.base import AuditTrail, PHIGuard, SecurityException


class TestPathTraversalProtection:
    """Test that process_batch validates paths to prevent directory traversal."""

    def test_valid_input_file_works(self, tmp_path):
        """Normal file path should work."""
        csv_in = tmp_path / "in.csv"
        csv_out = tmp_path / "out.csv"
        csv_in.write_text(
            "BAT-25,BAT-26,D2S123,D5S346,D17S250\n"
            "true,false,false,false,false\n",
            encoding="utf-8"
        )
        count = process_batch(str(csv_in), str(csv_out))
        assert count == 1

    def test_missing_input_file_raises(self, tmp_path):
        """Non-existent input file should raise FileNotFoundError."""
        csv_in = tmp_path / "nonexistent.csv"
        csv_out = tmp_path / "out.csv"
        with pytest.raises(FileNotFoundError):
            process_batch(str(csv_in), str(csv_out))

    def test_directory_as_input_raises(self, tmp_path):
        """Directory as input should raise ValueError."""
        csv_out = tmp_path / "out.csv"
        with pytest.raises((FileNotFoundError, ValueError)):
            process_batch(str(tmp_path), str(csv_out))

    def test_validate_safe_path_resolves(self, tmp_path):
        """_validate_safe_path should resolve to absolute path."""
        test_file = tmp_path / "test.csv"
        test_file.write_text("test", encoding="utf-8")
        result = _validate_safe_path(str(test_file), must_exist=True)
        assert result.is_absolute()
        assert result.exists()


class TestHMACKeyGeneration:
    """Test secure HMAC key handling."""

    def test_audit_trail_generates_key_when_env_unset(self):
        """When AUDIT_SECRET_KEY is not set, a random key should be generated."""
        # Ensure env var is not set
        old_val = os.environ.pop("AUDIT_SECRET_KEY", None)
        try:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", RuntimeWarning)
                trail = AuditTrail()
            # Key should be set (not empty)
            assert len(trail.secret_key) > 0
        finally:
            if old_val is not None:
                os.environ["AUDIT_SECRET_KEY"] = old_val

    def test_audit_trail_uses_env_var_when_set(self):
        """When AUDIT_SECRET_KEY is set, it should be used."""
        old_val = os.environ.get("AUDIT_SECRET_KEY")
        os.environ["AUDIT_SECRET_KEY"] = "test-secret-key-12345"
        try:
            trail = AuditTrail()
            assert trail.secret_key == b"test-secret-key-12345"
        finally:
            if old_val is not None:
                os.environ["AUDIT_SECRET_KEY"] = old_val
            else:
                del os.environ["AUDIT_SECRET_KEY"]

    def test_audit_trail_accepts_explicit_key(self):
        """Explicit key parameter should take precedence."""
        trail = AuditTrail(secret_key="my-explicit-key")
        assert trail.secret_key == b"my-explicit-key"

    def test_different_ephemeral_keys(self):
        """Two trails without env var should have different keys."""
        old_val = os.environ.pop("AUDIT_SECRET_KEY", None)
        try:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", RuntimeWarning)
                trail1 = AuditTrail()
                trail2 = AuditTrail()
            assert trail1.secret_key != trail2.secret_key
        finally:
            if old_val is not None:
                os.environ["AUDIT_SECRET_KEY"] = old_val


class TestPHIGuardEnhancements:
    """Test PHI guard pattern matching."""

    def test_mrn_pattern_detected(self):
        with pytest.raises(SecurityException):
            PHIGuard.assert_no_phi("Patient MRN-12345678")

    def test_ssn_pattern_detected(self):
        with pytest.raises(SecurityException):
            PHIGuard.assert_no_phi("SSN: 123-45-6789")

    def test_phone_pattern_detected(self):
        with pytest.raises(SecurityException):
            PHIGuard.assert_no_phi("Call 555-123-4567")

    def test_email_pattern_detected(self):
        with pytest.raises(SecurityException):
            PHIGuard.assert_no_phi("Email: patient@hospital.org")

    def test_clean_text_passes(self):
        """Non-PHI text should pass without exception."""
        PHIGuard.assert_no_phi("Specimen KEY-001 analytical result optimal")
        PHIGuard.assert_no_phi("BAT-25 marker stable, MSI assessment nominal")

    def test_redact_phi(self):
        """PHI should be redacted."""
        result = PHIGuard.redact_phi("Patient MRN-12345678 is stable")
        assert "MRN-12345678" not in result
        assert "REDACTED_IDENTIFIER" in result

    def test_empty_string_passes(self):
        """Empty string should pass."""
        PHIGuard.assert_no_phi("")
        PHIGuard.assert_no_phi(None)
