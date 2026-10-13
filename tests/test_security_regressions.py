import copy
from agents.base import AuditTrail, AuditLogger, SecurityException, PHIGuard
from agents.models import SystemTaskPayload
from agents.supervisor import SystemSupervisor


def test_hmac_verifies_signature_not_only_chain():
    ledger = AuditTrail(secret_key="test-secret-only")
    ledger.log("tester", "test", "EVENT", {"status": "SUCCESS"})
    assert ledger.verify_integrity()
    ledger.logs[0]["event_type"] = "FORGED"
    assert not ledger.verify_integrity()


def test_hmac_detects_bad_signature():
    ledger = AuditTrail(secret_key="test-secret-only")
    ledger.log("tester", "test", "EVENT", {"status": "SUCCESS"})
    ledger.logs[0]["current_hash"] = "0" * 64
    assert not ledger.verify_integrity()


def test_audit_reader_cannot_modify_internal_records():
    ledger = AuditTrail(secret_key="test-secret-only")
    ledger.log("tester", "test", "EVENT", {"status": "SUCCESS"})
    external = ledger.get_trail()
    external[0]["event_type"] = "TAMPERED"
    assert ledger.verify_integrity()
    assert ledger.get_trail()[0]["event_type"] == "EVENT"


def test_nested_identifiers_are_checked():
    supervisor = SystemSupervisor(model_provider="mock")
    payload = SystemTaskPayload(
        task_id="SAFE-1",
        target_identifier="KEY-1",
        primary_metric=12.0,
        attributes={"notes": "MRN-123456"}
    )
    import pytest
    with pytest.raises(SecurityException):
        supervisor.process_task(payload)


def test_public_audit_endpoint_does_not_disclose_trail():
    from agents.api import api_audit_logs
    result = api_audit_logs()
    assert set(result) == {"entry_count", "verified"}
    assert result["verified"]
