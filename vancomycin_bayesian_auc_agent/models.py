"""
Clinical Data Models & Enums for VancoBayes Sentinel: Bayesian 24-Hour AUC24/MIC Pharmacokinetic Agent.
Domain: Clinical Pharmacokinetics
Standard: ASHP / IDSA / SIDP Vancomycin Guidelines
"""
import datetime
from dataclasses import dataclass, field
from enum import Enum
from typing import Dict, List, Optional, Any


class UrgencyLevel(str, Enum):
    ROUTINE = "ROUTINE"
    ADVISORY = "ADVISORY"
    WARNING = "WARNING"
    STAT_CRITICAL = "STAT_CRITICAL"


class ClinicalIntegrityStatus(str, Enum):
    CONCORDANT_NORMAL = "CONCORDANT_NORMAL"
    DISCORDANCE_DETECTED = "DISCORDANCE_DETECTED"
    CRITICAL_ACTION_REQUIRED = "CRITICAL_ACTION_REQUIRED"


@dataclass
class ClinicalCasePayload:
    case_id: str
    patient_synthetic_id: str
    primary_metric: float
    secondary_metric: float
    status_flag: str
    is_stat: bool = False
    clinical_notes: str = ""
    biomarkers: Dict[str, Any] = field(default_factory=dict)
    timestamp: str = field(default_factory=lambda: datetime.datetime.now(datetime.timezone.utc).isoformat())


@dataclass
class AgentAlert:
    alert_id: str
    sub_agent: str
    urgency: UrgencyLevel
    title: str
    clinical_finding: str
    actionable_recommendation: str
    guideline_citation: str = "ASHP / IDSA / SIDP Vancomycin Guidelines"
    timestamp: str = field(default_factory=lambda: datetime.datetime.now(datetime.timezone.utc).isoformat())

    def to_dict(self) -> Dict[str, Any]:
        return {
            "alert_id": self.alert_id,
            "sub_agent": self.sub_agent,
            "urgency": self.urgency.value,
            "title": self.title,
            "clinical_finding": self.clinical_finding,
            "actionable_recommendation": self.actionable_recommendation,
            "guideline_citation": self.guideline_citation,
            "timestamp": self.timestamp,
        }
