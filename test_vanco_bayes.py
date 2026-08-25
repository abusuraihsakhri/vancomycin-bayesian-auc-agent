"""
Tests for Vancomycin AUC24/MIC Bayesian Dosing Calculator.
"""
import math
import pytest
from vanco_bayes import (
    calculate_crcl_cockcroft_gault,
    calculate_vancomycin_clearance,
    calculate_volume_of_distribution,
    calculate_elimination_rate_constant,
    calculate_half_life,
    calculate_auc24_intermittent,
    calculate_trough_steady_state,
    calculate_peak_steady_state,
    calculate_auc24_from_trough,
    calculate_loading_dose,
    adjust_dose_for_target_auc,
    interpret_auc_mic,
    interpret_trough,
    estimate_ke_from_two_levels,
    estimate_vd_from_level,
    bayesian_estimate,
    full_dosing_assessment,
    main,
    AUC_TARGET_MIN,
    AUC_TARGET_MAX,
    TROUGH_SERIOUS_MIN,
    TROUGH_SERIOUS_MAX,
)


# ============================================================================
# CrCl Tests
# ============================================================================

class TestCrCl:
    def test_normal_male(self):
        # ((140-65)*80)/(72*1.2) = 75*80/86.4 = 69.44
        crcl = calculate_crcl_cockcroft_gault(80, 65, 1.2, False)
        assert abs(crcl - 69.44) < 0.1

    def test_normal_female(self):
        crcl = calculate_crcl_cockcroft_gault(80, 65, 1.2, True)
        assert abs(crcl - 69.44 * 0.85) < 0.1

    def test_young_healthy(self):
        crcl = calculate_crcl_cockcroft_gault(70, 30, 0.9, False)
        expected = ((140 - 30) * 70) / (72 * 0.9)
        assert abs(crcl - expected) < 0.1

    def test_elderly(self):
        crcl = calculate_crcl_cockcroft_gault(60, 85, 1.8, False)
        expected = ((140 - 85) * 60) / (72 * 1.8)
        assert abs(crcl - expected) < 0.1

    def test_invalid_scr(self):
        with pytest.raises(ValueError):
            calculate_crcl_cockcroft_gault(70, 50, 0, False)

    def test_invalid_weight(self):
        with pytest.raises(ValueError):
            calculate_crcl_cockcroft_gault(0, 50, 1.0, False)


# ============================================================================
# Vancomycin Clearance Tests
# ============================================================================

class TestVancomycinClearance:
    def test_matzke_equation(self):
        # CL = CrCl * 0.689 + 3.66
        cl = calculate_vancomycin_clearance(100)
        expected = 100 * 0.689 + 3.66
        assert abs(cl - expected) < 0.01

    def test_zero_crcl(self):
        cl = calculate_vancomycin_clearance(0)
        assert abs(cl - 3.66) < 0.01

    def test_negative_crcl(self):
        with pytest.raises(ValueError):
            calculate_vancomycin_clearance(-10)


# ============================================================================
# Volume of Distribution Tests
# ============================================================================

class TestVolumeOfDistribution:
    def test_normal(self):
        vd = calculate_volume_of_distribution(80)
        assert abs(vd - 56.0) < 0.1  # 80 * 0.7

    def test_custom_vd_per_kg(self):
        vd = calculate_volume_of_distribution(80, vd_per_kg=0.9)
        assert abs(vd - 72.0) < 0.1

    def test_obesity_adjusted(self):
        vd = calculate_volume_of_distribution(
            100, obesity_adjusted=True,
            actual_weight_kg=100, ideal_body_weight_kg=70
        )
        # Vd = 0.7 * 70 + 0.4 * (100 - 70) = 49 + 12 = 61
        assert abs(vd - 61.0) < 0.1

    def test_invalid_weight(self):
        with pytest.raises(ValueError):
            calculate_volume_of_distribution(0)


# ============================================================================
# Ke and Half-life Tests
# ============================================================================

class TestKeAndHalfLife:
    def test_ke_calculation(self):
        # Ke = (CL_ml_min * 60 / 1000) / Vd
        ke = calculate_elimination_rate_constant(72.56, 56.0)
        expected = (72.56 * 60 / 1000) / 56.0
        assert abs(ke - expected) < 0.001

    def test_half_life(self):
        t_half = calculate_half_life(0.02)
        assert abs(t_half - 34.65) < 0.1  # 0.693 / 0.02

    def test_invalid_ke(self):
        with pytest.raises(ValueError):
            calculate_half_life(0)


# ============================================================================
# AUC24 Tests
# ============================================================================

class TestAUC24:
    def test_standard_calculation(self):
        # AUC24 = (Dose / CL_L_per_h) * (24 / tau)
        # CL_L_per_h = 72.56 * 60 / 1000 = 4.3536
        auc = calculate_auc24_intermittent(1000, 72.56, 12)
        expected = (1000 / (72.56 * 60 / 1000)) * (24 / 12)
        assert abs(auc - expected) < 1.0

    def test_different_intervals(self):
        auc_q24 = calculate_auc24_intermittent(1500, 72.56, 24)
        auc_q12 = calculate_auc24_intermittent(750, 72.56, 12)
        # Both should give same AUC24
        assert abs(auc_q24 - auc_q12) < 1.0

    def test_invalid_dose(self):
        with pytest.raises(ValueError):
            calculate_auc24_intermittent(0, 72.56, 12)

    def test_invalid_cl(self):
        with pytest.raises(ValueError):
            calculate_auc24_intermittent(1000, 0, 12)


# ============================================================================
# Trough Tests
# ============================================================================

class TestTrough:
    def test_steady_state_trough(self):
        trough = calculate_trough_steady_state(1000, 56.0, 0.0154, 12)
        assert trough > 0
        assert trough < 50  # Reasonable range

    def test_trough_increases_with_longer_interval(self):
        t_short = calculate_trough_steady_state(1000, 56.0, 0.0154, 8)
        t_long = calculate_trough_steady_state(1000, 56.0, 0.0154, 24)
        # Longer interval = lower trough
        assert t_short > t_long


# ============================================================================
# AUC from Trough Tests
# ============================================================================

class TestAUCFromTrough:
    def test_estimation(self):
        auc = calculate_auc24_from_trough(15.0, 0.0154, 12)
        assert auc > 0
        assert auc > 100  # Should be reasonable

    def test_invalid_trough(self):
        with pytest.raises(ValueError):
            calculate_auc24_from_trough(-1, 0.0154, 12)


# ============================================================================
# Loading Dose Tests
# ============================================================================

class TestLoadingDose:
    def test_standard(self):
        ld = calculate_loading_dose(25.0, 56.0)
        # LD = 56 * 25 / 1.0 = 1400, rounded to 1500
        assert ld == 1500

    def test_rounding(self):
        ld = calculate_loading_dose(20.0, 50.0)
        # LD = 50 * 20 = 1000
        assert ld == 1000

    def test_invalid_vd(self):
        with pytest.raises(ValueError):
            calculate_loading_dose(25.0, 0)


# ============================================================================
# Dose Adjustment Tests
# ============================================================================

class TestDoseAdjustment:
    def test_increase_for_low_auc(self):
        new_dose = adjust_dose_for_target_auc(1000, 300, 500)
        # New = 1000 * (500/300) = 1666.67 -> 1750
        assert new_dose > 1000

    def test_decrease_for_high_auc(self):
        new_dose = adjust_dose_for_target_auc(1000, 700, 500)
        # New = 1000 * (500/700) = 714.29 -> 750
        assert new_dose < 1000

    def test_no_change_at_target(self):
        new_dose = adjust_dose_for_target_auc(1000, 500, 500)
        assert new_dose == 1000

    def test_invalid_measured_auc(self):
        with pytest.raises(ValueError):
            adjust_dose_for_target_auc(1000, 0, 500)


# ============================================================================
# Interpretation Tests
# ============================================================================

class TestInterpretation:
    def test_auc_mic_therapeutic(self):
        result = interpret_auc_mic(500, 1.0)
        assert result["status"] == "THERAPEUTIC"
        assert result["auc_mic_ratio"] == 500.0

    def test_auc_mic_subtherapeutic(self):
        result = interpret_auc_mic(300, 1.0)
        assert result["status"] == "SUBTHERAPEUTIC"

    def test_auc_mic_supratherapeutic(self):
        result = interpret_auc_mic(700, 1.0)
        assert result["status"] == "SUPRATHERAPEUTIC"

    def test_trough_therapeutic(self):
        result = interpret_trough(17.5, "serious")
        assert result["status"] == "THERAPEUTIC"

    def test_trough_subtherapeutic(self):
        result = interpret_trough(8.0, "serious")
        assert result["status"] == "SUBTHERAPEUTIC"

    def test_trough_supratherapeutic(self):
        result = interpret_trough(25.0, "serious")
        assert result["status"] == "SUPRATHERAPEUTIC"


# ============================================================================
# Two-Level Ke Estimation Tests
# ============================================================================

class TestKeFromTwoLevels:
    def test_estimation(self):
        # Ke = ln(30/15) / 12 = 0.693/12 = 0.0578
        ke = estimate_ke_from_two_levels(30.0, 15.0, 12.0)
        assert abs(ke - 0.0578) < 0.001

    def test_invalid_levels(self):
        with pytest.raises(ValueError):
            estimate_ke_from_two_levels(15.0, 30.0, 12.0)  # Second higher


# ============================================================================
# Bayesian Estimation Tests
# ============================================================================

class TestBayesian:
    def test_basic_estimation(self):
        result = bayesian_estimate(
            measured_trough=15.0,
            measured_time_hours=12.0,
            dose_mg=1000,
            tau_hours=12.0,
            weight_kg=80,
            crcl_ml_min=70.0
        )
        assert "bayesian_vd" in result
        assert "bayesian_ke" in result
        assert "auc24" in result
        assert result["bayesian_vd"] > 0
        assert result["bayesian_ke"] > 0


# ============================================================================
# Full Assessment Tests
# ============================================================================

class TestFullAssessment:
    def test_basic_assessment(self):
        result = full_dosing_assessment(
            weight_kg=80, age_years=65, serum_creatinine_mg_dl=1.2,
            is_female=False, current_dose_mg=1000, tau_hours=12
        )
        assert "patient_parameters" in result
        assert "pk_parameters" in result
        assert "auc24_analysis" in result
        assert result["patient_parameters"]["crcl_ml_min"] > 0
        assert result["pk_parameters"]["vd_liters"] > 0

    def test_assessment_with_trough(self):
        result = full_dosing_assessment(
            weight_kg=80, age_years=65, serum_creatinine_mg_dl=1.2,
            is_female=False, current_dose_mg=1000, tau_hours=12,
            measured_trough=18.0
        )
        assert result["bayesian_estimation"] is not None
        assert result["dose_recommendation"] is not None

    def test_assessment_with_low_trough(self):
        result = full_dosing_assessment(
            weight_kg=80, age_years=65, serum_creatinine_mg_dl=1.2,
            is_female=False, current_dose_mg=750, tau_hours=12,
            measured_trough=8.0
        )
        assert result["dose_recommendation"]["action"] == "INCREASE"

    def test_assessment_with_high_trough(self):
        result = full_dosing_assessment(
            weight_kg=80, age_years=65, serum_creatinine_mg_dl=1.2,
            is_female=False, current_dose_mg=1500, tau_hours=12,
            measured_trough=25.0
        )
        assert result["dose_recommendation"]["action"] == "DECREASE"


# ============================================================================
# CLI Tests
# ============================================================================

class TestCLI:
    def test_crcl_command(self):
        ret = main(["crcl", "--weight", "80", "--age", "65", "--scr", "1.2"])
        assert ret == 0

    def test_auc_command(self):
        ret = main(["auc", "--dose", "1000", "--cl", "72.56", "--interval", "12"])
        assert ret == 0

    def test_loading_dose_command(self):
        ret = main(["loading-dose", "--weight", "80", "--target", "25"])
        assert ret == 0

    def test_assess_command(self):
        ret = main(["assess", "--weight", "80", "--age", "65", "--scr", "1.2",
                     "--dose", "1000", "--interval", "12"])
        assert ret == 0

    def test_bayesian_command(self):
        ret = main(["bayesian", "--trough", "15.0", "--dose", "1000",
                     "--interval", "12", "--weight", "80", "--crcl", "70"])
        assert ret == 0
