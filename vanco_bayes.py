#!/usr/bin/env python3
"""
Vancomycin AUC24/MIC Bayesian Dosing Calculator

Real pharmacokinetic calculations for vancomycin therapeutic drug monitoring.
Uses population PK parameters with Bayesian-informed individualization.

Key formulas:
- AUC24 = Dose / (CL * tau)  [intermittent infusion]
- CL = CrCl * 0.689 + 3.66  [Matzke equation]
- Vd = 0.7 L/kg (population), adjusted for patient factors
- Ke = CL / Vd
- t1/2 = 0.693 / Ke
- Trough estimation: Ctrough = (Dose/Vd) * (e^(-Ke*tau) / (1 - e^(-Ke*tau)))
- Bayesian update: minimize sum of ((predicted - observed)^2 / variance)

Target AUC24/MIC: 400-600 mg*h/L (MRSA infections per ASHP/IDSA 2020 guidelines)
Trough targets: 15-20 mg/L (serious), 10-15 mg/L (less serious)

Author: Dr. Abu Suraih Sakhri
License: MIT
"""

import argparse
import csv
import json
import math
import sys
from typing import Dict, Any, List, Optional, Tuple


# ============================================================================
# Population PK Parameters
# ============================================================================

# Default population parameters
DEFAULT_VD_PER_KG = 0.7       # L/kg - population Vd
DEFAULT_KE_POPULATION = 0.1   # h^-1 - population elimination rate constant
DEFAULT_CL_PER_CRCL_SLOPE = 0.689  # Matzke equation slope
DEFAULT_CL_INTERCEPT = 3.66   # mL/min - Matzke equation intercept

# Target ranges
AUC_TARGET_MIN = 400.0   # mg*h/L
AUC_TARGET_MAX = 600.0   # mg*h/L
AUC_TARGET_MID = 500.0   # mg*h/L - midpoint for dosing calculations

TROUGH_SERIOUS_MIN = 15.0  # mg/L
TROUGH_SERIOUS_MAX = 20.0  # mg/L
TROUGH_MODERATE_MIN = 10.0 # mg/L
TROUGH_MODERATE_MAX = 15.0 # mg/L

MIC_DEFAULT = 1.0  # mg/L - typical MRSA MIC


# ============================================================================
# Core PK Functions
# ============================================================================

def calculate_crcl_cockcroft_gault(
    weight_kg: float,
    age_years: int,
    serum_creatinine_mg_dl: float,
    is_female: bool = False
) -> float:
    """
    Calculate creatinine clearance using Cockcroft-Gault equation.
    
    CrCl = ((140 - age) * weight) / (72 * SCr) * (0.85 if female)
    
    Args:
        weight_kg: Patient weight in kg
        age_years: Patient age in years
        serum_creatinine_mg_dl: Serum creatinine in mg/dL
        is_female: Whether patient is female
        
    Returns:
        CrCl in mL/min
    """
    if serum_creatinine_mg_dl <= 0:
        raise ValueError("Serum creatinine must be positive")
    if weight_kg <= 0:
        raise ValueError("Weight must be positive")
    if age_years < 0:
        raise ValueError("Age must be non-negative")
    
    crcl = ((140 - age_years) * weight_kg) / (72 * serum_creatinine_mg_dl)
    if is_female:
        crcl *= 0.85
    return round(crcl, 2)


def calculate_vancomycin_clearance(crcl_ml_min: float) -> float:
    """
    Calculate vancomycin clearance using Matzke equation.
    
    CL_vanco (mL/min) = CrCl * 0.689 + 3.66
    
    Args:
        crcl_ml_min: Creatinine clearance in mL/min
        
    Returns:
        Vancomycin clearance in mL/min
    """
    if crcl_ml_min < 0:
        raise ValueError("CrCl must be non-negative")
    return round(crcl_ml_min * DEFAULT_CL_PER_CRCL_SLOPE + DEFAULT_CL_INTERCEPT, 2)


def calculate_volume_of_distribution(
    weight_kg: float,
    vd_per_kg: float = DEFAULT_VD_PER_KG,
    obesity_adjusted: bool = False,
    actual_weight_kg: Optional[float] = None,
    ideal_body_weight_kg: Optional[float] = None
) -> float:
    """
    Calculate volume of distribution.
    
    For obese patients: Vd = Vd_per_kg * IBW + 0.4 * (ABW - IBW)
    For normal: Vd = Vd_per_kg * weight
    
    Args:
        weight_kg: Weight to use for calculation (kg)
        vd_per_kg: Volume of distribution per kg (L/kg)
        obesity_adjusted: Whether to use adjusted body weight
        actual_weight_kg: Actual body weight (for obese adjustment)
        ideal_body_weight_kg: Ideal body weight (for obese adjustment)
        
    Returns:
        Volume of distribution in liters
    """
    if weight_kg <= 0:
        raise ValueError("Weight must be positive")
    
    if obesity_adjusted and actual_weight_kg and ideal_body_weight_kg:
        if actual_weight_kg > ideal_body_weight_kg:
            vd = vd_per_kg * ideal_body_weight_kg + 0.4 * (actual_weight_kg - ideal_body_weight_kg)
        else:
            vd = vd_per_kg * weight_kg
    else:
        vd = vd_per_kg * weight_kg
    
    return round(vd, 2)


def calculate_elimination_rate_constant(cl_ml_min: float, vd_liters: float) -> float:
    """
    Calculate elimination rate constant.
    
    Ke = CL / Vd  (converting CL from mL/min to L/h)
    Ke = (CL_ml_min * 60 / 1000) / Vd_L
    
    Args:
        cl_ml_min: Clearance in mL/min
        vd_liters: Volume of distribution in liters
        
    Returns:
        Ke in h^-1
    """
    if vd_liters <= 0:
        raise ValueError("Vd must be positive")
    cl_l_per_hour = cl_ml_min * 60 / 1000
    ke = cl_l_per_hour / vd_liters
    return round(ke, 4)


def calculate_half_life(ke: float) -> float:
    """
    Calculate elimination half-life.
    
    t1/2 = 0.693 / Ke
    
    Args:
        ke: Elimination rate constant in h^-1
        
    Returns:
        Half-life in hours
    """
    if ke <= 0:
        raise ValueError("Ke must be positive")
    return round(0.693 / ke, 2)


def calculate_auc24_intermittent(
    dose_mg: float,
    cl_ml_min: float,
    tau_hours: float,
    infusion_hours: float = 1.0
) -> float:
    """
    Calculate AUC24 for intermittent vancomycin infusion.
    
    AUC24 = (Dose / CL_L_per_h) * (24 / tau)
    Where CL_L_per_h = CL_ml_min * 60 / 1000
    
    More precisely for intermittent infusion:
    AUC_interval = Dose / CL (in consistent units)
    AUC24 = AUC_interval * (24 / tau)
    
    Args:
        dose_mg: Vancomycin dose in mg
        cl_ml_min: Vancomycin clearance in mL/min
        tau_hours: Dosing interval in hours
        infusion_hours: Infusion duration in hours (default 1)
        
    Returns:
        AUC24 in mg*h/L
    """
    if cl_ml_min <= 0:
        raise ValueError("Clearance must be positive")
    if dose_mg <= 0:
        raise ValueError("Dose must be positive")
    if tau_hours <= 0:
        raise ValueError("Dosing interval must be positive")
    
    cl_l_per_h = cl_ml_min * 60 / 1000
    auc_interval = dose_mg / cl_l_per_h
    auc24 = auc_interval * (24 / tau_hours)
    return round(auc24, 1)


def calculate_trough_steady_state(
    dose_mg: float,
    vd_liters: float,
    ke: float,
    tau_hours: float,
    infusion_hours: float = 1.0
) -> float:
    """
    Calculate steady-state trough concentration for intermittent infusion.
    
    Ctrough = (Dose / (Vd * Ke * T_inf)) * (1 - e^(-Ke*T_inf)) * 
              e^(-Ke*(tau - T_inf)) / (1 - e^(-Ke*tau))
    
    Simplified (for 1-hour infusion approximation):
    Ctrough ≈ Cmax * e^(-Ke*(tau - T_inf))
    where Cmax = (Dose/Vd) / (1 - e^(-Ke*tau))
    
    Args:
        dose_mg: Vancomycin dose in mg
        vd_liters: Volume of distribution in liters
        ke: Elimination rate constant in h^-1
        tau_hours: Dosing interval in hours
        infusion_hours: Infusion duration in hours
        
    Returns:
        Trough concentration in mg/L
    """
    if vd_liters <= 0 or ke <= 0 or tau_hours <= 0:
        raise ValueError("Vd, Ke, and tau must be positive")
    
    # Steady-state Cmax (end of infusion)
    cmax = (dose_mg / vd_liters) * (1 - math.exp(-ke * infusion_hours)) / (1 - math.exp(-ke * tau_hours))
    
    # Trough = Cmax * e^(-ke * (tau - infusion))
    trough = cmax * math.exp(-ke * (tau_hours - infusion_hours))
    return round(trough, 2)


def calculate_peak_steady_state(
    dose_mg: float,
    vd_liters: float,
    ke: float,
    tau_hours: float,
    infusion_hours: float = 1.0
) -> float:
    """
    Calculate steady-state peak concentration (end of infusion).
    
    Cmax_ss = (Dose / (Vd * Ke * T_inf)) * (1 - e^(-Ke*T_inf)) / (1 - e^(-Ke*tau))
    
    Args:
        dose_mg: Vancomycin dose in mg
        vd_liters: Volume of distribution in liters
        ke: Elimination rate constant in h^-1
        tau_hours: Dosing interval in hours
        infusion_hours: Infusion duration in hours
        
    Returns:
        Peak concentration in mg/L
    """
    if vd_liters <= 0 or ke <= 0 or tau_hours <= 0:
        raise ValueError("Vd, Ke, and tau must be positive")
    
    cmax = (dose_mg / (vd_liters * ke * infusion_hours)) * (1 - math.exp(-ke * infusion_hours)) / (1 - math.exp(-ke * tau_hours))
    return round(cmax, 2)


def calculate_auc24_from_trough(
    trough_mg_l: float,
    ke: float,
    tau_hours: float
) -> float:
    """
    Estimate AUC24 from trough concentration using log-trapezoidal method.
    
    AUC_interval ≈ Ctrough / Ke (simplified steady-state approximation)
    AUC24 = AUC_interval * (24 / tau)
    
    More precisely:
    AUC_interval = Ctrough * (e^(Ke*tau) - 1) / Ke
    AUC24 = AUC_interval * (24 / tau)
    
    Args:
        trough_mg_l: Measured trough concentration in mg/L
        ke: Elimination rate constant in h^-1
        tau_hours: Dosing interval in hours
        
    Returns:
        Estimated AUC24 in mg*h/L
    """
    if ke <= 0:
        raise ValueError("Ke must be positive")
    if trough_mg_l < 0:
        raise ValueError("Trough must be non-negative")
    
    # AUC over one interval from trough
    auc_interval = trough_mg_l * (math.exp(ke * tau_hours) - 1) / ke
    auc24 = auc_interval * (24 / tau_hours)
    return round(auc24, 1)


def calculate_loading_dose(
    target_concentration_mg_l: float,
    vd_liters: float,
    bioavailability: float = 1.0
) -> float:
    """
    Calculate loading dose.
    
    LD = Vd * Ctarget / F
    
    Args:
        target_concentration_mg_l: Target peak concentration in mg/L
        vd_liters: Volume of distribution in liters
        bioavailability: Bioavailability (1.0 for IV)
        
    Returns:
        Loading dose in mg
    """
    if vd_liters <= 0:
        raise ValueError("Vd must be positive")
    if target_concentration_mg_l <= 0:
        raise ValueError("Target concentration must be positive")
    if bioavailability <= 0 or bioavailability > 1:
        raise ValueError("Bioavailability must be between 0 and 1")
    
    ld = vd_liters * target_concentration_mg_l / bioavailability
    return round(ld / 250) * 250  # Round to nearest 250mg


def adjust_dose_for_target_auc(
    current_dose_mg: float,
    measured_auc24: float,
    target_auc24: float = AUC_TARGET_MID
) -> float:
    """
    Adjust dose to achieve target AUC24.
    
    New Dose = Current Dose * (Target AUC / Measured AUC)
    
    Args:
        current_dose_mg: Current vancomycin dose in mg
        measured_auc24: Measured AUC24 in mg*h/L
        target_auc24: Target AUC24 in mg*h/L
        
    Returns:
        Adjusted dose in mg (rounded to nearest 250mg)
    """
    if measured_auc24 <= 0:
        raise ValueError("Measured AUC must be positive")
    if current_dose_mg <= 0:
        raise ValueError("Current dose must be positive")
    
    new_dose = current_dose_mg * (target_auc24 / measured_auc24)
    return round(new_dose / 250) * 250  # Round to nearest 250mg


def interpret_auc_mic(auc24: float, mic: float = MIC_DEFAULT) -> Dict[str, Any]:
    """
    Interpret AUC24/MIC ratio per ASHP/IDSA guidelines.
    
    Target: AUC24/MIC 400-600 for MRSA
    
    Args:
        auc24: Calculated AUC24 in mg*h/L
        mic: Minimum inhibitory concentration in mg/L
        
    Returns:
        Dictionary with interpretation
    """
    if mic <= 0:
        raise ValueError("MIC must be positive")
    
    auc_mic_ratio = auc24 / mic
    
    if auc_mic_ratio < 400:
        status = "SUBTHERAPEUTIC"
        recommendation = "Increase dose or decrease interval. AUC24/MIC below target 400-600."
        risk = "Treatment failure risk - inadequate exposure"
    elif auc_mic_ratio <= 600:
        status = "THERAPEUTIC"
        recommendation = "Current regimen achieves target AUC24/MIC of 400-600."
        risk = "Within therapeutic range"
    else:
        status = "SUPRATHERAPEUTIC"
        recommendation = "Decrease dose or increase interval. AUC24/MIC above target 600."
        risk = "Nephrotoxicity risk increases with AUC24 >600-700"
    
    return {
        "auc24": auc24,
        "mic": mic,
        "auc_mic_ratio": round(auc_mic_ratio, 1),
        "status": status,
        "recommendation": recommendation,
        "risk": risk,
        "target_range": "400-600"
    }


def interpret_trough(trough_mg_l: float, infection_severity: str = "serious") -> Dict[str, Any]:
    """
    Interpret vancomycin trough concentration.
    
    Args:
        trough_mg_l: Trough concentration in mg/L
        infection_severity: "serious" (15-20) or "moderate" (10-15)
        
    Returns:
        Dictionary with interpretation
    """
    if infection_severity == "serious":
        low, high = TROUGH_SERIOUS_MIN, TROUGH_SERIOUS_MAX
        target_desc = "15-20 mg/L (serious infections)"
    else:
        low, high = TROUGH_MODERATE_MIN, TROUGH_MODERATE_MAX
        target_desc = "10-15 mg/L (less serious infections)"
    
    if trough_mg_l < low:
        status = "SUBTHERAPEUTIC"
        recommendation = "Increase dose or decrease interval."
    elif trough_mg_l <= high:
        status = "THERAPEUTIC"
        recommendation = "Trough within target range. Continue current regimen."
    else:
        status = "SUPRATHERAPEUTIC"
        recommendation = "Decrease dose or increase interval. Monitor for nephrotoxicity."
    
    return {
        "trough": trough_mg_l,
        "target_range": target_desc,
        "status": status,
        "recommendation": recommendation
    }


def estimate_ke_from_two_levels(
    level1_mg_l: float,
    level2_mg_l: float,
    time_between_hours: float
) -> float:
    """
    Estimate Ke from two measured serum levels.
    
    Ke = ln(C1/C2) / (t2 - t1)
    
    Args:
        level1_mg_l: First level (higher, e.g., peak) in mg/L
        level2_mg_l: Second level (lower, e.g., trough) in mg/L
        time_between_hours: Time between measurements in hours
        
    Returns:
        Estimated Ke in h^-1
    """
    if level1_mg_l <= 0 or level2_mg_l <= 0:
        raise ValueError("Both levels must be positive")
    if level2_mg_l >= level1_mg_l:
        raise ValueError("First level must be higher than second (declining phase)")
    if time_between_hours <= 0:
        raise ValueError("Time between levels must be positive")
    
    ke = math.log(level1_mg_l / level2_mg_l) / time_between_hours
    return round(ke, 4)


def estimate_vd_from_level(
    dose_mg: float,
    ke: float,
    peak_mg_l: float,
    time_after_dose_hours: float,
    infusion_hours: float = 1.0
) -> float:
    """
    Estimate individual Vd from a measured level.
    
    C = (Dose / (Vd * Ke * T_inf)) * (1 - e^(-Ke*T_inf)) * e^(-Ke*(t - T_inf))
    Solve for Vd
    
    Args:
        dose_mg: Dose administered in mg
        ke: Elimination rate constant in h^-1
        peak_mg_l: Measured concentration in mg/L
        time_after_dose_hours: Time after start of infusion in hours
        infusion_hours: Infusion duration in hours
        
    Returns:
        Estimated Vd in liters
    """
    if ke <= 0 or peak_mg_l <= 0:
        raise ValueError("Ke and measured level must be positive")
    
    t_after_infusion = max(0, time_after_dose_hours - infusion_hours)
    
    # C = (Dose / (Vd * Ke * T_inf)) * (1 - e^(-Ke*T_inf)) * e^(-Ke*t_after)
    # Vd = (Dose / (C * Ke * T_inf)) * (1 - e^(-Ke*T_inf)) * e^(-Ke*t_after)
    
    vd = (dose_mg / (peak_mg_l * ke * infusion_hours)) * (1 - math.exp(-ke * infusion_hours)) * math.exp(-ke * t_after_infusion)
    return round(vd, 2)


def bayesian_estimate(
    measured_trough: float,
    measured_time_hours: float,
    dose_mg: float,
    tau_hours: float,
    weight_kg: float,
    crcl_ml_min: float,
    pop_vd: Optional[float] = None,
    pop_ke: Optional[float] = None,
    infusion_hours: float = 1.0
) -> Dict[str, Any]:
    """
    Bayesian parameter estimation using population prior + measured level.
    
    Uses a simple grid search over Vd and Ke to minimize the weighted
    sum of squared differences between predicted and observed trough,
    penalized by distance from population parameters.
    
    Objective = ((predicted_C - observed_C)^2 / sigma_C^2) + 
                ((Vd - Vd_pop)^2 / sigma_Vd^2) + 
                ((Ke - Ke_pop)^2 / sigma_Ke^2)
    
    Args:
        measured_trough: Observed trough concentration in mg/L
        measured_time_hours: Time of trough relative to dose (hours after dose start)
        dose_mg: Vancomycin dose in mg
        tau_hours: Dosing interval in hours
        weight_kg: Patient weight in kg
        crcl_ml_min: Creatinine clearance in mL/min
        pop_vd: Population Vd override (liters)
        pop_ke: Population Ke override (h^-1)
        infusion_hours: Infusion duration in hours
        
    Returns:
        Dictionary with Bayesian-estimated PK parameters
    """
    # Population priors
    vd_pop = pop_vd if pop_vd else calculate_volume_of_distribution(weight_kg)
    cl_pop = calculate_vancomycin_clearance(crcl_ml_min)
    ke_pop = pop_ke if pop_ke else calculate_elimination_rate_constant(cl_pop, vd_pop)
    
    # Variance for Bayesian weighting
    sigma_c = 3.0    # mg/L - expected measurement error
    sigma_vd = vd_pop * 0.3  # 30% CV for Vd
    sigma_ke = ke_pop * 0.3  # 30% CV for Ke
    
    best_cost = float('inf')
    best_vd = vd_pop
    best_ke = ke_pop
    
    # Grid search (simplified Bayesian - not full MCMC)
    for vd_factor in [0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4]:
        for ke_factor in [0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5]:
            vd_test = vd_pop * vd_factor
            ke_test = ke_pop * ke_factor
            
            if vd_test <= 0 or ke_test <= 0:
                continue
            
            # Predict trough at measured time
            try:
                predicted_trough = calculate_trough_steady_state(
                    dose_mg, vd_test, ke_test, tau_hours, infusion_hours
                )
            except (ValueError, ZeroDivisionError):
                continue
            
            # Bayesian cost function
            cost = ((predicted_trough - measured_trough) ** 2 / sigma_c ** 2 +
                    (vd_test - vd_pop) ** 2 / sigma_vd ** 2 +
                    (ke_test - ke_pop) ** 2 / sigma_ke ** 2)
            
            if cost < best_cost:
                best_cost = cost
                best_vd = vd_test
                best_ke = ke_test
    
    # Calculate derived parameters from Bayesian estimates
    best_cl = best_ke * best_vd * 1000 / 60  # Convert L/h to mL/min
    best_half_life = calculate_half_life(best_ke)
    auc24 = calculate_auc24_intermittent(dose_mg, best_cl, tau_hours, infusion_hours)
    predicted_trough = calculate_trough_steady_state(dose_mg, best_vd, best_ke, tau_hours, infusion_hours)
    
    return {
        "bayesian_vd": round(best_vd, 2),
        "bayesian_ke": round(best_ke, 4),
        "bayesian_cl": round(best_cl, 2),
        "bayesian_half_life": best_half_life,
        "population_vd": round(vd_pop, 2),
        "population_ke": round(ke_pop, 4),
        "population_cl": round(cl_pop, 2),
        "predicted_trough": predicted_trough,
        "measured_trough": measured_trough,
        "auc24": auc24,
        "fit_cost": round(best_cost, 4)
    }


def full_dosing_assessment(
    weight_kg: float,
    age_years: int,
    serum_creatinine_mg_dl: float,
    is_female: bool,
    current_dose_mg: float,
    tau_hours: float,
    measured_trough: Optional[float] = None,
    mic: float = MIC_DEFAULT,
    infection_severity: str = "serious",
    infusion_hours: float = 1.0
) -> Dict[str, Any]:
    """
    Complete vancomycin dosing assessment.
    
    Args:
        weight_kg: Patient weight in kg
        age_years: Patient age in years
        serum_creatinine_mg_dl: Serum creatinine in mg/dL
        is_female: Whether patient is female
        current_dose_mg: Current vancomycin dose in mg
        tau_hours: Dosing interval in hours
        measured_trough: Optional measured trough in mg/L
        mic: MIC of organism in mg/L
        infection_severity: "serious" or "moderate"
        infusion_hours: Infusion duration in hours
        
    Returns:
        Complete dosing assessment dictionary
    """
    # Step 1: Calculate CrCl
    crcl = calculate_crcl_cockcroft_gault(weight_kg, age_years, serum_creatinine_mg_dl, is_female)
    
    # Step 2: Calculate vancomycin clearance
    cl_vanco = calculate_vancomycin_clearance(crcl)
    
    # Step 3: Calculate Vd
    vd = calculate_volume_of_distribution(weight_kg)
    
    # Step 4: Calculate Ke and half-life
    ke = calculate_elimination_rate_constant(cl_vanco, vd)
    half_life = calculate_half_life(ke)
    
    # Step 5: Calculate predicted steady-state concentrations
    predicted_trough = calculate_trough_steady_state(current_dose_mg, vd, ke, tau_hours, infusion_hours)
    predicted_peak = calculate_peak_steady_state(current_dose_mg, vd, ke, tau_hours, infusion_hours)
    
    # Step 6: Calculate AUC24
    auc24 = calculate_auc24_intermittent(current_dose_mg, cl_vanco, tau_hours)
    
    # Step 7: Interpret AUC24/MIC
    auc_mic_interp = interpret_auc_mic(auc24, mic)
    
    # Step 8: If measured trough available, do Bayesian estimation
    bayesian = None
    if measured_trough is not None:
        bayesian = bayesian_estimate(
            measured_trough, tau_hours, current_dose_mg, tau_hours,
            weight_kg, crcl, infusion_hours=infusion_hours
        )
        # Use Bayesian AUC24 if available
        auc24 = bayesian["auc24"]
        auc_mic_interp = interpret_auc_mic(auc24, mic)
    
    # Step 9: Dose adjustment recommendation
    dose_adjustment = None
    if measured_trough is not None:
        # Use measured trough for dose adjustment
        trough_interp = interpret_trough(measured_trough, infection_severity)
        if trough_interp["status"] == "SUBTHERAPEUTIC":
            target_trough = TROUGH_SERIOUS_MIN if infection_severity == "serious" else TROUGH_MODERATE_MIN
            dose_ratio = target_trough / measured_trough if measured_trough > 0 else 1.5
            new_dose = round(current_dose_mg * dose_ratio / 250) * 250
            dose_adjustment = {
                "action": "INCREASE",
                "current_dose": current_dose_mg,
                "recommended_dose": new_dose,
                "interval_hours": tau_hours,
                "reason": f"Trough {measured_trough} mg/L is below target"
            }
        elif trough_interp["status"] == "SUPRATHERAPEUTIC":
            target_trough = TROUGH_SERIOUS_MAX if infection_severity == "serious" else TROUGH_MODERATE_MAX
            dose_ratio = target_trough / measured_trough if measured_trough > 0 else 0.7
            new_dose = round(current_dose_mg * dose_ratio / 250) * 250
            dose_adjustment = {
                "action": "DECREASE",
                "current_dose": current_dose_mg,
                "recommended_dose": new_dose,
                "interval_hours": tau_hours,
                "reason": f"Trough {measured_trough} mg/L is above target"
            }
        else:
            dose_adjustment = {
                "action": "NO_CHANGE",
                "current_dose": current_dose_mg,
                "recommended_dose": current_dose_mg,
                "interval_hours": tau_hours,
                "reason": "Trough within therapeutic range"
            }
    else:
        # Use AUC-based adjustment
        if auc_mic_interp["status"] == "SUBTHERAPEUTIC":
            new_dose = adjust_dose_for_target_auc(current_dose_mg, auc24, AUC_TARGET_MID)
            dose_adjustment = {
                "action": "INCREASE",
                "current_dose": current_dose_mg,
                "recommended_dose": new_dose,
                "interval_hours": tau_hours,
                "reason": f"AUC24/MIC {auc_mic_interp['auc_mic_ratio']} below target 400-600"
            }
        elif auc_mic_interp["status"] == "SUPRATHERAPEUTIC":
            new_dose = adjust_dose_for_target_auc(current_dose_mg, auc24, AUC_TARGET_MID)
            dose_adjustment = {
                "action": "DECREASE",
                "current_dose": current_dose_mg,
                "recommended_dose": new_dose,
                "interval_hours": tau_hours,
                "reason": f"AUC24/MIC {auc_mic_interp['auc_mic_ratio']} above target 400-600"
            }
        else:
            dose_adjustment = {
                "action": "NO_CHANGE",
                "current_dose": current_dose_mg,
                "recommended_dose": current_dose_mg,
                "interval_hours": tau_hours,
                "reason": "AUC24/MIC within target range"
            }
    
    # Step 10: Loading dose recommendation
    loading_dose = None
    if measured_trough is not None and measured_trough < 10:
        target_peak = 25.0  # mg/L
        loading_dose = calculate_loading_dose(target_peak, vd)
    
    return {
        "patient_parameters": {
            "weight_kg": weight_kg,
            "age_years": age_years,
            "serum_creatinine_mg_dl": serum_creatinine_mg_dl,
            "is_female": is_female,
            "crcl_ml_min": crcl
        },
        "pk_parameters": {
            "vd_liters": vd,
            "cl_vanco_ml_min": cl_vanco,
            "ke_per_hour": ke,
            "half_life_hours": half_life
        },
        "current_regimen": {
            "dose_mg": current_dose_mg,
            "interval_hours": tau_hours,
            "infusion_hours": infusion_hours
        },
        "predicted_levels": {
            "trough_mg_l": predicted_trough,
            "peak_mg_l": predicted_peak
        },
        "auc24_analysis": auc_mic_interp,
        "bayesian_estimation": bayesian,
        "dose_recommendation": dose_adjustment,
        "loading_dose_mg": loading_dose,
        "disclaimer": "FOR EDUCATIONAL/RESEARCH USE ONLY. Not a substitute for clinical pharmacist review."
    }


def main(argv=None):
    """CLI entry point for vancomycin dosing calculator."""
    parser = argparse.ArgumentParser(
        prog="vanco-bayes",
        description="Vancomycin AUC24/MIC Bayesian Dosing Calculator"
    )
    subparsers = parser.add_subparsers(dest="command", required=True)
    
    # --- Assess command ---
    assess_parser = subparsers.add_parser("assess", help="Full dosing assessment")
    assess_parser.add_argument("--weight", type=float, required=True, help="Weight in kg")
    assess_parser.add_argument("--age", type=int, required=True, help="Age in years")
    assess_parser.add_argument("--scr", type=float, required=True, help="Serum creatinine mg/dL")
    assess_parser.add_argument("--female", action="store_true", help="Patient is female")
    assess_parser.add_argument("--dose", type=float, required=True, help="Current dose in mg")
    assess_parser.add_argument("--interval", type=float, required=True, help="Dosing interval in hours")
    assess_parser.add_argument("--trough", type=float, help="Measured trough mg/L")
    assess_parser.add_argument("--mic", type=float, default=MIC_DEFAULT, help="MIC mg/L (default 1.0)")
    assess_parser.add_argument("--severity", choices=["serious", "moderate"], default="serious",
                               help="Infection severity (default: serious)")
    
    # --- AUC command ---
    auc_parser = subparsers.add_parser("auc", help="Calculate AUC24 from dose and clearance")
    auc_parser.add_argument("--dose", type=float, required=True, help="Dose in mg")
    auc_parser.add_argument("--cl", type=float, required=True, help="Clearance in mL/min")
    auc_parser.add_argument("--interval", type=float, required=True, help="Dosing interval in hours")
    
    # --- CrCl command ---
    crcl_parser = subparsers.add_parser("crcl", help="Calculate creatinine clearance")
    crcl_parser.add_argument("--weight", type=float, required=True, help="Weight in kg")
    crcl_parser.add_argument("--age", type=int, required=True, help="Age in years")
    crcl_parser.add_argument("--scr", type=float, required=True, help="Serum creatinine mg/dL")
    crcl_parser.add_argument("--female", action="store_true", help="Patient is female")
    
    # --- Loading dose command ---
    ld_parser = subparsers.add_parser("loading-dose", help="Calculate loading dose")
    ld_parser.add_argument("--weight", type=float, required=True, help="Weight in kg")
    ld_parser.add_argument("--target", type=float, default=25.0, help="Target peak mg/L (default 25)")
    
    # --- Bayesian command ---
    bayes_parser = subparsers.add_parser("bayesian", help="Bayesian parameter estimation")
    bayes_parser.add_argument("--trough", type=float, required=True, help="Measured trough mg/L")
    bayes_parser.add_argument("--dose", type=float, required=True, help="Dose in mg")
    bayes_parser.add_argument("--interval", type=float, required=True, help="Dosing interval hours")
    bayes_parser.add_argument("--weight", type=float, required=True, help="Weight in kg")
    bayes_parser.add_argument("--crcl", type=float, required=True, help="CrCl mL/min")
    
    args = parser.parse_args(argv)
    
    if args.command == "assess":
        result = full_dosing_assessment(
            weight_kg=args.weight,
            age_years=args.age,
            serum_creatinine_mg_dl=args.scr,
            is_female=args.female,
            current_dose_mg=args.dose,
            tau_hours=args.interval,
            measured_trough=args.trough,
            mic=args.mic,
            infection_severity=args.severity
        )
        print(json.dumps(result, indent=2))
    
    elif args.command == "auc":
        auc24 = calculate_auc24_intermittent(args.dose, args.cl, args.interval)
        interp = interpret_auc_mic(auc24)
        print(json.dumps({"auc24": auc24, "interpretation": interp}, indent=2))
    
    elif args.command == "crcl":
        crcl = calculate_crcl_cockcroft_gault(args.weight, args.age, args.scr, args.female)
        cl_vanco = calculate_vancomycin_clearance(crcl)
        print(json.dumps({"crcl_ml_min": crcl, "vancomycin_cl_ml_min": cl_vanco}, indent=2))
    
    elif args.command == "loading-dose":
        vd = calculate_volume_of_distribution(args.weight)
        ld = calculate_loading_dose(args.target, vd)
        print(json.dumps({"vd_liters": vd, "loading_dose_mg": ld}, indent=2))
    
    elif args.command == "bayesian":
        result = bayesian_estimate(
            measured_trough=args.trough,
            measured_time_hours=args.interval,
            dose_mg=args.dose,
            tau_hours=args.interval,
            weight_kg=args.weight,
            crcl_ml_min=args.crcl
        )
        print(json.dumps(result, indent=2))
    
    return 0


if __name__ == "__main__":
    sys.exit(main())
