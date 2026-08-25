# Vancomycin AUC24/MIC Bayesian Dosing Calculator

Real pharmacokinetic calculator for vancomycin therapeutic drug monitoring using population PK with Bayesian individualization.

## Clinical Background

Vancomycin dosing is guided by the **AUC24/MIC ratio** per ASHP/IDSA 2020 consensus guidelines:
- **Target AUC24/MIC: 400-600 mg·h/L** for MRSA infections
- Trough-only monitoring is no longer recommended as primary target
- Bayesian estimation provides more accurate individual PK parameters

## Key Formulas

| Parameter | Formula |
|-----------|---------|
| CrCl (Cockcroft-Gault) | `((140 - age) × weight) / (72 × SCr) × 0.85 if female` |
| Vancomycin CL (Matzke) | `CL = CrCl × 0.689 + 3.66` (mL/min) |
| Volume of distribution | `Vd = 0.7 L/kg` (population) |
| Elimination rate constant | `Ke = CL / Vd` |
| Half-life | `t1/2 = 0.693 / Ke` |
| AUC24 (intermittent) | `AUC24 = (Dose / CL) × (24 / τ)` |
| Loading dose | `LD = Vd × Ctarget / F` |
| Dose adjustment | `New Dose = Current × (Target AUC / Measured AUC)` |

## Installation

```bash
# No dependencies required - Python 3.8+ stdlib only
cd vancomycin-bayesian-auc-agent
```

## Usage

### Full Dosing Assessment
```bash
python cli.py assess --weight 80 --age 65 --scr 1.2 --dose 1000 --interval 12
python cli.py assess --weight 80 --age 65 --scr 1.2 --dose 1000 --interval 12 --trough 18.5
python cli.py assess --weight 80 --age 65 --scr 1.2 --dose 1000 --interval 12 --trough 8.0 --mic 1.0 --severity serious
```

### Calculate AUC24
```bash
python cli.py auc --dose 1000 --cl 72.56 --interval 12
```

### Calculate Creatinine Clearance
```bash
python cli.py crcl --weight 80 --age 65 --scr 1.2
python cli.py crcl --weight 60 --age 70 --scr 1.5 --female
```

### Calculate Loading Dose
```bash
python cli.py loading-dose --weight 80 --target 25
```

### Bayesian Parameter Estimation
```bash
python cli.py bayesian --trough 15.0 --dose 1000 --interval 12 --weight 80 --crcl 72.56
```

## Output Format

All commands output JSON. Example assessment output:
```json
{
  "patient_parameters": {
    "weight_kg": 80.0,
    "age_years": 65,
    "serum_creatinine_mg_dl": 1.2,
    "crcl_ml_min": 69.44
  },
  "pk_parameters": {
    "vd_liters": 56.0,
    "cl_vanco_ml_min": 51.57,
    "ke_per_hour": 0.0154,
    "half_life_hours": 45.0
  },
  "auc24_analysis": {
    "auc24": 518.4,
    "auc_mic_ratio": 518.4,
    "status": "THERAPEUTIC"
  }
}
```

## Tests

```bash
python -m pytest test_vanco_bayes.py -v
```

## Disclaimer

**FOR EDUCATIONAL AND RESEARCH USE ONLY.** This calculator is not a substitute for clinical pharmacist review and institutional dosing protocols. All dosing decisions should be made by qualified healthcare professionals.

## References

- Rybak MJ, et al. Therapeutic monitoring of vancomycin for serious MRSA infections. *Am J Health Syst Pharm*. 2020;77(11):835-864.
- Matzke GR, et al. Vancomycin pharmacokinetics. *Clin Pharm*. 1984;3:216-222.
- Cockcroft DW, Gault MH. Prediction of creatinine clearance from serum creatinine. *Nephron*. 1976;16(1):31-41.

## License

MIT License. See [LICENSE](LICENSE) for details.
