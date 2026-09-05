# Vancomycin Bayesian AUC Agent

> **Domain:** Clinical Pharmacology & Precision Pharmacotherapy  
> **Reference Guidelines & Standards:** `CPIC Guidelines & FDA Table of Pharmacogenomic Biomarkers`

<div align="center">

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
![Python](https://img.shields.io/badge/Python-3.10%20%7C%203.11%20%7C%203.12-3776AB.svg?logo=python&logoColor=white)
![FastAPI](https://img.shields.io/badge/FastAPI-0.111-009688.svg?logo=fastapi&logoColor=white)
![Audit Trail](https://img.shields.io/badge/Audit-HMAC--SHA256_Tamper--Evident-brightgreen.svg)
![Zero-PHI Guard](https://img.shields.io/badge/Guard-Zero--PHI_Outbound-blue.svg)
![Docker](https://img.shields.io/badge/Docker-Ready-2496ED.svg?logo=docker&logoColor=white)

</div>

---

## 📖 What It Does

**Vancomycin Bayesian AUC Agent** is an advanced analytical and computational platform implementing Bayesian 24h AUC24/MIC (400-600) Precision Pharmacokinetic Engine.

---

## ⚙️ Key Capabilities & Algorithmic Modules

### 🔬 Analytical Functions

- **`calculate_crcl_cockcroft_gault()`**: Calculate creatinine clearance using Cockcroft-Gault equation.

CrCl = ((140 - age) * weight) / (72 * SCr) * (0.85 if female)

Args:
    weight_kg: Patient weight in kg
    age_years: Patient age in years
    serum_creatinine_mg_dl: Serum creatinine in mg/dL
    is_female: Whether patient is female
    
Returns:
    CrCl in mL/min
- **`calculate_vancomycin_clearance()`**: Calculate vancomycin clearance using Matzke equation.

CL_vanco (mL/min) = CrCl * 0.689 + 3.66

Args:
    crcl_ml_min: Creatinine clearance in mL/min
    
Returns:
    Vancomycin clearance in mL/min
- **`calculate_volume_of_distribution()`**: Calculate volume of distribution.

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
- **`calculate_elimination_rate_constant()`**: Calculate elimination rate constant.

Ke = CL / Vd  (converting CL from mL/min to L/h)
Ke = (CL_ml_min * 60 / 1000) / Vd_L

Args:
    cl_ml_min: Clearance in mL/min
    vd_liters: Volume of distribution in liters
    
Returns:
    Ke in h^-1
- **`calculate_half_life()`**: Calculate elimination half-life.

t1/2 = 0.693 / Ke

Args:
    ke: Elimination rate constant in h^-1
    
Returns:
    Half-life in hours

---

## 📐 Mathematical Formulation & Logic

```text
  Key formulas:
  Calculate creatinine clearance using Cockcroft-Gault equation.
  Calculate vancomycin clearance using Matzke equation.
  Calculate volume of distribution.
  Calculate elimination rate constant.
```

---

## 💻 CLI Quickstart & Usage

### 1. Guided Interactive Mode
```bash
python cli.py
```

### 2. Direct Parameterized Evaluation
```bash
python cli.py --input data.csv
```

### Parameter Reference
- `--interactive`: Launch guided terminal interactive wizard.
- `--input <path>`: Evaluate input from JSON or CSV specification.
- `--json`: Output deterministic structured results in JSON format.

### Input Data Schema

| Field | Description | Requirement |
|:------|:------------|:------------|
| `case_id` | Parameter / observation metric | Required |
| `patient_synthetic_id` | Parameter / observation metric | Required |
| `metric_primary` | Parameter / observation metric | Required |
| `metric_secondary` | Parameter / observation metric | Required |
| `is_stat` | Parameter / observation metric | Required |
| `status_flag` | Parameter / observation metric | Required |

---

## 🛡️ Security & Enterprise Architecture

* **Zero-PHI Outbound Interceptor:** Active AST and regex inspection blocking SSNs, MRNs, phone numbers, and patient identifiers.
* **Tamper-Evident HMAC-SHA256 Audit Trail:** Chained, cryptographically signed logs for every evaluation and state transition.
* **Air-Gapped LLM Reasoning Adapter:** Agnostic integration for local Ollama instances (`llama3`, `mistral`), Claude 3.5 Sonnet, GPT-4o, and deterministic test mocks.
* **Active Learning Bayesian Calibration:** Dynamic tracker updating worker reliability weights and monitoring Brier calibration drift.
* **FastAPI & Prometheus Telemetry:** Exposes OpenAPI 3.1 REST endpoints and operational Prometheus metrics (`/metrics`).

---

## 🧪 Testing & Verification

Run the automated test suite:

```bash
pytest -v
```

Execute high-throughput batch simulation benchmarks:

```bash
python simulator.py --tasks 1000
```

---

## 🐳 Container Deployment

```bash
docker build -t vancomycin-bayesian-auc-agent .
docker run -p 8000:8000 vancomycin-bayesian-auc-agent
```

For production deployments, set a stable audit secret key:

```bash
docker run -e AUDIT_SECRET_KEY=$(openssl rand -hex 32) -p 8000:8000 vancomycin-bayesian-auc-agent
```
