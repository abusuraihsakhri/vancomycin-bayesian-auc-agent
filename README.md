# Vancomycin Bayesian AUC Agent

Research-only single-compartment vancomycin PK explorer. **Not clinically validated: never use to select a patient dose. Do not enter patient identifiers or PHI.**

## Browser application

The static, dependency-free application in `web/` runs entirely in the browser, with optional one-level penalized parameter fitting, an AUC exposure comparison, and local JSON export. No personal data is uploaded or stored.

```bash
python -m http.server 8001 --directory web
# Open http://localhost:8001/
```

GitHub Pages uses `.github/workflows/pages.yml` when Pages is configured for GitHub Actions. A public live URL must be verified before publication here.

## Python calculator

Python 3.10+:

```bash
python -m pip install -e '.[api,test]'
vanco-bayes assess --weight 80 --age 65 --scr 1.2 --dose 1000 --interval 12
vanco-bayes auc --dose 1000 --cl 72.56 --interval 12
vanco-bayes bayesian --trough 15 --dose 1000 --interval 12 --weight 80 --crcl 70
```

The assessment returns `CLINICIAN_REVIEW` and **no automated replacement dose**. Standalone legacy loading-dose and trough helpers remain for backward compatibility but are not clinically actionable.

## Clinical model limits

The implementation uses Cockcroft–Gault CrCl, an empirical clearance relation, 0.7 L/kg population distribution volume, constant-rate IV infusion, steady state, and a coarse weighted grid fit. This is not a validated clinical Bayesian model; it lacks credible intervals and support for complex dosing histories, unstable renal function, dialysis, obesity, or pediatric dosing.

The ASHP/IDSA/PIDS/SIDP 2020 consensus AUC24 target of 400–600 mg·h/L concerns serious MRSA infection **assuming a broth-microdilution MIC of 1 mg/L**. Trough-only targeting of 15–20 mg/L is not recommended for serious MRSA. [Consensus guideline](https://www.idsociety.org/practice-guideline/vancomycin/).

## Legacy agent services

`agents/` and `vancomycin_bayesian_auc_agent/` provide separate generic threshold-audit interfaces; these **do not implement Bayesian vancomycin AUC** and should not be treated as clinical guideline validation.

```bash
vancomycin-bayesian-auc-agent audit
vancomycin-bayesian-auc-agent serve --host 127.0.0.1 --port 8000
docker compose up --build
```

The Docker service exposes the unauthenticated generic audit API, not the browser calculator. Do not expose it publicly or send it protected health information. Its regex identifier guard is not comprehensive de-identification.

## Tests

```bash
python -m pytest -q
node --test web/calculator.test.cjs
python -m compileall -q agents vancomycin_bayesian_auc_agent vanco_bayes.py enrichment.py
```

CI tests Python 3.10–3.12, Node browser-model checks and the package build. Browser data remain local; manual JSON export is the only file output. This repository is [MIT licensed](LICENSE).
