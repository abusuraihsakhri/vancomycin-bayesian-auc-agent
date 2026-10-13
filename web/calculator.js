"use strict";
// Research-only, single-compartment vancomycin PK simulation.
// No patient data are transmitted or persisted by this module.
(function (root, factory) {
  const api = factory();
  if (typeof module === "object" && module.exports) module.exports = api;
  else root.VancoPK = api;
})(typeof globalThis !== "undefined" ? globalThis : this, function () {
  function finitePositive(value, name) {
    const n = Number(value);
    if (value === "" || value === null || value === undefined || !Number.isFinite(n) || n <= 0) {
      throw new Error(name + " must be a positive finite number.");
    }
    return n;
  }

  function concentration(dose, volume, ke, interval, infusion, time) {
    if (!(infusion > 0 && infusion <= interval && time >= 0 && time <= interval)) {
      throw new Error("Infusion and sample time must be within the dosing interval.");
    }
    const cmax = dose * (-Math.expm1(-ke * infusion)) /
      (volume * ke * infusion * (-Math.expm1(-ke * interval)));
    const cmin = cmax * Math.exp(-ke * (interval - infusion));
    if (time <= infusion) {
      return cmin * Math.exp(-ke * time) +
        dose * (-Math.expm1(-ke * time)) / (volume * ke * infusion);
    }
    return cmax * Math.exp(-ke * (time - infusion));
  }

  function calculate(input) {
    const weight = finitePositive(input.weight, "Weight");
    const age = finitePositive(input.age, "Age");
    const scr = finitePositive(input.scr, "Serum creatinine");
    const dose = finitePositive(input.dose, "Dose");
    const interval = finitePositive(input.interval, "Dosing interval");
    const infusion = finitePositive(input.infusion, "Infusion duration");
    const mic = finitePositive(input.mic, "MIC");
    if (age < 18 || age > 120) throw new Error("This illustrative model supports ages 18–120 only.");
    if (infusion > interval) throw new Error("Infusion cannot exceed the interval.");
    const crcl = (140 - age) * weight / (72 * scr) * (input.female ? 0.85 : 1);
    const populationCL = crcl * 0.689 + 3.66; // mL/min
    const populationVd = 0.7 * weight; // L
    const populationKe = populationCL * 0.06 / populationVd; // h^-1
    let volume = populationVd;
    let ke = populationKe;
    let fit = null;

    const hasLevel = input.level !== "" && input.level !== undefined &&
      input.level !== null;
    if (hasLevel) {
      const level = finitePositive(input.level, "Measured concentration");
      if (input.sampleTime === "" || input.sampleTime === undefined ||
          input.sampleTime === null || !Number.isFinite(Number(input.sampleTime))) {
        throw new Error("Sample time is required when a measured concentration is provided.");
      }
      const time = Number(input.sampleTime);
      if (time < 0 || time > interval) throw new Error("Sample time must lie within the interval.");
      const sigmaLevel = 3.0;
      const sigmaVd = populationVd * 0.3;
      const sigmaKe = populationKe * 0.3;
      let bestCost = Infinity;
      for (const vf of [0.6, 0.7, 0.8, 0.9, 1, 1.1, 1.2, 1.3, 1.4]) {
        for (const kf of [0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.1, 1.2, 1.3, 1.4, 1.5]) {
          const vd = populationVd * vf;
          const candidateKe = populationKe * kf;
          const predicted = concentration(dose, vd, candidateKe, interval, infusion, time);
          const cost = Math.pow((predicted - level) / sigmaLevel, 2) +
            Math.pow((vd - populationVd) / sigmaVd, 2) +
            Math.pow((candidateKe - populationKe) / sigmaKe, 2);
          if (cost < bestCost) {
            volume = vd;
            ke = candidateKe;
            bestCost = cost;
          }
        }
      }
      fit = {
        observed_mg_l: level,
        time_after_infusion_start_h: time,
        predicted_mg_l: concentration(dose, volume, ke, interval, infusion, time),
        weighted_fit_cost: bestCost
      };
    }
    const clearance = volume * ke / 0.06; // mL/min
    const auc = (dose / (volume * ke)) * (24 / interval);
    const peak = concentration(dose, volume, ke, interval, infusion, infusion);
    const trough = concentration(dose, volume, ke, interval, infusion, interval);
    const targetApplicable = Math.abs(mic - 1) < 1e-12;
    return {
      mode: hasLevel ? "illustrative penalized grid fit (steady state)" : "population estimate",
      creatinine_clearance_ml_min: crcl,
      clearance_ml_min: clearance,
      volume_l: volume,
      elimination_constant_per_h: ke,
      half_life_h: Math.log(2) / ke,
      auc24_mg_h_l: auc,
      auc24_mic: auc / mic,
      steady_state_peak_mg_l: peak,
      steady_state_trough_mg_l: trough,
      guideline_comparison: targetApplicable ?
        (auc < 400 ? "Below 400 mg·h/L" : auc > 600 ? "Above 600 mg·h/L" : "Between 400 and 600 mg·h/L") :
        "No direct 400–600 interpretation: MIC differs from assumed 1 mg/L",
      fit: fit,
      disclaimer: "Illustrative steady-state one-compartment model only; not validated for patient care or dosing decisions."
    };
  }

  return { calculate, concentration };
});
