"use strict";
(function () {
  const form = document.getElementById("pk-form");
  const errors = document.getElementById("errors");
  const results = document.getElementById("results");
  const empty = document.getElementById("empty");
  const exportButton = document.getElementById("export");
  let report = null;
  const field = (id) => document.getElementById(id).value;
  const setText = (id, text) => { document.getElementById(id).textContent = text; };
  const fmt = (value, digits) => Number(value).toFixed(digits);

  form.addEventListener("submit", function (event) {
    event.preventDefault();
    errors.textContent = "";
    results.hidden = true;
    empty.hidden = false;
    exportButton.disabled = true;
    report = null;
    try {
      const input = {
        weight: field("weight"), age: field("age"), scr: field("scr"),
        female: field("female") === "true", dose: field("dose"),
        interval: field("interval"), infusion: field("infusion"), mic: field("mic"),
        level: field("level"), sampleTime: field("sampleTime")
      };
      report = window.VancoPK.calculate(input);
      setText("mode", "Mode: " + report.mode);
      setText("auc", fmt(report.auc24_mg_h_l, 1));
      setText("ratio", fmt(report.auc24_mic, 1));
      setText("cl", fmt(report.clearance_ml_min * 0.06, 2));
      setText("vd", fmt(report.volume_l, 1));
      setText("peak", fmt(report.steady_state_peak_mg_l, 2));
      setText("trough", fmt(report.steady_state_trough_mg_l, 2));
      setText("half", fmt(report.half_life_h, 1));
      setText("crcl", fmt(report.creatinine_clearance_ml_min, 1));
      setText("interpretation", report.guideline_comparison +
        ". This is a mathematical comparison, not a treatment recommendation.");
      setText("fit", report.fit ?
        ("Observed level: " + fmt(report.fit.observed_mg_l, 2) + " mg/L at " +
          fmt(report.fit.time_after_infusion_start_h, 1) +
          " hours; model-predicted: " + fmt(report.fit.predicted_mg_l, 2) +
          " mg/L. Fit is not clinically validated.") :
        "No concentration supplied; all estimates use population assumptions.");
      empty.hidden = true;
      results.hidden = false;
      exportButton.disabled = false;
    } catch (error) {
      errors.textContent = error instanceof Error ? error.message : "Calculation failed.";
    }
  });

  exportButton.addEventListener("click", function () {
    if (!report) return;
    const blob = new Blob([JSON.stringify(report, null, 2)], {type: "application/json"});
    const url = URL.createObjectURL(blob);
    const anchor = document.createElement("a");
    anchor.href = url;
    anchor.download = "vancomycin-pk-research-report.json";
    document.body.appendChild(anchor);
    anchor.click();
    anchor.remove();
    window.setTimeout(function () { URL.revokeObjectURL(url); }, 0);
  });
})();
