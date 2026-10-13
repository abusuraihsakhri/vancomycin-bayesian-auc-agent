"use strict";
const test = require("node:test");
const assert = require("node:assert/strict");
const fs = require("node:fs");
const path = require("node:path");
const {calculate, concentration} = require("./calculator.js");

const scenario = {
  weight: 80, age: 65, scr: 1.2, female: false,
  dose: 1000, interval: 12, infusion: 1, mic: 1,
  level: "", sampleTime: ""
};

test("population AUC follows total daily dose / clearance", () => {
  const result = calculate(scenario);
  const expected = scenario.dose / (result.clearance_ml_min * 0.06) * 24 / scenario.interval;
  assert.ok(Math.abs(result.auc24_mg_h_l - expected) < 1e-8);
  assert.ok(result.steady_state_peak_mg_l > result.steady_state_trough_mg_l);
  assert.equal(result.fit, null);
});

test("steady-state concentration is continuous across interval boundary", () => {
  const dose = 1000, vd = 50, ke = 0.1, tau = 12, infusion = 1;
  const preDose = concentration(dose, vd, ke, tau, infusion, tau);
  const postPriorDose = concentration(dose, vd, ke, tau, infusion, 0);
  assert.ok(Math.abs(preDose - postPriorDose) < 1e-10);
  assert.ok(concentration(dose, vd, ke, tau, infusion, infusion) > preDose);
});

test("Bayesian toy fit uses actual sampling time", () => {
  const a = calculate({...scenario, level: 15, sampleTime: 2});
  const b = calculate({...scenario, level: 15, sampleTime: 12});
  assert.notDeepEqual(
    [a.volume_l, a.elimination_constant_per_h],
    [b.volume_l, b.elimination_constant_per_h]
  );
  assert.equal(a.fit.time_after_infusion_start_h, 2);
});

test("no numeric input may be missing, nonfinite or outside regimen", () => {
  for (const value of ["", null, "Infinity", "NaN", -1, 0]) {
    assert.throws(() => calculate({...scenario, dose: value}), /positive finite/);
  }
  assert.throws(() => calculate({...scenario, infusion: 13}), /Infusion/);
  assert.throws(() => calculate({...scenario, level: 15, sampleTime: ""}), /Sample time/);
  assert.throws(() => calculate({...scenario, level: 15, sampleTime: 13}), /Sample time/);
  assert.throws(() => calculate({...scenario, age: 12}), /18/);
});

test("reference comparison is restricted to assumed MIC 1", () => {
  const result = calculate({...scenario, mic: 2});
  assert.match(result.guideline_comparison, /MIC differs/);
});

test("static UI does not fabricate audit outcomes or fetch API responses", () => {
  const html = fs.readFileSync(path.join(__dirname, "index.html"), "utf8");
  const ui = fs.readFileSync(path.join(__dirname, "ui.js"), "utf8");
  assert.doesNotMatch(html + ui, /fetch\s*\(/);
  assert.doesNotMatch(html + ui, /HMAC_SHA256_e83f|DOSSIER-SIM/);
  assert.match(html, /<label\s+for="weight">/);
  assert.match(html, /research and education only/i);
});
