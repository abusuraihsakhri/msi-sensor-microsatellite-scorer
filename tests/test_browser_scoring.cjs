"use strict";
const test = require("node:test");
const assert = require("node:assert/strict");
const {scoreMarkers, interpretIHC, combined, parseCSV, batchCSV} = require("../web/scoring.js");
const markers = {"BAT-25":true,"BAT-26":true,"D2S123":false,"D5S346":false,"D17S250":false};

test("Bethesda threshold: MSI-H, MSI-L and MSS", () => {
  assert.equal(scoreMarkers(markers).msi_status, "MSI-H");
  assert.equal(scoreMarkers({...markers,"BAT-26":false}).msi_status, "MSI-L");
  assert.equal(scoreMarkers({...markers,"BAT-26":false,"BAT-25":false}).msi_status, "MSS");
});
test("strict marker validation", () => {
  assert.throws(() => scoreMarkers({...markers,unknown:true}));
  assert.throws(() => scoreMarkers({...markers,"BAT-25":"true"}));
});
test("MLH1/PMS2 loss and methylation", () => {
  const result=interpretIHC({MLH1:false,MSH2:true,MSH6:true,PMS2:false},null,true);
  assert.equal(result.mmr_status,"dMMR");
  assert.equal(result.suspected_lynch_syndrome,false);
});
test("discordant result is not claimed concordant", () => {
  const r=combined(markers,{MLH1:true,MSH2:true,MSH6:true,PMS2:true});
  assert.equal(r.concordance,"Discordant");
  assert.equal(r.final_status,"Discordant (PCR/IHC)");
});
test("RFC-style quoted CSV fields", () => {
  assert.deepEqual(parseCSV('A,B\r\n"x,y","line\nbreak"\r\n'),[["A","B"],["x,y","line\nbreak"]]);
});
test("CSV processing rejects invalid or missing states", () => {
  const header="BAT-25,BAT-26,D2S123,D5S346,D17S250\n";
  assert.throws(()=>batchCSV(header+"true,wrong,false,false,false\n"),/invalid/);
  assert.throws(()=>batchCSV(header+"true,,false,false,false\n"),/invalid/);
  assert.equal(batchCSV(header+"true,true,false,false,false\n").count,1);
});
test("CSV output escapes spreadsheet formulas", () => {
  const r=batchCSV("BAT-25,BAT-26,D2S123,D5S346,D17S250,note\ntrue,true,false,false,false,=2+2\n");
  assert.match(r.csv, /'=2\+2/);
});
