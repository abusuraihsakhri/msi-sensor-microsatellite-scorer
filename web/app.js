"use strict";

document.addEventListener("DOMContentLoaded", () => {
  const form = document.getElementById("score-form");
  const result = document.getElementById("result");
  const ihcEnabled = document.getElementById("ihc-enabled");
  const ihcPanel = document.getElementById("ihc-panel");
  const batchFile = document.getElementById("batch-file");
  const batchInfo = document.getElementById("batch-info");
  const download = document.getElementById("download");
  let resultCSV = "";

  const showResult = (heading, detail, isError = false) => {
    result.textContent = heading + "\n\n" + detail;
    result.classList.toggle("error", isError);
    result.focus();
  };
  const updateIHC = () => {
    ihcPanel.hidden = !ihcEnabled.checked;
    ihcPanel.querySelectorAll("select").forEach(el => {
      el.disabled = !ihcEnabled.checked;
    });
  };
  ihcEnabled.addEventListener("change", updateIHC);
  updateIHC();

  form.addEventListener("submit", event => {
    event.preventDefault();
    try {
      const names = ["BAT-25", "BAT-26", "D2S123", "D5S346", "D17S250"];
      const markers = Object.fromEntries(names.map(name =>
        [name, form.elements.namedItem(name).value === "unstable"]));
      let proteins = null;
      let braf = null;
      let methylation = null;
      if (ihcEnabled.checked) {
        const proteinsNames = ["MLH1", "MSH2", "MSH6", "PMS2"];
        proteins = Object.fromEntries(proteinsNames.map(name =>
          [name, form.elements.namedItem(name).value === "retained"]));
        const optional = name => {
          const val = form.elements.namedItem(name).value;
          return val === "" ? null : val === "present";
        };
        braf = optional("BRAF_V600E");
        methylation = optional("MLH1_PROMOTER_METHYLATION");
      }
      const outcome = MSIScorer.combined(markers, proteins, braf, methylation);
      const parts = [
        "PCR (Bethesda five-marker panel): " + outcome.msi_pcr_result.msi_status,
        "Unstable loci: " + outcome.msi_pcr_result.unstable_count + "/5 (" + outcome.msi_pcr_result.percent_unstable + "%)"
      ];
      if (outcome.mmr_ihc_result) {
        parts.push("MMR IHC: " + outcome.mmr_ihc_result.mmr_status);
        parts.push("Lost proteins: " + (outcome.mmr_ihc_result.deficient_proteins.join(", ") || "None"));
        parts.push("Pattern: " + outcome.mmr_ihc_result.etiology);
        parts.push("PCR/IHC concordance: " + outcome.concordance);
      }
      if (outcome.concordance === "Discordant") {
        parts.push("The assays disagree. Further pathological review and confirmatory evaluation are necessary.");
      }
      showResult("Classification: " + outcome.final_status, parts.join("\n"));
    } catch (err) {
      showResult("Could not classify", err.message, true);
    }
  });

  form.addEventListener("reset", () => {
    setTimeout(() => {
      updateIHC();
      showResult("No assessment yet", "Select five PCR marker states and choose Calculate.");
    }, 0);
  });

  batchFile.addEventListener("change", async () => {
    resultCSV = "";
    download.disabled = true;
    batchInfo.textContent = "";
    const file = batchFile.files && batchFile.files[0];
    if (!file) return;
    try {
      if (file.size > 2000000) throw Error("CSV exceeds the 2 MB limit");
      const parsed = MSIScorer.batchCSV(await file.text());
      resultCSV = parsed.csv;
      batchInfo.textContent = parsed.count + " record(s) validated and classified locally. Download the processed CSV.";
      download.disabled = false;
    } catch (err) {
      batchInfo.textContent = "CSV error: " + err.message;
    }
  });

  download.addEventListener("click", () => {
    if (!resultCSV) return;
    const blob = new Blob([resultCSV], {type: "text/csv;charset=utf-8"});
    const url = URL.createObjectURL(blob);
    const link = document.createElement("a");
    link.href = url;
    link.download = "msi-results.csv";
    document.body.append(link);
    link.click();
    link.remove();
    setTimeout(() => URL.revokeObjectURL(url), 1000);
  });
});
