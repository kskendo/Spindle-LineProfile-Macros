// ============================================================
// Macro: Combined RAW + NORMALIZED averaging (siCTRL vs siTPX2)
// ============================================================
// Description:
//   One-pass version of the averaging pipeline. Reads a folder of
//   RAW _lineprofile.csv files only (no need to pre-run
//   Macro_normalize_perSpindle.ijm or create a separate normalized
//   folder) and produces BOTH:
//     - RAW averaged profiles (mean/SEM/N per condition) - for stats
//     - NORMALIZED averaged profiles (per-spindle 0-1, then
//       mean/SEM/N per condition) - for visualization
//   in a single run, both Prism-ready.
//
// Usage:
//   1. Place all RAW _lineprofile.csv files (all replicates, both
//      conditions) into a single folder.
//   2. Run via Plugins > Macros > Run.
//   3. Select that folder. Two output CSVs are saved in the same
//      folder.
//
// Input:
//   - Folder of _lineprofile.csv files (standard 7-column structure:
//     Image, NormPos_um, RealDist_um, Ch1_raw, Ch1_bgsub, Ch2_raw, Ch2_bgsub)
//   - Filenames must contain "siCTRL"/"siControl" or "siTPX2"
//     (case-insensitive, either spelling of control accepted)
//   - Each CSV must have exactly 161 valid data rows
//   - Files with "norm" in the name are skipped automatically
//     (in case a mixed folder is used by mistake)
//
// Output (saved in same folder):
//   - averaged_profiles_RAW_siCTRL_siTPX2.csv / .xls
//   - averaged_profiles_NORM_siCTRL_siTPX2.csv / .xls
//   - Four live preview plots (RAW Ch1, RAW Ch2, NORM Ch1, NORM Ch2)
//
// Output CSV columns (one row per NormPos_um bin, 161 rows):
//   NormPos_um,
//   siCTRL_Ch1_mean, siCTRL_Ch1_SEM, siCTRL_n,
//   siTPX2_Ch1_mean, siTPX2_Ch1_SEM, siTPX2_n,
//   siCTRL_Ch2_mean, siCTRL_Ch2_SEM, siCTRL_n,
//   siTPX2_Ch2_mean, siTPX2_Ch2_SEM, siTPX2_n
//
// Notes:
//   - Per-spindle normalization is min-max over that image's own
//     161 points, done in-memory here - identical logic to
//     Macro_normalize_perSpindle.ijm, just without writing
//     intermediate files to disk.
//   - N in both outputs = number of individual cells/images pooled
//     for that condition, across all replicates.
//   - SEM uses Bessel's correction (n-1).
//   - Files with fewer than 161 valid rows are skipped with a log message.
//
// Dependencies: Fiji/ImageJ (https://fiji.sc)
//
// Authors: Kristjana Skendo, Meraldi Lab, University of Geneva
// ============================================================

macro "Combined RAW+NORM averaging - siCTRL siTPX2" {

    dir = getDirectory("Choose the folder containing RAW _lineprofile.csv files");
    if (dir=="") exit("No folder selected.");

    nPts = 161;

    x = newArray(nPts);

    // RAW accumulators
    rCtrl1_sum = newArray(nPts); rCtrl1_sumsq = newArray(nPts);
    rCtrl2_sum = newArray(nPts); rCtrl2_sumsq = newArray(nPts);
    rTpx21_sum = newArray(nPts); rTpx21_sumsq = newArray(nPts);
    rTpx22_sum = newArray(nPts); rTpx22_sumsq = newArray(nPts);

    // NORM accumulators
    nCtrl1_sum = newArray(nPts); nCtrl1_sumsq = newArray(nPts);
    nCtrl2_sum = newArray(nPts); nCtrl2_sumsq = newArray(nPts);
    nTpx21_sum = newArray(nPts); nTpx21_sumsq = newArray(nPts);
    nTpx22_sum = newArray(nPts); nTpx22_sumsq = newArray(nPts);

    ctrl_n = 0;
    tpx2_n = 0;
    totalFilesUsed = 0;

    print("\\Clear");
    print("Combined RAW+NORM averaging in: " + dir);
    print("Searching filenames for: siCTRL / siControl, siTPX2");
    print("--------------------------------------------------");

    list = getFileList(dir);

    for (f = 0; f < list.length; f++) {

        name = list[f];

        if (endsWith(name, "/") || endsWith(name, "\\")) continue;

        lower = toLowerCase(name);

        if (indexOf(lower, "norm") != -1) {
            print("Skipped (looks like a normalized file): " + name);
            continue;
        }

        if (!(endsWith(lower, ".csv") || endsWith(lower, ".txt") || endsWith(lower, ".xls")))
            continue;

        cond = "";
        if (indexOf(lower, "sitpx2") != -1)
            cond = "sitpx2";
        else if (indexOf(lower, "sictrl") != -1 || indexOf(lower, "sicontrol") != -1)
            cond = "sictrl";
        else
            continue;

        path = dir + name;
        text = File.openAsString(path);

        if (lengthOf(text) == 0) {
            print("Skipped empty/unreadable file: " + path);
            continue;
        }

        lines = split(text, "\n");

        // First pass: collect this file's own values so we can
        // per-spindle normalize before accumulating into NORM sums.
        fx  = newArray(nPts);
        fv1 = newArray(nPts);
        fv2 = newArray(nPts);
        valid = 0;

        for (i = 1; i < lines.length; i++) {

            line = replace(lines[i], "\r", "");
            if (lengthOf(line) == 0) continue;

            if (indexOf(line, "\t") != -1)
                cols = split(line, "\t");
            else
                cols = split(line, ",");

            if (cols.length < 7) continue;

            idx = valid;
            if (idx >= nPts) break;

            xpos = parseFloat(cols[1]);
            v1   = parseFloat(cols[4]);
            v2   = parseFloat(cols[6]);

            if (isNaN(xpos) || isNaN(v1) || isNaN(v2))
                continue;

            fx[idx]  = xpos;
            fv1[idx] = v1;
            fv2[idx] = v2;

            valid++;
        }

        if (valid != nPts) {
            print("Skipped (expected " + nPts + " valid rows, got " + valid + "): " + path);
            continue;
        }

        // per-spindle min-max normalization for this file
        min1 = fv1[0]; max1 = fv1[0];
        min2 = fv2[0]; max2 = fv2[0];
        for (i = 1; i < nPts; i++) {
            if (fv1[i] < min1) min1 = fv1[i];
            if (fv1[i] > max1) max1 = fv1[i];
            if (fv2[i] < min2) min2 = fv2[i];
            if (fv2[i] > max2) max2 = fv2[i];
        }
        range1 = max1 - min1;
        range2 = max2 - min2;
        if (range1 == 0) range1 = 1;
        if (range2 == 0) range2 = 1;

        for (i = 0; i < nPts; i++) {

            x[i] = fx[i];

            v1n = (fv1[i] - min1) / range1;
            v2n = (fv2[i] - min2) / range2;

            if (cond == "sictrl") {
                rCtrl1_sum[i] += fv1[i]; rCtrl1_sumsq[i] += fv1[i]*fv1[i];
                rCtrl2_sum[i] += fv2[i]; rCtrl2_sumsq[i] += fv2[i]*fv2[i];
                nCtrl1_sum[i] += v1n;    nCtrl1_sumsq[i] += v1n*v1n;
                nCtrl2_sum[i] += v2n;    nCtrl2_sumsq[i] += v2n*v2n;
            } else if (cond == "sitpx2") {
                rTpx21_sum[i] += fv1[i]; rTpx21_sumsq[i] += fv1[i]*fv1[i];
                rTpx22_sum[i] += fv2[i]; rTpx22_sumsq[i] += fv2[i]*fv2[i];
                nTpx21_sum[i] += v1n;    nTpx21_sumsq[i] += v1n*v1n;
                nTpx22_sum[i] += v2n;    nTpx22_sumsq[i] += v2n*v2n;
            }
        }

        totalFilesUsed++;
        if (cond == "sictrl") ctrl_n++;
        if (cond == "sitpx2") tpx2_n++;
        print("Used: " + path + "   [" + cond + "]");
    }

    if (totalFilesUsed == 0)
        exit("No matching files found.");

    // ---- compute mean/SEM for RAW and NORM ----
    rCtrl1_mean = newArray(nPts); rCtrl1_sem = newArray(nPts);
    rCtrl2_mean = newArray(nPts); rCtrl2_sem = newArray(nPts);
    rTpx21_mean = newArray(nPts); rTpx21_sem = newArray(nPts);
    rTpx22_mean = newArray(nPts); rTpx22_sem = newArray(nPts);

    nCtrl1_mean = newArray(nPts); nCtrl1_sem = newArray(nPts);
    nCtrl2_mean = newArray(nPts); nCtrl2_sem = newArray(nPts);
    nTpx21_mean = newArray(nPts); nTpx21_sem = newArray(nPts);
    nTpx22_mean = newArray(nPts); nTpx22_sem = newArray(nPts);

    for (i = 0; i < nPts; i++) {

        if (ctrl_n > 0) {
            rCtrl1_mean[i] = rCtrl1_sum[i]/ctrl_n; rCtrl2_mean[i] = rCtrl2_sum[i]/ctrl_n;
            nCtrl1_mean[i] = nCtrl1_sum[i]/ctrl_n; nCtrl2_mean[i] = nCtrl2_sum[i]/ctrl_n;
        }
        if (tpx2_n > 0) {
            rTpx21_mean[i] = rTpx21_sum[i]/tpx2_n; rTpx22_mean[i] = rTpx22_sum[i]/tpx2_n;
            nTpx21_mean[i] = nTpx21_sum[i]/tpx2_n; nTpx22_mean[i] = nTpx22_sum[i]/tpx2_n;
        }

        if (ctrl_n > 1) {
            m = rCtrl1_sum[i]/ctrl_n; vr = (rCtrl1_sumsq[i]-ctrl_n*m*m)/(ctrl_n-1); if (vr<0) vr=0; rCtrl1_sem[i]=sqrt(vr)/sqrt(ctrl_n);
            m = rCtrl2_sum[i]/ctrl_n; vr = (rCtrl2_sumsq[i]-ctrl_n*m*m)/(ctrl_n-1); if (vr<0) vr=0; rCtrl2_sem[i]=sqrt(vr)/sqrt(ctrl_n);
            m = nCtrl1_sum[i]/ctrl_n; vr = (nCtrl1_sumsq[i]-ctrl_n*m*m)/(ctrl_n-1); if (vr<0) vr=0; nCtrl1_sem[i]=sqrt(vr)/sqrt(ctrl_n);
            m = nCtrl2_sum[i]/ctrl_n; vr = (nCtrl2_sumsq[i]-ctrl_n*m*m)/(ctrl_n-1); if (vr<0) vr=0; nCtrl2_sem[i]=sqrt(vr)/sqrt(ctrl_n);
        }

        if (tpx2_n > 1) {
            m = rTpx21_sum[i]/tpx2_n; vr = (rTpx21_sumsq[i]-tpx2_n*m*m)/(tpx2_n-1); if (vr<0) vr=0; rTpx21_sem[i]=sqrt(vr)/sqrt(tpx2_n);
            m = rTpx22_sum[i]/tpx2_n; vr = (rTpx22_sumsq[i]-tpx2_n*m*m)/(tpx2_n-1); if (vr<0) vr=0; rTpx22_sem[i]=sqrt(vr)/sqrt(tpx2_n);
            m = nTpx21_sum[i]/tpx2_n; vr = (nTpx21_sumsq[i]-tpx2_n*m*m)/(tpx2_n-1); if (vr<0) vr=0; nTpx21_sem[i]=sqrt(vr)/sqrt(tpx2_n);
            m = nTpx22_sum[i]/tpx2_n; vr = (nTpx22_sumsq[i]-tpx2_n*m*m)/(tpx2_n-1); if (vr<0) vr=0; nTpx22_sem[i]=sqrt(vr)/sqrt(tpx2_n);
        }
    }

    header = "NormPos_um," +
             "siCTRL_Ch1_mean,siCTRL_Ch1_SEM,siCTRL_n," +
             "siTPX2_Ch1_mean,siTPX2_Ch1_SEM,siTPX2_n," +
             "siCTRL_Ch2_mean,siCTRL_Ch2_SEM,siCTRL_n," +
             "siTPX2_Ch2_mean,siTPX2_Ch2_SEM,siTPX2_n\n";

    outRaw = header;
    outNorm = header;

    for (i = 0; i < nPts; i++) {
        outRaw += x[i] + "," +
                  rCtrl1_mean[i] + "," + rCtrl1_sem[i] + "," + ctrl_n + "," +
                  rTpx21_mean[i] + "," + rTpx21_sem[i] + "," + tpx2_n + "," +
                  rCtrl2_mean[i] + "," + rCtrl2_sem[i] + "," + ctrl_n + "," +
                  rTpx22_mean[i] + "," + rTpx22_sem[i] + "," + tpx2_n + "\n";

        outNorm += x[i] + "," +
                   nCtrl1_mean[i] + "," + nCtrl1_sem[i] + "," + ctrl_n + "," +
                   nTpx21_mean[i] + "," + nTpx21_sem[i] + "," + tpx2_n + "," +
                   nCtrl2_mean[i] + "," + nCtrl2_sem[i] + "," + ctrl_n + "," +
                   nTpx22_mean[i] + "," + nTpx22_sem[i] + "," + tpx2_n + "\n";
    }

    File.saveString(outRaw,  dir + "averaged_profiles_RAW_siCTRL_siTPX2.csv");
    File.saveString(outNorm, dir + "averaged_profiles_NORM_siCTRL_siTPX2.csv");

    // ---- RAW Results table + xls ----
    if (isOpen("Results")) { selectWindow("Results"); run("Close"); }
    for (i = 0; i < nPts; i++) {
        setResult("NormPos_um", i, x[i]);
        setResult("siCTRL_Ch1_mean", i, rCtrl1_mean[i]); setResult("siCTRL_Ch1_SEM", i, rCtrl1_sem[i]); setResult("siCTRL_n", i, ctrl_n);
        setResult("siTPX2_Ch1_mean", i, rTpx21_mean[i]); setResult("siTPX2_Ch1_SEM", i, rTpx21_sem[i]); setResult("siTPX2_n", i, tpx2_n);
        setResult("siCTRL_Ch2_mean", i, rCtrl2_mean[i]); setResult("siCTRL_Ch2_SEM", i, rCtrl2_sem[i]);
        setResult("siTPX2_Ch2_mean", i, rTpx22_mean[i]); setResult("siTPX2_Ch2_SEM", i, rTpx22_sem[i]);
    }
    updateResults();
    saveAs("Results", dir + "averaged_profiles_RAW_siCTRL_siTPX2.xls");

    // ---- NORM Results table + xls ----
    if (isOpen("Results")) { selectWindow("Results"); run("Close"); }
    for (i = 0; i < nPts; i++) {
        setResult("NormPos_um", i, x[i]);
        setResult("siCTRL_Ch1_mean", i, nCtrl1_mean[i]); setResult("siCTRL_Ch1_SEM", i, nCtrl1_sem[i]); setResult("siCTRL_n", i, ctrl_n);
        setResult("siTPX2_Ch1_mean", i, nTpx21_mean[i]); setResult("siTPX2_Ch1_SEM", i, nTpx21_sem[i]); setResult("siTPX2_n", i, tpx2_n);
        setResult("siCTRL_Ch2_mean", i, nCtrl2_mean[i]); setResult("siCTRL_Ch2_SEM", i, nCtrl2_sem[i]);
        setResult("siTPX2_Ch2_mean", i, nTpx22_mean[i]); setResult("siTPX2_Ch2_SEM", i, nTpx22_sem[i]);
    }
    updateResults();
    saveAs("Results", dir + "averaged_profiles_NORM_siCTRL_siTPX2.xls");

    // ---- plots ----
    Plot.create("RAW Channel 1 averaged profiles", "Normalized position (um)", "Intensity (bgsub)");
    if (ctrl_n > 0) { Plot.setColor("black"); Plot.add("line", x, rCtrl1_mean); }
    if (tpx2_n > 0) { Plot.setColor("red");   Plot.add("line", x, rTpx21_mean); }
    Plot.show();

    Plot.create("RAW Channel 2 averaged profiles", "Normalized position (um)", "Intensity (bgsub)");
    if (ctrl_n > 0) { Plot.setColor("black"); Plot.add("line", x, rCtrl2_mean); }
    if (tpx2_n > 0) { Plot.setColor("red");   Plot.add("line", x, rTpx22_mean); }
    Plot.show();

    Plot.create("NORM Channel 1 averaged profiles", "Normalized position (um)", "Intensity (0-1 norm)");
    if (ctrl_n > 0) { Plot.setColor("black"); Plot.add("line", x, nCtrl1_mean); }
    if (tpx2_n > 0) { Plot.setColor("red");   Plot.add("line", x, nTpx21_mean); }
    Plot.show();

    Plot.create("NORM Channel 2 averaged profiles", "Normalized position (um)", "Intensity (0-1 norm)");
    if (ctrl_n > 0) { Plot.setColor("black"); Plot.add("line", x, nCtrl2_mean); }
    if (tpx2_n > 0) { Plot.setColor("red");   Plot.add("line", x, nTpx22_mean); }
    Plot.show();

    print("--------------------------------------------------");
    print("DONE");
    print("Files used = " + totalFilesUsed);
    print("siCTRL n = " + ctrl_n);
    print("siTPX2 n = " + tpx2_n);
    print("Saved: " + dir + "averaged_profiles_RAW_siCTRL_siTPX2.csv");
    print("Saved: " + dir + "averaged_profiles_NORM_siCTRL_siTPX2.csv");
}
