// ============================================================
// Macro: Average NORMALIZED line profiles (siCTRL vs siTPX2) - Prism-ready
// ============================================================
// Description:
//   Averages per-spindle NORMALIZED (0-1) line profile CSVs
//   across replicates, grouped by condition (siCTRL / siTPX2),
//   and outputs one row per position bin with mean, SEM, and N
//   per condition per channel — ready to paste straight into a
//   Prism XY table (one X column, then Mean/SEM/N per data set).
//
//   Use the NORMALIZED output for visualization/plotting shape
//   comparisons between conditions (not for statistics — use the
//   RAW averaging macro for that).
//
// Usage:
//   1. First run Macro_normalize_perSpindle.ijm on your folder of
//      RAW _lineprofile.csv files to generate _lineprofile_norm.csv
//      files (one per image, saved alongside the raw ones).
//   2. Place all _lineprofile_norm.csv files for one experiment,
//      across all replicates, into a single folder.
//   3. Run via Plugins > Macros > Run.
//   4. Select that folder. Output is saved in the same folder.
//
// Input:
//   - Folder of _lineprofile_norm.csv files
//   - Filenames must contain "siCTRL" or "siTPX2" (case-insensitive)
//   - Each CSV must have exactly 161 valid data rows
//   - Columns read: col[1]=NormPos_um, col[4]=Ch1 (0-1 normalized),
//     col[6]=Ch2 (0-1 normalized)
//
// Output (saved in same folder):
//   - averaged_profiles_NORM_siCTRL_siTPX2.csv
//   - averaged_profiles_NORM_siCTRL_siTPX2.xls (same as Results table)
//   - Two live preview plots (Ch1 and Ch2) shown in Fiji
//
// Output CSV columns (one row per NormPos_um bin, 161 rows):
//   NormPos_um,
//   siCTRL_Ch1_mean, siCTRL_Ch1_SEM, siCTRL_n,
//   siTPX2_Ch1_mean, siTPX2_Ch1_SEM, siTPX2_n,
//   siCTRL_Ch2_mean, siCTRL_Ch2_SEM, siCTRL_n,
//   siTPX2_Ch2_mean, siTPX2_Ch2_SEM, siTPX2_n
//
// Notes:
//   - SEM calculated using Bessel's correction (n-1)
//   - Files with fewer than 161 valid rows are skipped with a log message
//   - Files not matching "siCTRL" or "siTPX2" are skipped
//   - To import into Prism: paste NormPos_um as X, then the
//     matching Mean/SEM/N triplet for each data set (Ch1 siCTRL,
//     Ch1 siTPX2, etc.) as separate Y data sets in "Mean, SEM, N" format
//
// Dependencies: Fiji/ImageJ (https://fiji.sc)
//
// Authors: Kristjana Skendo, Meraldi Lab, University of Geneva
// ============================================================

macro "Average NORMALIZED line profiles - siCTRL siTPX2" {

    dir = getDirectory("Choose the folder containing NORMALIZED _lineprofile_norm.csv files");
    if (dir=="") exit("No folder selected.");

    nPts = 161;

    x = newArray(nPts);

    ctrl1_sum = newArray(nPts); ctrl1_sumsq = newArray(nPts);
    ctrl2_sum = newArray(nPts); ctrl2_sumsq = newArray(nPts);

    tpx21_sum = newArray(nPts); tpx21_sumsq = newArray(nPts);
    tpx22_sum = newArray(nPts); tpx22_sumsq = newArray(nPts);

    ctrl_n = 0;
    tpx2_n = 0;
    totalFilesUsed = 0;

    print("\\Clear");
    print("Averaging NORMALIZED profiles in: " + dir);
    print("Searching filenames for: siCTRL, siTPX2");
    print("--------------------------------------------------");

    list = getFileList(dir);

    for (f = 0; f < list.length; f++) {

        name = list[f];

        if (endsWith(name, "/") || endsWith(name, "\\")) continue;

        lower = toLowerCase(name);

        // Only accept files that look like normalized output
        if (indexOf(lower, "norm") == -1) {
            print("Skipped (does not look like a normalized file): " + name);
            continue;
        }

        if (!(endsWith(lower, ".csv") || endsWith(lower, ".txt") || endsWith(lower, ".xls")))
            continue;

        cond = "";
        if (indexOf(lower, "sitpx2") != -1)
            cond = "sitpx2";
        else if (indexOf(lower, "sictrl") != -1)
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

            x[idx] = xpos;

            if (cond == "sictrl") {
                ctrl1_sum[idx] += v1; ctrl1_sumsq[idx] += v1*v1;
                ctrl2_sum[idx] += v2; ctrl2_sumsq[idx] += v2*v2;
            } else if (cond == "sitpx2") {
                tpx21_sum[idx] += v1; tpx21_sumsq[idx] += v1*v1;
                tpx22_sum[idx] += v2; tpx22_sumsq[idx] += v2*v2;
            }

            valid++;
        }

        if (valid == nPts) {
            totalFilesUsed++;
            if (cond == "sictrl") ctrl_n++;
            if (cond == "sitpx2") tpx2_n++;
            print("Used: " + path + "   [" + cond + "]");
        } else {
            print("Skipped (expected " + nPts + " valid rows, got " + valid + "): " + path);
        }
    }

    if (totalFilesUsed == 0)
        exit("No matching files found.");

    ctrl1_mean = newArray(nPts); ctrl1_sem = newArray(nPts);
    ctrl2_mean = newArray(nPts); ctrl2_sem = newArray(nPts);
    tpx21_mean = newArray(nPts); tpx21_sem = newArray(nPts);
    tpx22_mean = newArray(nPts); tpx22_sem = newArray(nPts);

    for (i = 0; i < nPts; i++) {

        if (ctrl_n > 0) { ctrl1_mean[i] = ctrl1_sum[i]/ctrl_n; ctrl2_mean[i] = ctrl2_sum[i]/ctrl_n; }
        if (tpx2_n > 0) { tpx21_mean[i] = tpx21_sum[i]/tpx2_n; tpx22_mean[i] = tpx22_sum[i]/tpx2_n; }

        if (ctrl_n > 1) {
            mean = ctrl1_sum[i]/ctrl_n;
            variance = (ctrl1_sumsq[i] - ctrl_n*mean*mean)/(ctrl_n-1);
            if (variance < 0) variance = 0;
            ctrl1_sem[i] = sqrt(variance)/sqrt(ctrl_n);
            mean = ctrl2_sum[i]/ctrl_n;
            variance = (ctrl2_sumsq[i] - ctrl_n*mean*mean)/(ctrl_n-1);
            if (variance < 0) variance = 0;
            ctrl2_sem[i] = sqrt(variance)/sqrt(ctrl_n);
        }

        if (tpx2_n > 1) {
            mean = tpx21_sum[i]/tpx2_n;
            variance = (tpx21_sumsq[i] - tpx2_n*mean*mean)/(tpx2_n-1);
            if (variance < 0) variance = 0;
            tpx21_sem[i] = sqrt(variance)/sqrt(tpx2_n);
            mean = tpx22_sum[i]/tpx2_n;
            variance = (tpx22_sumsq[i] - tpx2_n*mean*mean)/(tpx2_n-1);
            if (variance < 0) variance = 0;
            tpx22_sem[i] = sqrt(variance)/sqrt(tpx2_n);
        }
    }

    out = "NormPos_um," +
          "siCTRL_Ch1_mean,siCTRL_Ch1_SEM,siCTRL_n," +
          "siTPX2_Ch1_mean,siTPX2_Ch1_SEM,siTPX2_n," +
          "siCTRL_Ch2_mean,siCTRL_Ch2_SEM,siCTRL_n," +
          "siTPX2_Ch2_mean,siTPX2_Ch2_SEM,siTPX2_n\n";

    for (i = 0; i < nPts; i++) {
        out += x[i] + "," +
               ctrl1_mean[i] + "," + ctrl1_sem[i] + "," + ctrl_n + "," +
               tpx21_mean[i] + "," + tpx21_sem[i] + "," + tpx2_n + "," +
               ctrl2_mean[i] + "," + ctrl2_sem[i] + "," + ctrl_n + "," +
               tpx22_mean[i] + "," + tpx22_sem[i] + "," + tpx2_n + "\n";
    }

    File.saveString(out, dir + "averaged_profiles_NORM_siCTRL_siTPX2.csv");

    if (isOpen("Results")) { selectWindow("Results"); run("Close"); }

    for (i = 0; i < nPts; i++) {
        setResult("NormPos_um", i, x[i]);
        setResult("siCTRL_Ch1_mean", i, ctrl1_mean[i]); setResult("siCTRL_Ch1_SEM", i, ctrl1_sem[i]); setResult("siCTRL_n", i, ctrl_n);
        setResult("siTPX2_Ch1_mean", i, tpx21_mean[i]); setResult("siTPX2_Ch1_SEM", i, tpx21_sem[i]); setResult("siTPX2_n", i, tpx2_n);
        setResult("siCTRL_Ch2_mean", i, ctrl2_mean[i]); setResult("siCTRL_Ch2_SEM", i, ctrl2_sem[i]);
        setResult("siTPX2_Ch2_mean", i, tpx22_mean[i]); setResult("siTPX2_Ch2_SEM", i, tpx22_sem[i]);
    }
    updateResults();
    saveAs("Results", dir + "averaged_profiles_NORM_siCTRL_siTPX2.xls");

    Plot.create("NORM Channel 1 averaged profiles", "Normalized position (um)", "Intensity (0-1 norm)");
    if (ctrl_n > 0) { Plot.setColor("black"); Plot.add("line", x, ctrl1_mean); }
    if (tpx2_n > 0) { Plot.setColor("red");   Plot.add("line", x, tpx21_mean); }
    Plot.show();

    Plot.create("NORM Channel 2 averaged profiles", "Normalized position (um)", "Intensity (0-1 norm)");
    if (ctrl_n > 0) { Plot.setColor("black"); Plot.add("line", x, ctrl2_mean); }
    if (tpx2_n > 0) { Plot.setColor("red");   Plot.add("line", x, tpx22_mean); }
    Plot.show();

    print("--------------------------------------------------");
    print("DONE (NORMALIZED)");
    print("Files used = " + totalFilesUsed);
    print("siCTRL n = " + ctrl_n);
    print("siTPX2 n = " + tpx2_n);
    print("Saved: " + dir + "averaged_profiles_NORM_siCTRL_siTPX2.csv");
    print("Saved: " + dir + "averaged_profiles_NORM_siCTRL_siTPX2.xls");
}
