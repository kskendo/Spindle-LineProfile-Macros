// ============================================================
// Macro: Spindle line profile averager (DMSO / CTRL / 15NM)
// ============================================================
// Description:
//   Takes a folder of per-image line profile CSV files and computes
//   the grand mean and SEM per bin across three conditions:
//   DMSO, CTRL, and 15NM (taxol). Condition is detected automatically
//   from the filename (case-insensitive).
//
//   Run this macro TWICE — once on the raw folder and once on the
//   normalized folder — to produce both outputs:
//     - Raw averaged profiles (for statistical analysis)
//     - Normalized averaged profiles (for visualization)
//
// Usage:
//   1. Place all _lineprofile.csv (or _lineprofile_norm.csv) files
//      for one experiment in a single folder.
//   2. Run via Plugins > Macros > Run.
//   3. Select the folder containing the CSV files.
//   4. Output is saved in the same folder.
//
// Input:
//   - Folder of CSV files from macro_lineprofile-new_2channels.ijm
//     (raw: _lineprofile.csv) OR from Macro_normalize_perSpindle.ijm
//     (normalized: _lineprofile_norm.csv)
//   - Filenames must contain "DMSO", "CTRL", or "15NM" (case-insensitive)
//   - Each CSV must have exactly 161 valid data rows
//   - Columns read: col[1]=NormPos_um, col[4]=Ch1_bgsub, col[6]=Ch2_bgsub
//
// Output (saved in same folder):
//   - averaged_profiles_DMSO_CTRL_15NM.csv — mean and SEM per condition
//   - averaged_profiles_DMSO_CTRL_15NM.xls — same as Results table
//   - Two live preview plots (Ch1 and Ch2) shown in Fiji
//
// Output CSV columns:
//   NormPos_um,
//   DMSO_Ch1_mean, DMSO_Ch1_SEM, DMSO_n,
//   CTRL_Ch1_mean, CTRL_Ch1_SEM, CTRL_n,
//   15NM_Ch1_mean, 15NM_Ch1_SEM, 15NM_n,
//   DMSO_Ch2_mean, DMSO_Ch2_SEM, DMSO_n,
//   CTRL_Ch2_mean, CTRL_Ch2_SEM, CTRL_n,
//   15NM_Ch2_mean, 15NM_Ch2_SEM, 15NM_n
//
// Notes:
//   - SEM calculated using Bessel's correction (n-1)
//   - Files with fewer than 161 valid rows are skipped with a log message
//   - Files not matching any condition keyword are skipped
//
// Dependencies: Fiji/ImageJ (https://fiji.sc)
//
// Authors: Kristjana Skendo, Meraldi Lab, University of Geneva
// Publication: Skendo et al., 2025 (Cell Reports, title to be confirmed)
// Repository: https://github.com/kskendo/Spindle-LineProfile-Macros
// ============================================================

macro "Average line profiles from one folder - DMSO CTRL 15NM" {

    dir = getDirectory("Choose the folder containing all profile files");
    if (dir=="") exit("No folder selected.");

    nPts = 161;

    x = newArray(nPts);

    // DMSO
    dmso1_sum   = newArray(nPts); dmso1_sumsq   = newArray(nPts);
    dmso2_sum   = newArray(nPts); dmso2_sumsq   = newArray(nPts);

    // CTRL
    ctrl1_sum   = newArray(nPts); ctrl1_sumsq   = newArray(nPts);
    ctrl2_sum   = newArray(nPts); ctrl2_sumsq   = newArray(nPts);

    // 15NM
    nm151_sum   = newArray(nPts); nm151_sumsq   = newArray(nPts);
    nm152_sum   = newArray(nPts); nm152_sumsq   = newArray(nPts);

    dmso_n = 0;
    ctrl_n = 0;
    nm15_n = 0;
    totalFilesUsed = 0;

    print("\\Clear");
    print("Scanning folder: " + dir);
    print("Searching filenames for: DMSO, CTRL, 15NM");
    print("Case-insensitive, anywhere in filename");
    print("--------------------------------------------------");

    list = getFileList(dir);

    for (f = 0; f < list.length; f++) {

        name = list[f];

        if (endsWith(name, "/") || endsWith(name, "\\"))
            continue;

        lower = toLowerCase(name);

        // Accept CSV, TXT, and Fiji-style XLS text exports
        if (!(endsWith(lower, ".csv") || endsWith(lower, ".txt") || endsWith(lower, ".xls")))
            continue;

        cond = "";

        if (indexOf(lower, "dmso") != -1)
            cond = "dmso";
        else if (indexOf(lower, "ctrl") != -1)
            cond = "ctrl";
        else if (indexOf(lower, "15nm") != -1)
            cond = "15nm";
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

            xpos = parseFloat(cols[1]);  // NormPos_um
            v1   = parseFloat(cols[4]);  // Ch1_bgsub (raw) or Ch1_norm (normalized)
            v2   = parseFloat(cols[6]);  // Ch2_bgsub (raw) or Ch2_norm (normalized)

            if (isNaN(xpos) || isNaN(v1) || isNaN(v2))
                continue;

            x[idx] = xpos;

            if (cond == "dmso") {
                dmso1_sum[idx] += v1; dmso1_sumsq[idx] += v1*v1;
                dmso2_sum[idx] += v2; dmso2_sumsq[idx] += v2*v2;
            } else if (cond == "ctrl") {
                ctrl1_sum[idx] += v1; ctrl1_sumsq[idx] += v1*v1;
                ctrl2_sum[idx] += v2; ctrl2_sumsq[idx] += v2*v2;
            } else if (cond == "15nm") {
                nm151_sum[idx] += v1; nm151_sumsq[idx] += v1*v1;
                nm152_sum[idx] += v2; nm152_sumsq[idx] += v2*v2;
            }

            valid++;
        }

        if (valid == nPts) {
            totalFilesUsed++;
            if (cond == "dmso") dmso_n++;
            if (cond == "ctrl") ctrl_n++;
            if (cond == "15nm") nm15_n++;
            print("Used: " + path + "   [" + cond + "]");
        } else {
            print("Skipped (expected " + nPts + " valid rows, got " + valid + "): " + path);
        }
    }

    if (totalFilesUsed == 0)
        exit("No matching files found.");

    dmso1_mean = newArray(nPts); dmso1_sem = newArray(nPts);
    dmso2_mean = newArray(nPts); dmso2_sem = newArray(nPts);
    ctrl1_mean = newArray(nPts); ctrl1_sem = newArray(nPts);
    ctrl2_mean = newArray(nPts); ctrl2_sem = newArray(nPts);
    nm151_mean = newArray(nPts); nm151_sem = newArray(nPts);
    nm152_mean = newArray(nPts); nm152_sem = newArray(nPts);

    for (i = 0; i < nPts; i++) {

        if (dmso_n > 0) { dmso1_mean[i] = dmso1_sum[i]/dmso_n; dmso2_mean[i] = dmso2_sum[i]/dmso_n; }
        if (ctrl_n > 0) { ctrl1_mean[i] = ctrl1_sum[i]/ctrl_n; ctrl2_mean[i] = ctrl2_sum[i]/ctrl_n; }
        if (nm15_n > 0) { nm151_mean[i] = nm151_sum[i]/nm15_n; nm152_mean[i] = nm152_sum[i]/nm15_n; }

        if (dmso_n > 1) {
            mean = dmso1_sum[i]/dmso_n;
            variance = (dmso1_sumsq[i] - dmso_n*mean*mean)/(dmso_n-1);
            if (variance < 0) variance = 0;
            dmso1_sem[i] = sqrt(variance)/sqrt(dmso_n);
            mean = dmso2_sum[i]/dmso_n;
            variance = (dmso2_sumsq[i] - dmso_n*mean*mean)/(dmso_n-1);
            if (variance < 0) variance = 0;
            dmso2_sem[i] = sqrt(variance)/sqrt(dmso_n);
        }

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

        if (nm15_n > 1) {
            mean = nm151_sum[i]/nm15_n;
            variance = (nm151_sumsq[i] - nm15_n*mean*mean)/(nm15_n-1);
            if (variance < 0) variance = 0;
            nm151_sem[i] = sqrt(variance)/sqrt(nm15_n);
            mean = nm152_sum[i]/nm15_n;
            variance = (nm152_sumsq[i] - nm15_n*mean*mean)/(nm15_n-1);
            if (variance < 0) variance = 0;
            nm152_sem[i] = sqrt(variance)/sqrt(nm15_n);
        }
    }

    out = "NormPos_um," +
          "DMSO_Ch1_mean,DMSO_Ch1_SEM,DMSO_n," +
          "CTRL_Ch1_mean,CTRL_Ch1_SEM,CTRL_n," +
          "15NM_Ch1_mean,15NM_Ch1_SEM,15NM_n," +
          "DMSO_Ch2_mean,DMSO_Ch2_SEM,DMSO_n," +
          "CTRL_Ch2_mean,CTRL_Ch2_SEM,CTRL_n," +
          "15NM_Ch2_mean,15NM_Ch2_SEM,15NM_n\n";

    for (i = 0; i < nPts; i++) {
        out += x[i] + "," +
               dmso1_mean[i] + "," + dmso1_sem[i] + "," + dmso_n + "," +
               ctrl1_mean[i] + "," + ctrl1_sem[i] + "," + ctrl_n + "," +
               nm151_mean[i] + "," + nm151_sem[i] + "," + nm15_n + "," +
               dmso2_mean[i] + "," + dmso2_sem[i] + "," + dmso_n + "," +
               ctrl2_mean[i] + "," + ctrl2_sem[i] + "," + ctrl_n + "," +
               nm152_mean[i] + "," + nm152_sem[i] + "," + nm15_n + "\n";
    }

    File.saveString(out, dir + "averaged_profiles_DMSO_CTRL_15NM.csv");

    if (isOpen("Results")) { selectWindow("Results"); run("Close"); }

    for (i = 0; i < nPts; i++) {
        setResult("NormPos_um", i, x[i]);
        setResult("DMSO_Ch1_mean", i, dmso1_mean[i]); setResult("DMSO_Ch1_SEM", i, dmso1_sem[i]); setResult("DMSO_n", i, dmso_n);
        setResult("CTRL_Ch1_mean", i, ctrl1_mean[i]); setResult("CTRL_Ch1_SEM", i, ctrl1_sem[i]); setResult("CTRL_n", i, ctrl_n);
        setResult("15NM_Ch1_mean", i, nm151_mean[i]); setResult("15NM_Ch1_SEM", i, nm151_sem[i]); setResult("15NM_n", i, nm15_n);
        setResult("DMSO_Ch2_mean", i, dmso2_mean[i]); setResult("DMSO_Ch2_SEM", i, dmso2_sem[i]);
        setResult("CTRL_Ch2_mean", i, ctrl2_mean[i]); setResult("CTRL_Ch2_SEM", i, ctrl2_sem[i]);
        setResult("15NM_Ch2_mean", i, nm152_mean[i]); setResult("15NM_Ch2_SEM", i, nm152_sem[i]);
    }
    updateResults();
    saveAs("Results", dir + "averaged_profiles_DMSO_CTRL_15NM.xls");

    Plot.create("Channel 1 averaged profiles", "Normalized position (um)", "Ch1 intensity");
    if (dmso_n > 0) { Plot.setColor("black"); Plot.add("line", x, dmso1_mean); }
    if (ctrl_n > 0) { Plot.setColor("blue");  Plot.add("line", x, ctrl1_mean); }
    if (nm15_n > 0) { Plot.setColor("red");   Plot.add("line", x, nm151_mean); }
    Plot.show();

    Plot.create("Channel 2 averaged profiles", "Normalized position (um)", "Ch2 intensity");
    if (dmso_n > 0) { Plot.setColor("black"); Plot.add("line", x, dmso2_mean); }
    if (ctrl_n > 0) { Plot.setColor("blue");  Plot.add("line", x, ctrl2_mean); }
    if (nm15_n > 0) { Plot.setColor("red");   Plot.add("line", x, nm152_mean); }
    Plot.show();

    print("--------------------------------------------------");
    print("DONE");
    print("Files used = " + totalFilesUsed);
    print("DMSO n = " + dmso_n);
    print("CTRL  n = " + ctrl_n);
    print("15NM  n = " + nm15_n);
    print("Saved: " + dir + "averaged_profiles_DMSO_CTRL_15NM.csv");
    print("Saved: " + dir + "averaged_profiles_DMSO_CTRL_15NM.xls");
}
