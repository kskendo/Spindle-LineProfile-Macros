// ============================================================
// Macro: Combined RAW + NORMALIZED averaging
// HAUS6 / HURP / Poles+Kinetochores - DMSO vs Importazole
// ============================================================
// Description:
//   Position-by-position averaging for the 3-channel
//   HAUS6/HURP/PolesKT dataset. Reads a folder of RAW
//   _lineprofile.csv files (9-column structure: Ch1=HAUS6,
//   Ch2=HURP, Ch3=Poles+kinetochores) and produces BOTH:
//     - RAW averaged profiles (mean/SEM/N per condition)  - stats
//     - NORMALIZED averaged profiles (per-spindle 0-1, then
//       mean/SEM/N per condition) - visualization
//   for DMSO vs Importazole, all three channels, in one run.
//
//   Same condition-matching logic as
//   PerCell_HAUS6_HURP_PolesKT_Summary.ijm: "dmso" -> DMSO,
//   "importazole" or "ipz" -> Importazole.
//
// Usage:
//   1. Place all RAW _lineprofile.csv files (both conditions, all
//      replicates) into a single folder. Remove any stray output
//      files first.
//   2. Run via Plugins > Macros > Run.
//   3. Select that folder. Two output CSVs are saved in the same
//      folder.
//
// Input:
//   - Folder of _lineprofile.csv files with the 9-column structure:
//     Image, NormPos_um, RealDist_um,
//     Ch1_raw, Ch1_bgsub, Ch2_raw, Ch2_bgsub, Ch3_raw, Ch3_bgsub
//   - Each CSV must have exactly 161 valid data rows
//   - Files with "norm" in the name, or pipeline output files
//     (averaged_profiles_*, percell_*), are skipped automatically
//
// Output (saved in same folder):
//   - averaged_profiles_RAW_HAUS6_HURP_PolesKT.csv
//   - averaged_profiles_NORM_HAUS6_HURP_PolesKT.csv
//   (no live preview plots - build comparison plots in Prism from
//   the CSV, same as the multi-condition tubulin pipeline)
//
// Output CSV columns (one row per NormPos_um bin, 161 rows),
// grouped by channel (all Ch1 columns, then Ch2, then Ch3):
//   NormPos_um,
//   DMSO_Ch1_mean, DMSO_Ch1_SEM, DMSO_n,
//   Importazole_Ch1_mean, Importazole_Ch1_SEM, Importazole_n,
//   DMSO_Ch2_mean, DMSO_Ch2_SEM, DMSO_n,
//   Importazole_Ch2_mean, Importazole_Ch2_SEM, Importazole_n,
//   DMSO_Ch3_mean, DMSO_Ch3_SEM, DMSO_n,
//   Importazole_Ch3_mean, Importazole_Ch3_SEM, Importazole_n
//
// Notes:
//   - N = number of individual cells pooled for that condition,
//     across all replicates.
//   - SEM uses Bessel's correction (n-1).
//   - This position-by-position averaging is for VISUALIZATION
//     only. For statistics, use one row per cell
//     (PerCell_HAUS6_HURP_PolesKT_Summary.ijm output) with an
//     unpaired t-test (or Welch's), not this table.
//
// Dependencies: Fiji/ImageJ (https://fiji.sc)
//
// Authors: Kristjana Skendo, Meraldi Lab, University of Geneva
// ============================================================

macro "Combined RAW+NORM averaging - HAUS6/HURP/PolesKT" {

    dir = getDirectory("Choose the folder containing RAW _lineprofile.csv files");
    if (dir=="") exit("No folder selected.");

    nPts = 161;

    condNames = newArray("DMSO","Importazole");
    nCond = condNames.length;
    totalSize = nCond * nPts;

    x = newArray(nPts);

    // RAW accumulators (flattened: index = condIndex*nPts + i), one set per channel
    rSum1 = newArray(totalSize); rSumsq1 = newArray(totalSize);
    rSum2 = newArray(totalSize); rSumsq2 = newArray(totalSize);
    rSum3 = newArray(totalSize); rSumsq3 = newArray(totalSize);

    // NORM accumulators
    nmSum1 = newArray(totalSize); nmSumsq1 = newArray(totalSize);
    nmSum2 = newArray(totalSize); nmSumsq2 = newArray(totalSize);
    nmSum3 = newArray(totalSize); nmSumsq3 = newArray(totalSize);

    condN = newArray(nCond);

    totalFilesUsed = 0;

    print("\\Clear");
    print("HAUS6/HURP/PolesKT RAW+NORM averaging in: " + dir);
    print("Conditions: DMSO, Importazole");
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
        if (indexOf(lower, "averaged_profiles") != -1 || indexOf(lower, "percell") != -1) {
            print("Skipped (pipeline output file, not a cell): " + name);
            continue;
        }
        if (!(endsWith(lower, ".csv") || endsWith(lower, ".txt") || endsWith(lower, ".xls")))
            continue;

        cond = classifyCondition(lower);
        if (cond == "") {
            print("Skipped (condition not recognized in filename): " + name);
            continue;
        }

        condIndex = -1;
        for (ci = 0; ci < nCond; ci++) {
            if (condNames[ci] == cond) { condIndex = ci; ci = nCond; }
        }
        if (condIndex == -1) {
            print("Skipped (condition label mismatch): " + name);
            continue;
        }

        path = dir + name;
        text = File.openAsString(path);
        if (lengthOf(text) == 0) {
            print("Skipped empty/unreadable file: " + path);
            continue;
        }

        lines = split(text, "\n");

        fx  = newArray(nPts);
        fv1 = newArray(nPts);
        fv2 = newArray(nPts);
        fv3 = newArray(nPts);
        valid = 0;

        for (i = 1; i < lines.length; i++) {

            line = replace(lines[i], "\r", "");
            if (lengthOf(line) == 0) continue;

            if (indexOf(line, "\t") != -1)
                cols = split(line, "\t");
            else
                cols = split(line, ",");

            if (cols.length < 9) continue;

            idx = valid;
            if (idx >= nPts) break;

            xpos = parseFloat(cols[1]);
            v1   = parseFloat(cols[4]);
            v2   = parseFloat(cols[6]);
            v3   = parseFloat(cols[8]);

            if (isNaN(xpos) || isNaN(v1) || isNaN(v2) || isNaN(v3))
                continue;

            fx[idx]  = xpos;
            fv1[idx] = v1;
            fv2[idx] = v2;
            fv3[idx] = v3;

            valid++;
        }

        if (valid != nPts) {
            print("Skipped (expected " + nPts + " valid rows, got " + valid + "): " + path);
            continue;
        }

        // per-spindle min-max normalization for this file, all 3 channels
        min1 = fv1[0]; max1 = fv1[0];
        min2 = fv2[0]; max2 = fv2[0];
        min3 = fv3[0]; max3 = fv3[0];
        for (i = 1; i < nPts; i++) {
            if (fv1[i] < min1) min1 = fv1[i];
            if (fv1[i] > max1) max1 = fv1[i];
            if (fv2[i] < min2) min2 = fv2[i];
            if (fv2[i] > max2) max2 = fv2[i];
            if (fv3[i] < min3) min3 = fv3[i];
            if (fv3[i] > max3) max3 = fv3[i];
        }
        range1 = max1 - min1; if (range1 == 0) range1 = 1;
        range2 = max2 - min2; if (range2 == 0) range2 = 1;
        range3 = max3 - min3; if (range3 == 0) range3 = 1;

        for (i = 0; i < nPts; i++) {

            x[i] = fx[i];
            flatIdx = condIndex * nPts + i;

            v1n = (fv1[i] - min1) / range1;
            v2n = (fv2[i] - min2) / range2;
            v3n = (fv3[i] - min3) / range3;

            rSum1[flatIdx] += fv1[i]; rSumsq1[flatIdx] += fv1[i]*fv1[i];
            rSum2[flatIdx] += fv2[i]; rSumsq2[flatIdx] += fv2[i]*fv2[i];
            rSum3[flatIdx] += fv3[i]; rSumsq3[flatIdx] += fv3[i]*fv3[i];

            nmSum1[flatIdx] += v1n; nmSumsq1[flatIdx] += v1n*v1n;
            nmSum2[flatIdx] += v2n; nmSumsq2[flatIdx] += v2n*v2n;
            nmSum3[flatIdx] += v3n; nmSumsq3[flatIdx] += v3n*v3n;
        }

        condN[condIndex]++;
        totalFilesUsed++;
        print("Used: " + name + "   [" + cond + "]");
    }

    if (totalFilesUsed == 0)
        exit("No matching files found.");

    print("--------------------------------------------------");
    for (ci = 0; ci < nCond; ci++) {
        print(condNames[ci] + ": n = " + condN[ci]);
    }

    // ---- build headers (grouped by channel: all Ch1, then Ch2, then Ch3) ----
    headerRaw  = "NormPos_um";
    headerNorm = "NormPos_um";
    for (chLabel = 1; chLabel <= 3; chLabel++) {
        for (ci = 0; ci < nCond; ci++) {
            headerRaw  = headerRaw  + "," + condNames[ci] + "_Ch" + chLabel + "_mean," +
                         condNames[ci] + "_Ch" + chLabel + "_SEM," + condNames[ci] + "_n";
            headerNorm = headerNorm + "," + condNames[ci] + "_Ch" + chLabel + "_mean," +
                         condNames[ci] + "_Ch" + chLabel + "_SEM," + condNames[ci] + "_n";
        }
    }
    outRaw  = headerRaw  + "\n";
    outNorm = headerNorm + "\n";

    // ---- compute mean/SEM per bin and build rows ----
    for (i = 0; i < nPts; i++) {

        rowRaw  = "" + x[i];
        rowNorm = "" + x[i];

        for (chLabel = 1; chLabel <= 3; chLabel++) {
            for (ci = 0; ci < nCond; ci++) {

                flatIdx = ci * nPts + i;
                n = condN[ci];

                if (chLabel == 1) {
                    rSum = rSum1[flatIdx]; rSumsq = rSumsq1[flatIdx];
                    nmSum = nmSum1[flatIdx]; nmSumsq = nmSumsq1[flatIdx];
                } else if (chLabel == 2) {
                    rSum = rSum2[flatIdx]; rSumsq = rSumsq2[flatIdx];
                    nmSum = nmSum2[flatIdx]; nmSumsq = nmSumsq2[flatIdx];
                } else {
                    rSum = rSum3[flatIdx]; rSumsq = rSumsq3[flatIdx];
                    nmSum = nmSum3[flatIdx]; nmSumsq = nmSumsq3[flatIdx];
                }

                rMean = 0; rSem = 0; nmMean = 0; nmSem = 0;
                if (n > 0) {
                    rMean  = rSum  / n;
                    nmMean = nmSum / n;
                }
                if (n > 1) {
                    m = rSum / n;
                    vr = (rSumsq - n*m*m) / (n-1); if (vr < 0) vr = 0;
                    rSem = sqrt(vr) / sqrt(n);

                    m = nmSum / n;
                    vr = (nmSumsq - n*m*m) / (n-1); if (vr < 0) vr = 0;
                    nmSem = sqrt(vr) / sqrt(n);
                }

                rowRaw  = rowRaw  + "," + rMean  + "," + rSem  + "," + n;
                rowNorm = rowNorm + "," + nmMean + "," + nmSem + "," + n;
            }
        }

        outRaw  = outRaw  + rowRaw  + "\n";
        outNorm = outNorm + rowNorm + "\n";
    }

    File.saveString(outRaw,  dir + "averaged_profiles_RAW_HAUS6_HURP_PolesKT.csv");
    File.saveString(outNorm, dir + "averaged_profiles_NORM_HAUS6_HURP_PolesKT.csv");

    print("--------------------------------------------------");
    print("DONE");
    print("Files used = " + totalFilesUsed);
    print("Saved: " + dir + "averaged_profiles_RAW_HAUS6_HURP_PolesKT.csv");
    print("Saved: " + dir + "averaged_profiles_NORM_HAUS6_HURP_PolesKT.csv");
}

function classifyCondition(lower) {
    if (indexOf(lower, "importazole") != -1 || indexOf(lower, "ipz") != -1)
        return "Importazole";
    if (indexOf(lower, "dmso") != -1)
        return "DMSO";
    return "";
}
