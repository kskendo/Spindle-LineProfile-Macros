// ============================================================
// Macro: Per-cell full summary - length, amount, and shape
// ============================================================
// Description:
//   One row per cell, combining everything discussed:
//     - SpindleLength_um            : true pole-to-pole distance
//     - MeanCh1 / AUC_Ch1           : amount, image channel 1
//     - MeanCh2 / AUC_Ch2           : amount, image channel 2
//     - TubulinDipDepth             : SHAPE descriptor for the
//                                     tubulin channel (Ch2) - how
//                                     deep the central dip is,
//                                     relative to the two flanking
//                                     peaks, on that cell's own
//                                     per-spindle normalized (0-1)
//                                     profile. Large value = clear
//                                     two-peak/dip shape (like
//                                     siCTRL). Near zero/negative =
//                                     flat plateau (like siTPX2).
//     - KTPeakProminence             : SHAPE descriptor for the
//                                     poles+kinetochore channel
//                                     (Ch1) - how much taller the
//                                     central peak is than the
//                                     outer baseline, on the same
//                                     0-1 normalized profile. Large
//                                     value = sharp, well-clustered
//                                     central peak (like siCTRL).
//                                     Small value = broad/blunted
//                                     peak (like siTPX2).
//
//   Both shape descriptors are computed on the per-spindle
//   normalized (0-1) profile, calculated in memory - so they
//   capture SHAPE independent of that cell's absolute intensity,
//   matching how the shape was visualized in the normalized plots.
//
//   Fixed anatomical windows used (in NormPos_um):
//     center window  : -1 to +1 um
//     flank windows  : -5 to -3 um  and  +3 to +5 um  (tubulin)
//     edge windows   : -8 to -6 um  and  +6 to +8 um  (poles+KT)
//
//   TubulinDipDepth   = mean(flank windows, Ch2) - mean(center window, Ch2)
//   KTPeakProminence  = mean(center window, Ch1) - mean(edge windows, Ch1)
//
// Usage:
//   1. Place all RAW _lineprofile.csv files (all replicates, both
//      conditions) into a single folder. Remove any stray output
//      files (e.g. averaged_profiles_*.csv/.xls) from that folder
//      first - the macro will otherwise try to read them as cells.
//   2. Run via Plugins > Macros > Run.
//   3. Select that folder. One summary CSV is saved in the same
//      folder.
//
// Input:
//   - Folder of _lineprofile.csv files (standard 7-column structure)
//   - Filenames must contain "siCTRL"/"siControl" or "siTPX2"
//   - Each CSV must have exactly 161 valid data rows
//   - Files with "norm" in the name are skipped automatically
//
// Output (saved in same folder):
//   - percell_full_summary.csv
//   - Same table also pushed to the Results window / .xls
//
// Output CSV columns (one row per cell):
//   Filename, Condition, SpindleLength_um,
//   MeanCh1, AUC_Ch1, MeanCh2, AUC_Ch2,
//   TubulinDipDepth, KTPeakProminence
//
// Dependencies: Fiji/ImageJ (https://fiji.sc)
//
// Authors: Kristjana Skendo, Meraldi Lab, University of Geneva
// ============================================================

macro "Per-cell full summary - length, amount, shape" {

    dir = getDirectory("Choose the folder containing RAW _lineprofile.csv files");
    if (dir=="") exit("No folder selected.");

    nPts = 161;

    if (isOpen("Results")) { selectWindow("Results"); run("Close"); }

    print("\\Clear");
    print("Per-cell full summary for: " + dir);
    print("--------------------------------------------------");

    list = getFileList(dir);
    row = 0;

    out = "Filename,Condition,SpindleLength_um,MeanCh1,AUC_Ch1,MeanCh2,AUC_Ch2,TubulinDipDepth,KTPeakProminence\n";

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

        cond = "";
        if (indexOf(lower, "sitpx2") != -1)
            cond = "siTPX2";
        else if (indexOf(lower, "sictrl") != -1 || indexOf(lower, "sicontrol") != -1)
            cond = "siCTRL";
        else
            continue;

        path = dir + name;
        text = File.openAsString(path);

        if (lengthOf(text) == 0) {
            print("Skipped empty/unreadable file: " + path);
            continue;
        }

        lines = split(text, "\n");

        rd = newArray(nPts);
        fp = newArray(nPts); // NormPos_um
        v1 = newArray(nPts);
        v2 = newArray(nPts);
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

            normPos  = parseFloat(cols[1]);
            realDist = parseFloat(cols[2]);
            ch1      = parseFloat(cols[4]);
            ch2      = parseFloat(cols[6]);

            if (isNaN(normPos) || isNaN(realDist) || isNaN(ch1) || isNaN(ch2))
                continue;

            fp[idx] = normPos;
            rd[idx] = realDist;
            v1[idx] = ch1;
            v2[idx] = ch2;

            valid++;
        }

        if (valid != nPts) {
            print("Skipped (expected " + nPts + " valid rows, got " + valid + "): " + path);
            continue;
        }

        // ---- length ----
        spindleLength = rd[nPts-1] - rd[0];

        // ---- amount: mean + AUC (true um) for both channels ----
        sum1 = 0; sum2 = 0;
        for (i = 0; i < nPts; i++) { sum1 += v1[i]; sum2 += v2[i]; }
        mean1 = sum1 / nPts;
        mean2 = sum2 / nPts;

        auc1 = 0; auc2 = 0;
        for (i = 0; i < nPts-1; i++) {
            dx = rd[i+1] - rd[i];
            auc1 += (v1[i] + v1[i+1]) / 2.0 * dx;
            auc2 += (v2[i] + v2[i+1]) / 2.0 * dx;
        }

        // ---- per-spindle 0-1 normalization (in memory) ----
        min1 = v1[0]; max1 = v1[0];
        min2 = v2[0]; max2 = v2[0];
        for (i = 1; i < nPts; i++) {
            if (v1[i] < min1) min1 = v1[i];
            if (v1[i] > max1) max1 = v1[i];
            if (v2[i] < min2) min2 = v2[i];
            if (v2[i] > max2) max2 = v2[i];
        }
        range1 = max1 - min1; if (range1 == 0) range1 = 1;
        range2 = max2 - min2; if (range2 == 0) range2 = 1;

        n1 = newArray(nPts);
        n2 = newArray(nPts);
        for (i = 0; i < nPts; i++) {
            n1[i] = (v1[i] - min1) / range1;
            n2[i] = (v2[i] - min2) / range2;
        }

        // ---- shape windows (on NormPos_um, fixed anatomical positions) ----
        centerSum1 = 0; centerCount1 = 0;
        centerSum2 = 0; centerCount2 = 0;
        flankSum2  = 0; flankCount2  = 0;
        edgeSum1   = 0; edgeCount1   = 0;

        for (i = 0; i < nPts; i++) {

            p = fp[i];

            // center window: -1 to +1 um
            if (p >= -1 && p <= 1) {
                centerSum1 += n1[i]; centerCount1++;
                centerSum2 += n2[i]; centerCount2++;
            }

            // flank windows for tubulin (Ch2): -5..-3 and 3..5 um
            if ((p >= -5 && p <= -3) || (p >= 3 && p <= 5)) {
                flankSum2 += n2[i]; flankCount2++;
            }

            // edge windows for poles+KT (Ch1): -8..-6 and 6..8 um
            if ((p >= -8 && p <= -6) || (p >= 6 && p <= 8)) {
                edgeSum1 += n1[i]; edgeCount1++;
            }
        }

        centerMean1 = centerSum1 / centerCount1;
        centerMean2 = centerSum2 / centerCount2;
        flankMean2  = flankSum2  / flankCount2;
        edgeMean1   = edgeSum1   / edgeCount1;

        tubulinDipDepth  = flankMean2 - centerMean2;
        ktPeakProminence = centerMean1 - edgeMean1;

        out += name + "," + cond + "," + spindleLength + "," +
               mean1 + "," + auc1 + "," + mean2 + "," + auc2 + "," +
               tubulinDipDepth + "," + ktPeakProminence + "\n";

        setResult("Filename", row, name);
        setResult("Condition", row, cond);
        setResult("SpindleLength_um", row, spindleLength);
        setResult("MeanCh1", row, mean1);
        setResult("AUC_Ch1", row, auc1);
        setResult("MeanCh2", row, mean2);
        setResult("AUC_Ch2", row, auc2);
        setResult("TubulinDipDepth", row, tubulinDipDepth);
        setResult("KTPeakProminence", row, ktPeakProminence);
        row++;

        print("Used: " + name + "   [" + cond + "]");
    }

    if (row == 0)
        exit("No matching files found.");

    File.saveString(out, dir + "percell_full_summary.csv");
    updateResults();
    saveAs("Results", dir + "percell_full_summary.xls");

    print("--------------------------------------------------");
    print("DONE");
    print("Cells summarized = " + row);
    print("Saved: " + dir + "percell_full_summary.csv");
    print("Saved: " + dir + "percell_full_summary.xls");
}
