// ============================================================
// Macro: Per-spindle normalizer
// ============================================================
// Description:
//   Takes a folder of per-image line profile CSV files produced
//   by "macro_lineprofile-new_2channels.ijm" and generates a
//   normalized version of each file.
//
//   For each CSV file, Ch1_bgsub and Ch2_bgsub are min-max
//   normalized independently to [0, 1] across the 161 bins
//   of that spindle. All other columns are preserved unchanged.
//
//   This allows you to run Macro_averaging_.ijm separately on:
//     - The raw folder      → averaged raw profiles (for statistics)
//     - The norm folder     → averaged normalized profiles (for visualization)
//
// Usage:
//   1. Run "macro_lineprofile-new_2channels.ijm" on all images
//      to generate per-image _lineprofile.csv files.
//   2. Place all _lineprofile.csv files in one folder.
//   3. Run this macro via Plugins > Macros > Run.
//   4. Select the folder containing the _lineprofile.csv files.
//   5. A _lineprofile_norm.csv is saved alongside each input file.
//
// Input:
//   - Folder of CSV files from macro_lineprofile-new_2channels.ijm
//   - Each CSV must have 7 columns:
//     Image, NormPos_um, RealDist_um, Ch1_raw, Ch1_bgsub, Ch2_raw, Ch2_bgsub
//   - Only files ending in _lineprofile.csv are processed
//
// Output:
//   - <imagename>_lineprofile_norm.csv — same 7-column structure,
//     with Ch1_bgsub and Ch2_bgsub replaced by per-spindle
//     min-max normalized values (0 to 1)
//
// Notes:
//   - Normalization is per-file (per spindle), not global
//   - Ch1_raw and Ch2_raw columns are preserved unchanged
//   - Files not ending in _lineprofile.csv are skipped
//   - Already-normalized files (_lineprofile_norm.csv) are skipped
//
// Dependencies: Fiji/ImageJ (https://fiji.sc)
//
// Authors: Kristjana Skendo, Meraldi Lab, University of Geneva
// Publication: Skendo et al., 2025 (Cell Reports, title to be confirmed)
// Repository: https://github.com/kskendo/Spindle-LineProfile-Macros
// ============================================================

macro "Per-spindle min-max normalizer" {

    dir = getDirectory("Choose folder containing _lineprofile.csv files");
    if (dir == "") exit("No folder selected.");

    list = getFileList(dir);
    processed = 0;
    skipped = 0;

    print("\\Clear");
    print("Scanning: " + dir);
    print("--------------------------------------------------");

    for (f = 0; f < list.length; f++) {

        name = list[f];
        lower = toLowerCase(name);

        // Skip non-CSV, skip already-normalized files, skip averaged outputs
        if (!endsWith(lower, ".csv")) continue;
        if (indexOf(lower, "_norm") != -1) { print("Skipped (already normalized): " + name); continue; }
        if (!endsWith(lower, "_lineprofile.csv")) { print("Skipped (not a lineprofile file): " + name); continue; }

        path = dir + name;
        text = File.openAsString(path);
        if (lengthOf(text) == 0) { print("Skipped (empty): " + name); continue; }

        text = replace(text, "\r", "");
        lines = split(text, "\n");

        // --- Read header ---
        header = lines[0];

        // --- Read data rows ---
        nRows = 0;
        for (i = 1; i < lines.length; i++) {
            line = trim(lines[i]);
            if (line == "") continue;
            cols = split(line, ",");
            if (cols.length < 7) continue;
            nRows++;
        }

        if (nRows == 0) { print("Skipped (no valid rows): " + name); skipped++; continue; }

        // --- Allocate arrays ---
        imgCol    = newArray(nRows);
        xpos      = newArray(nRows);
        realDist  = newArray(nRows);
        ch1raw    = newArray(nRows);
        ch1bgsub  = newArray(nRows);
        ch2raw    = newArray(nRows);
        ch2bgsub  = newArray(nRows);

        // --- Parse rows ---
        row = 0;
        for (i = 1; i < lines.length; i++) {
            line = trim(lines[i]);
            if (line == "") continue;
            cols = split(line, ",");
            if (cols.length < 7) continue;

            // store image name as string index (just keep original)
            imgCol[row]   = 0; // placeholder — we rebuild from cols[0]
            xpos[row]     = parseFloat(cols[1]);
            realDist[row] = parseFloat(cols[2]);
            ch1raw[row]   = parseFloat(cols[3]);
            ch1bgsub[row] = parseFloat(cols[4]);
            ch2raw[row]   = parseFloat(cols[5]);
            ch2bgsub[row] = parseFloat(cols[6]);
            row++;
        }

        // --- Per-spindle min-max normalization of bgsub columns ---
        min1 = ch1bgsub[0]; max1 = ch1bgsub[0];
        min2 = ch2bgsub[0]; max2 = ch2bgsub[0];

        for (i = 1; i < row; i++) {
            if (ch1bgsub[i] < min1) min1 = ch1bgsub[i];
            if (ch1bgsub[i] > max1) max1 = ch1bgsub[i];
            if (ch2bgsub[i] < min2) min2 = ch2bgsub[i];
            if (ch2bgsub[i] > max2) max2 = ch2bgsub[i];
        }

        range1 = max1 - min1;
        range2 = max2 - min2;

        // --- Re-read image name from first data row for output ---
        firstLine = "";
        for (i = 1; i < lines.length; i++) {
            line = trim(lines[i]);
            if (line == "") continue;
            firstLine = line;
            break;
        }
        firstCols = split(firstLine, ",");
        imgName = firstCols[0];

        // --- Build output CSV (same 7-column structure) ---
        out = "Image,NormPos_um,RealDist_um,Ch1_raw,Ch1_bgsub,Ch2_raw,Ch2_bgsub\n";

        for (i = 0; i < row; i++) {
            if (range1 == 0) n1 = 0;
            else n1 = (ch1bgsub[i] - min1) / range1;

            if (range2 == 0) n2 = 0;
            else n2 = (ch2bgsub[i] - min2) / range2;

            out += imgName + "," +
                   xpos[i]     + "," +
                   realDist[i] + "," +
                   ch1raw[i]   + "," +
                   n1          + "," +
                   ch2raw[i]   + "," +
                   n2          + "\n";
        }

        // --- Save normalized file ---
        outName = replace(name, "_lineprofile.csv", "_lineprofile_norm.csv");
        File.saveString(out, dir + outName);
        print("Saved: " + outName);
        processed++;
    }

    print("--------------------------------------------------");
    print("DONE — " + processed + " files normalized, " + skipped + " skipped.");
    showMessage("Done", processed + " files normalized.\nNorm files saved alongside originals.");
}
