// ============================================================
// Macro: Spindle pole-to-pole line profile extractor
// Adapted for: channel 2 (poles + kinetochores) & channel 3 (alpha-tubulin)
// ============================================================
// Description:
//   Extracts fluorescence intensity profiles along the pole-to-pole
//   axis of the mitotic spindle, reading from image CHANNELS 2 and 3:
//     Channel 2 = combined poles + kinetochores marker
//     Channel 3 = alpha-tubulin
//   Channel 1 (e.g. DNA/DAPI), if present, is ignored.
//
//   The user draws a line as wide as the spindle from pole to pole,
//   then draws a background rectangle. Background mean is subtracted
//   per channel. The profile is resampled to 161 bins spanning
//   -8 to +8 um, centered on the spindle midpoint.
//
// Usage:
//   1. Open a multi-channel sum projection image in Fiji/ImageJ.
//      Set line width to match spindle width BEFORE running
//      (Edit > Options > Line Width, or double-click the line tool).
//   2. Run via Plugins > Macros > Run.
//   3. Draw a straight line from pole to pole, then click OK.
//   4. Draw a rectangular background ROI outside the spindle, click OK.
//   5. Choose an output folder. A CSV is saved automatically.
//
// Input:
//   - Multi-channel sum projection image (.tif, .tiff, or .ome.tif)
//     with at least 3 channels
//   - Image must have pixel scale set in micrometers
//   - Channel 2: poles + kinetochores marker
//   - Channel 3: alpha-tubulin
//   (Channel 1 is ignored by this macro.)
//
// Output:
//   - <imagename>_lineprofile.csv with columns:
//     Image, NormPos_um, RealDist_um,
//     Ch1_raw, Ch1_bgsub, Ch2_raw, Ch2_bgsub
//
//   NOTE: to stay compatible with Macro_normalize_perSpindle.ijm and
//   the averaging macro (which expect a generic 7-column structure),
//   the CSV columns are still named Ch1/Ch2 — but the data itself
//   comes from image channels 2 (poles+kinetochores -> CSV "Ch1")
//   and 3 (alpha-tubulin -> CSV "Ch2"). Keep that mapping in mind
//   when interpreting downstream averaged output.
//
// Column definitions:
//   NormPos_um  — position along spindle axis centered at 0
//                 (range: -8 to +8 um, 161 bins, 0.1 um steps)
//   RealDist_um — actual distance from spindle center in um
//   Ch1_*       — poles + kinetochores marker (image channel 2)
//   Ch2_*       — alpha-tubulin (image channel 3)
//
// Notes:
//   - Use this output directly for statistical analysis
//   - Run Macro_normalize_perSpindle.ijm on this output to generate
//     per-spindle normalized files for visualization
//
// Dependencies: Fiji/ImageJ (https://fiji.sc)
//
// Authors: Kristjana Skendo, Meraldi Lab, University of Geneva
// Adapted from: https://github.com/kskendo/Spindle-LineProfile-Macros
// ============================================================

macro "Poles+KT / aTub spindle line profile (ch2 & ch3)" {

    if (nImages == 0) exit("Open an image first.");

    outDir = getDirectory("Choose folder to save CSV files");
    if (outDir == "") exit("No output folder selected.");

    title = getTitle();

    nNorm    = 161;
    normMin  = -8.0;
    normMax  =  8.0;
    normStep = (normMax - normMin) / (nNorm - 1);

    Stack.getDimensions(width, height, channels, slices, frames);
    if (channels < 3) exit("Need at least 3 channels (poles+KT = ch2, tubulin = ch3).");

    Stack.getPosition(c0, z0, t0);
    if (frames >= 1) Stack.setPosition(c0, z0, 1);

    getPixelSize(unit, pixelWidth, pixelHeight);
    if (unit == "pixel" || unit == "pixels")
        exit("Set image scale in micrometers first (Image > Properties).");

    run("Select None");
    waitForUser("Draw pole-to-pole line, then click OK.");

    if (selectionType() != 5)
        exit("Draw a straight line.");

    getSelectionCoordinates(xp, yp);
    x1 = xp[0]; y1 = yp[0];
    x2 = xp[1]; y2 = yp[1];

    dx = x2 - x1;
    dy = y2 - y1;
    lineLengthPx = sqrt(dx*dx + dy*dy);
    lineLengthUm = lineLengthPx * pixelWidth;

    run("Select None");
    waitForUser("Draw background rectangle, then click OK.");

    if (selectionType() != 0)
        exit("Draw rectangular background ROI.");

    // Channel 2 = poles + kinetochores marker
    Stack.setPosition(2, z0, 1);
    getStatistics(area1, bg1);

    // Channel 3 = alpha-tubulin
    Stack.setPosition(3, z0, 1);
    getStatistics(area2, bg2);

    makeLine(x1, y1, x2, y2);

    Stack.setPosition(2, z0, 1);
    prof1 = getProfile();   // poles + kinetochores

    Stack.setPosition(3, z0, 1);
    prof2 = getProfile();   // alpha-tubulin

    nRaw = prof1.length;

    fileBase = replace(title, ".ome", "");
    fileBase = replace(fileBase, ".tif", "");
    fileBase = replace(fileBase, ".tiff", "");

    out = "Image,NormPos_um,RealDist_um,Ch1_raw,Ch1_bgsub,Ch2_raw,Ch2_bgsub\n";

    for (j = 0; j < nNorm; j++) {

        normPos = normMin + j * normStep;
        frac    = (normPos - normMin) / (normMax - normMin);
        rawPos  = frac * (nRaw - 1);

        i0 = floor(rawPos);
        i1 = floor(rawPos + 1);

        if (i0 < 0) i0 = 0;
        if (i1 >= nRaw) i1 = nRaw - 1;

        if (i0 == i1) {
            val1 = prof1[i0];
            val2 = prof2[i0];
        } else {
            w    = rawPos - i0;
            val1 = prof1[i0]*(1-w) + prof1[i1]*w;
            val2 = prof2[i0]*(1-w) + prof2[i1]*w;
        }

        corr1    = val1 - bg1;
        corr2    = val2 - bg2;
        realDist = frac * lineLengthUm - lineLengthUm/2.0;

        out += fileBase + "," + normPos + "," + realDist + "," +
               val1 + "," + corr1 + "," + val2 + "," + corr2 + "\n";
    }

    File.saveString(out, outDir + fileBase + "_lineprofile.csv");
    print("Saved: " + fileBase + "_lineprofile.csv");
    showMessage("Done", "Profile saved to:\n" + outDir + fileBase + "_lineprofile.csv");
}
