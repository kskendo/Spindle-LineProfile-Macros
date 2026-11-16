// ============================================================
// Macro: Spindle pole-to-pole line profile extractor
// HAUS6 / HURP / Poles+Kinetochores - interactive channel assignment
// ============================================================
// Description:
//   Extracts pole-to-pole line profiles for THREE channels at once
//   (HAUS6, HURP, and the poles+kinetochores marker), asking you
//   to confirm which physical image channel is which for every
//   image - same reasoning as the tubulin pipeline's interactive
//   macro, extended to a third sampled channel.
//
//   Regardless of which physical channel you assign to each role,
//   the OUTPUT CSV always uses this fixed convention:
//     Ch1 = HAUS6
//     Ch2 = HURP
//     Ch3 = Poles + kinetochores marker
//
// Usage:
//   1. Open a 3-channel sum projection image in Fiji/ImageJ.
//      Set line width to match spindle width BEFORE running
//      (Edit > Options > Line Width, or double-click the line tool).
//   2. Run via Plugins > Macros > Run.
//   3. A dialog asks you to assign each role to a channel number -
//      confirm/select correctly for THIS image.
//   4. Draw a straight line from pole to pole, click OK.
//   5. Draw a rectangular background ROI outside the spindle, click OK.
//   6. Choose an output folder. A CSV is saved automatically.
//
// Input:
//   - 3-channel sum projection image (.tif, .tiff, or .ome.tif)
//   - Image must have pixel scale set in micrometers
//
// Output:
//   - <imagename>_lineprofile.csv, 9-column structure:
//     Image, NormPos_um, RealDist_um,
//     Ch1_raw, Ch1_bgsub, Ch2_raw, Ch2_bgsub, Ch3_raw, Ch3_bgsub
//     (Ch1 = HAUS6, Ch2 = HURP, Ch3 = Poles+kinetochores,
//     regardless of which physical channel number you assigned to
//     each role)
//
// Notes:
//   - This is a DIFFERENT CSV structure (9 columns, 3 channels)
//     from the earlier tubulin/PRC1 pipeline (7 columns, 2
//     channels) - use the matching per-cell summary macro for this
//     dataset, not the older one.
//   - The channel assignment you make is printed to the Log
//     window for traceability.
//
// Dependencies: Fiji/ImageJ (https://fiji.sc)
//
// Authors: Kristjana Skendo, Meraldi Lab, University of Geneva
// ============================================================

macro "Spindle line profile - HAUS6/HURP/PolesKT (interactive)" {

    if (nImages == 0) exit("Open an image first.");

    title = getTitle();

    Stack.getDimensions(width, height, channels, slices, frames);
    if (channels < 3) exit("Need at least 3 channels (HAUS6, HURP, Poles+kinetochores).");

    channelList = newArray(channels);
    for (c = 1; c <= channels; c++) channelList[c-1] = "Channel " + c;

    Dialog.create("Assign channels - " + title);
    Dialog.addMessage("Confirm which channel is which for THIS image:");
    Dialog.addChoice("HAUS6 is:", channelList, channelList[0]);
    Dialog.addChoice("HURP is:", channelList, channelList[minOf(1, channels-1)]);
    Dialog.addChoice("Poles + kinetochores marker is:", channelList, channelList[minOf(2, channels-1)]);
    Dialog.show();

    haus6Choice = Dialog.getChoice();
    hurpChoice  = Dialog.getChoice();
    ktChoice    = Dialog.getChoice();

    haus6Ch = parseInt(substring(haus6Choice, 8));
    hurpCh  = parseInt(substring(hurpChoice, 8));
    ktCh    = parseInt(substring(ktChoice, 8));

    if (haus6Ch == hurpCh || haus6Ch == ktCh || hurpCh == ktCh)
        exit("All three roles must be assigned to different channels. Please re-run and pick distinct channels.");

    print("--------------------------------------------------");
    print("Image: " + title);
    print("Assigned: HAUS6 = channel " + haus6Ch + ", HURP = channel " + hurpCh +
          ", Poles+kinetochores = channel " + ktCh);

    outDir = getDirectory("Choose folder to save CSV files");
    if (outDir == "") exit("No output folder selected.");

    nNorm    = 161;
    normMin  = -8.0;
    normMax  =  8.0;
    normStep = (normMax - normMin) / (nNorm - 1);

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

    // background, for all three channels
    Stack.setPosition(haus6Ch, z0, 1);
    getStatistics(area1, bg1);

    Stack.setPosition(hurpCh, z0, 1);
    getStatistics(area2, bg2);

    Stack.setPosition(ktCh, z0, 1);
    getStatistics(area3, bg3);

    makeLine(x1, y1, x2, y2);

    Stack.setPosition(haus6Ch, z0, 1);
    prof1 = getProfile();   // HAUS6

    Stack.setPosition(hurpCh, z0, 1);
    prof2 = getProfile();   // HURP

    Stack.setPosition(ktCh, z0, 1);
    prof3 = getProfile();   // Poles + kinetochores

    nRaw = prof1.length;

    fileBase = replace(title, ".ome", "");
    fileBase = replace(fileBase, ".tif", "");
    fileBase = replace(fileBase, ".tiff", "");

    out = "Image,NormPos_um,RealDist_um,Ch1_raw,Ch1_bgsub,Ch2_raw,Ch2_bgsub,Ch3_raw,Ch3_bgsub\n";

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
            val3 = prof3[i0];
        } else {
            w    = rawPos - i0;
            val1 = prof1[i0]*(1-w) + prof1[i1]*w;
            val2 = prof2[i0]*(1-w) + prof2[i1]*w;
            val3 = prof3[i0]*(1-w) + prof3[i1]*w;
        }

        corr1    = val1 - bg1;
        corr2    = val2 - bg2;
        corr3    = val3 - bg3;
        realDist = frac * lineLengthUm - lineLengthUm/2.0;

        out += fileBase + "," + normPos + "," + realDist + "," +
               val1 + "," + corr1 + "," + val2 + "," + corr2 + "," + val3 + "," + corr3 + "\n";
    }

    File.saveString(out, outDir + fileBase + "_lineprofile.csv");
    print("Saved: " + fileBase + "_lineprofile.csv");
    showMessage("Done", "Profile saved to:\n" + outDir + fileBase + "_lineprofile.csv" +
                 "\n\nChannel assignment used:\nHAUS6 = ch" + haus6Ch + ", HURP = ch" + hurpCh +
                 ", Poles+KT = ch" + ktCh);
}
