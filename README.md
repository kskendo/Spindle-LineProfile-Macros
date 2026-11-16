# Spindle-LineProfile-Macros

**Repository:** https://github.com/kskendo/Spindle-LineProfile-Macros

Fiji/ImageJ macros for pole-to-pole fluorescence intensity line profiles of the mitotic spindle, measured on sum projection images, and for per-cell parameters derived from those profiles.

Developed for:
> Skendo K, Harry D, Castrogiovanni C, Meraldi P. Spindle assembly factors form and stabilize bridging microtubules. *Cell Reports* (in revision).

---

## Contents

| Folder | Used for | Pipeline |
|--------|----------|----------|
| `v1_general_two_channel/` | General two-channel pipeline; used for α-tubulin and HAUS6 profiles after taxol (Figures 4N, 4O) | Two-channel extraction → per-spindle normalization → averaging by condition |
| `v2_folder_pipeline/` | siCTRL vs siTPX2 (Figure 4), α-tubulin after all single and double depletions (Figure S5), TPX2/HAUS6/HURP/NuSAP/Kid after importazole (Figure 6) | Extraction with channel assignment → RAW + normalized averaging → per-cell parameters |

Both versions use the same principle: a line is drawn from pole to pole, a square background ROI is drawn outside the spindle, the background mean is subtracted per channel, and the profile is resampled to 161 points from pole to pole, centred at the spindle midpoint.

---

## Requirements

- [Fiji/ImageJ](https://fiji.sc)
- Multi-channel sum projection images (`.tif`, `.tiff` or `.ome.tif`) with the pixel size set in µm (Image > Properties)

---

## What changed from v1 to v2

The measurement is identical (pole-to-pole line, background square, background subtraction, 161 points from −8 to +8, same CSV columns). v2 differs in:

| | v1 | v2 |
|---|---|---|
| Channels read | Image channels 1 and 2 | Channels 2 and 3 (poles + kinetochores, α-tubulin), or three channels assigned in a dialog (HAUS6, HURP, poles + kinetochores) |
| Averaging | Normalize first, then average RAW and normalized files in two runs | RAW and normalized averages in one run |
| Per-cell parameters | — | Spindle length, mean intensity, AUC, tubulin dip depth, kinetochore peak prominence |

## v2 — folder pipeline (TPX2, importazole, depletion panel)

### Step 1 — Extract one profile per cell (`1_extract/`)

| Macro | Channels |
|-------|----------|
| `macro_lineprofile_polesKT_tubulin_ch2ch3.ijm` | Image channel 2 = poles + kinetochores, channel 3 = α-tubulin (channel 1, e.g. DAPI, ignored) |
| `macro_lineprofile_HAUS6_HURP_PolesKT_interactive.ijm` | Three channels assigned in a dialog (HAUS6, HURP, poles + kinetochores); the assignment is printed to the Log |

1. Open a sum projection. Set the line width to 100 pixels (double-click the line tool) and keep it the same for all cells.
2. Run the macro (**Plugins > Macros > Run**).
3. Draw a straight line from pole to pole, then OK.
4. Draw a rectangular background ROI outside the spindle, then OK.
5. Choose the output folder. `<image>_lineprofile.csv` is saved.

**Output columns:** `Image, NormPos_um, RealDist_um, Ch1_raw, Ch1_bgsub, Ch2_raw, Ch2_bgsub[, Ch3_raw, Ch3_bgsub]`

- `NormPos_um`: position on the normalized pole-to-pole axis, −8 to +8 (161 points, 0 = spindle midpoint). The two poles are always at −8 and +8, whatever the spindle length.
- `RealDist_um`: true distance from the spindle midpoint in µm. Spindle length = last − first value.
- `*_bgsub`: background-subtracted intensities. Use these for all quantifications.

### Step 2 — Average profiles per condition (`2_average/`)

| Macro | Dataset |
|-------|---------|
| `Combined_RAW_NORM_averaging_siCTRL_siTPX2.ijm` | siCTRL vs siTPX2, two channels; RAW and normalized averages in one run |
| `Combined_RAW_NORM_averaging_HAUS6_HURP_PolesKT.ijm` | DMSO vs importazole, three channels |
| `Average_RAW_lineprofiles_siCTRL_siTPX2.ijm`, `Average_NORM_lineprofiles_siCTRL_siTPX2.ijm` | Same as the combined macro, run separately on RAW or on normalized files |

1. Put all RAW `_lineprofile.csv` files of one experiment (all replicates, all conditions) in one folder. Conditions are read from the file names.
2. Run the macro and select the folder.

**Output:** `averaged_profiles_RAW_*.csv` and `averaged_profiles_NORM_*.csv` — mean, SEM and n per position and condition, ready to paste into a Prism XY table (Mean, SEM, N). Normalization is min–max per cell (0–1).

Position-by-position averages are for plotting only. The 161 points of one cell are not independent, so statistics are done on one value per cell or per replicate (Step 3).

### Step 3 — Per-cell parameters (`3_per_cell/PerCell_FullSummary.ijm`)

One row per cell:

| Parameter | Calculated from | Definition |
|-----------|-----------------|------------|
| `SpindleLength_um` | `RealDist_um` | Pole-to-pole distance |
| `MeanCh1`, `MeanCh2` | RAW profile | Mean background-subtracted intensity of the profile |
| `AUC_Ch1`, `AUC_Ch2` | RAW profile over `RealDist_um` | Area under the profile (total amount) |
| `TubulinDipDepth` | Normalized tubulin profile | Mean of the flanking windows (−5 to −3 and +3 to +5) minus mean of the central window (−1 to +1) |
| `KTPeakProminence` | Normalized poles/kinetochores profile | Mean of the central window (−1 to +1) minus mean of the edge windows (−8 to −6 and +6 to +8) |

Per-cell values are then averaged per replicate, and replicate means are used for statistics.

---

## v1 — general two-channel pipeline (used for taxol)

| Step | Macro | What it does |
|------|-------|--------------|
| 1 | `macro_lineprofile-new_2channels.ijm` | Extracts background-subtracted profiles for two channels (line as wide as the spindle, background square behind the poles). One CSV per image. |
| 2 | `Macro_normalize_perSpindle.ijm` | Min–max normalizes each profile (0–1) per spindle and channel. |
| 3 | `Average_line_profiles_DMSO_CTRL_15NM_one_folder.ijm` | Averages a folder of CSVs by condition (file names containing `DMSO`, `CTRL` or `15NM`); run once on RAW and once on normalized files. |

---

## Citation

If you use these macros, please cite:

> Skendo K et al. Spindle assembly factors form and stabilize bridging microtubules. *Cell Reports*. DOI to be added upon publication.

## License

MIT License — see `LICENSE`.

## Contact

Kristjana Skendo, Meraldi Lab, Department of Cell Physiology and Metabolism, University of Geneva
Lead contact: patrick.meraldi@unige.ch
