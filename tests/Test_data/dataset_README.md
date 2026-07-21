# Example dataset — QCI Lab logo (two-photon polymerization)

Five in-line (Gabor) holograms of a test target fabricated by two-photon
polymerization, recorded with a lensless holographic microscope. The dataset is
intended as a working example for [QCI HoLLoApp](https://github.com/QCI-LAB/QCI_HoLLoApp):
it exercises the full processing chain — background normalization, autofocusing,
shift correction, and iterative phase retrieval.

---

## Sample

The imaged structure is a custom-made **QCI Lab logo**, fabricated with **two-photon polymerization direct laser writing (Photonic Professional GT2 from Nanoscribe GmbH)**, with IP-S polymeric photoresist on a glass substrate. The structure is a predominantly phase object: it is nearly transparent at the illumination wavelength and produces weak amplitude contrast, which makes it a representative test case for phase-retrieval reconstruction rather than simple intensity imaging.
 

---

## Acquisition

| Parameter | Value |
|---|---|
| Configuration | Lensless in-line (Gabor) holography |
| Camera | Allied Vision Alvium 1800 U-2050m (monochrome) |
| Pixel size | 2.4 µm |
| Illumination wavelength | 0.561 µm |
| Number of holograms | 5 |


### Files

| File | Sample–sensor distance *z* [µm] | Notes |
|---|---|---|
| `Photo1.mat` | 3148 | |
| `Photo2.mat` | 4162 | |
| `Photo3.mat` | 5186 | |
| `Photo4.mat` | 6218 | |
| `Photo5.mat` | 7250 | |

The holograms were recorded at different sample–sensor distances (multi-height
acquisition), which provides the diversity required by the iterative reconstruction.
Distances are nominal values read from the translation stage; the app refines them
during autofocusing.

---

## Data format

Files are stored as raw camera output, without background subtraction, flat-field
correction, or normalization. Normalization is applied by the application at load
time, so the files can also be used to verify that step independently.

---

## Using the dataset with QCI HoLLoApp

1. **Load Data → Load Holograms**, and select this folder.
2. In the **HoloPrep** window, enable the high-pass filter (background division) and
   add all five holograms, then click **Load**.
3. Set the acquisition parameters in the *Default values* panel:
   - Pixel Size: `2.4` µm
   - Wavelength: `0.561` µm
   - Height: the nominal *z* of the first hologram
4. Set **Mode** to *Multi-height* and enter the per-hologram distances in the table.
5. In the **Focus** tab, refine each distance manually or with **Auto focus**
   (dark-focus metric). Enable *Pass forward* to propagate a corrected distance to the
   remaining holograms.
6. In the **Shift** tab, run **Auto correct** in *Shift and Scaling* mode to compensate
   for lateral drift and magnification change between heights.
7. In the **Reconstruction** tab, run the reconstruction. Suggested starting point:
   σ = `2`, iterations = `5`. Inspect the result in *Phase* mode — the logo is a phase
   object and will be most visible there.

Expected result: a legible QCI Lab logo in the reconstructed phase, free of the twin-image
artefact visible in a single back-propagated hologram.

---

## Reconstruction method

The application implements Iterative Gabor Averaging (IGA), which combines
Gerchberg–Saxton iteration with Gabor averaging. The σ parameter balances twin-image
suppression against shot-noise suppression: lower σ removes the twin image more
effectively, higher σ suppresses noise. See:

> M. Rogalski, P. Arcab, E. Wdowiak, J. Á. Picazo-Bueno, V. Micó, M. Józwik,
> M. Trusiak, *Hybrid iterating-averaging low photon budget Gabor holographic
> microscopy* (2024).

---

## License and citation

This dataset is licensed under Creative Commons Attribution 4.0 International (CC BY 4.0). See LICENSE.md in this folder for the full text.

In short: you are free to share and adapt the data for any purpose, including commercially, as long as you give appropriate credit, indicate if changes were made, and link back to the license. See https://creativecommons.org/licenses/by/4.0/ for the human-readable summary.

If you use this dataset, please cite it as:

 Wdowiak, E., Rogalski, M. & Trusiak, M. (2026). QCI lab Two-photon polymerization logo [Data set]. Zenodo. https://doi.org/10.5281/zenodo.21473964

If you also use the processing software, please cite it separately — see CITATION.cff in the QCI HoLLoApp repository.

---

## Contact

Quantitative and Computational Imaging Laboratory (QCI Lab)
Institute of Micromechanics and Photonics, Warsaw University of Technology

> contact emails: emilia.wdowiak.dokt@pw.edu.pl, maciej.trusiak@pw.edu.pl.
