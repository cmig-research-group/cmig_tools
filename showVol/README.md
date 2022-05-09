# showVol
MR volume viewer and annotator. These instructions pertain to use with ABCD data, but showVol can visualize any 3d volumes. The ROI outlines will only be valid for ABCD data as these have been registered to the ABCD atlas.

# Setup

## Get code and prerendered atlas images

### 1) showVol code

If you don't have the code yet, get the latest code from github, currently:

```git clone git@github.com:cmig-research-group/cmig_tools.git```

To use showVol, you'll need to add `cmig_tools` and its sub-directories to your Matlab path.

### 2) atlas data files

You'll need prerendered data and anatomical ROIs stored in a showVolData directory. These are atlas-specific, with the current available being ABCD1_cor10 and ABCD2_cor10.

The files for showVolData can be downloaded from here: https://drive.google.com/drive/folders/1Uq04VEKRkgQ5SlWDfV9T4KXbwNuCY0nR?usp=sharing
Generally, you should download the most recently dated .tgz files. first download and expand the ```showVolData_yyyy-mm-dd.tgz```. If you are interested in plotting fiber orientation distributions (and why not!) also download the showVolData_FOD files and extract onto the ```showVolData``` directory which should then contain three subdirectories: Atlas, atlas_dspace, and FOD.

#### Path Configuration
Paths are saved in ```~/abcdConfig.json```
The first time you run showVol it will prompt you to edit the path for "showVolData" to point to the directory where you uncompressed it.

```
# ABCD Configuration: Set these paths to point to different data resorces necessary for FEMA (abcd-sync) or visualization tools (showVolData).
{
"data":
	{"abcd_sync":
		"/your/path/to/abcd-sync",
	"showVolData":
		"/your/path/to/showVolData"}    #<-------replace this with your showVolData directory path
}
```

# Usage

## Example code

```showVol/ABCD_anatomy_visualization_example.m``` gives a starting-point example you can run to test everything out. *Do make a copy* outside of the showVol directory so you can play with it without affecting the public git archive of showVol.

## Visualizing statistical maps from FEMA

```FEMA/DEMOS/FEMA_showVol_demo.m```

Has code to load and colorize your statistical maps, which can be adapted and added to your script.

# ABCD Atlas 

```showVol``` currently uses the ABCD2 atalas by default.

## ROI visualizations

For ABCD volumes, several atlases are available for visualizing ROIs: aparc, aseg, fiber, Pauli, thalamus
In addition, registered T1 and T2 and rendered rgb images for fiber orientation, FA, etc are available.

## FOD

Prerendered FOD (fiber orientation distribution).

For more information, see ```utils/loadPrerenderedIfNeeded.m```.

