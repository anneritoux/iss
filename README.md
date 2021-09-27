The README of the master branch gives instructions of how to run the pipeline and view the results in 2D. Most of this still applies in 3D, but there are a few small differences listed below:

* You cannot save 4D tif files, so the tiles are saved differently. The file named as ```o.FileBase{r}_tT_cC.tif``` contains all the z plane images for round r, tile T, colour channel C.
* The figure number of ```o.plot``` is now 93454, so when calling ```iss_view_prob``` or ```iss_view_codes```, [the figure number will need to be changed to this](https://github.com/jduffield65/iss/blob/11c185d41534e46a4662a6089c50cb8fddb4b89f/bridge_process_template.m#L144-L145).
* The ```o.plot``` figure now has a scrollbar which allows you to see how the genes are distributed at each z coordinate. To concatenate genes from other planes onto a single plane, use [```o.PlotZThick```](https://github.com/jduffield65/iss/blob/11c185d41534e46a4662a6089c50cb8fddb4b89f/bridge_process_template.m#L134). This is defined such that all genes with z coordinate within ```+/-o.PlotZThick``` of the current z coordinate (chosen by the scrollbar) are shown. After changing this, ```iss_change_plot(o)``` needs to be run for it to take effect.
* Apart from these, there are of course many instances when 2D coordinates/shifts/transforms are converted into 3D.

# System requirements

- Matlab signal processing toolbox
- Matlab image processing toolbox