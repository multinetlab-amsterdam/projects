## Spin Test ##

The scripts SpinTest.m, cortical2surface.m and SpinPermuFS.m aim to perform the spin test in volumetric data. 

The main goal of the spin test is to guarantee the statistical significance of the spatial distribution of cortical data.
Correlations between brain maps may in part be due to the higher chance that neighboring regions display comparable variations.
The spin test generates null models in order to check if the spatial overlap between two maps is significantly higher than any random rotation of the data in the cortical surface.
These scripts were developed based on the [spin test toolbox](https://github.com/spin-test/spin-test), expanding on, and adapting, their previous work on suface data to parcellated cortical data.
Contact person: @bern-maciel on github.

### Scripts ###
* SpinTest.m is the main script of this project and requires cortical2surface.m and SpinPermuFS.m. (Note: requires nearestneighbour.m from the [spin test toolbox](https://github.com/spin-test/spin-test))
* cortical2surface.m maps atlas-parecellated volumetric cortical data onto the fsaverage5 space (10242 vertices on the cortical surface) by inputing the volumetric value of each region onto the surface vertex correspondent to the geometric centroid. (Currently only supports cortical data parcellated according to the [AAL](https://www.gin.cnrs.fr/en/tools/aal/) and [BNA](https://atlas.brainnetome.org/) atlases)
* SpinPermuFS.m is an updated version of the script with the same name from the [spin test toolbox](https://github.com/spin-test/spin-test) to work smoothly with our input data.

### Notes on usage ###
* Please read carefully the headers of each function.
* Follow the file structure indicated in the scripts for a smoother execution:
  * ./atlas/ directory with annotation .annot and .label files for the BNA (available for download [here](https://atlas.brainnetome.org/download.html) on BN_Atlas_freesurfer.zip/fsaverage/label/) and AAL parcellations. These can be easily expanded upon by editing the first section of cortical2surface.m.
  * ./data/ directory with the .csv files of the input data with 2 columns represented the data to be tested parcellated under the same atlas.
* Requires previous instalation of FreeSurfer software and the environment variable FREESURFER_HOME.
* It was developed in a Unix environment, might not work smoothly for windows users.
