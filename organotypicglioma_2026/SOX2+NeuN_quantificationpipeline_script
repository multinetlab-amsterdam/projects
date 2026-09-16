//QuPath Quantification Pipeline for SOX2 and NeuN IHC staining
//NOTE: this script documents the steps taken for analysis in QuPath. 
//The script is NOT intended to be re-run directly on new images without adjustment.
//The dataset showed staining/intensity heterogeneity across images, so threshold values were set manually per image, as explained per step below.

//STEP 1: Set image type
setImageType('FLUORESCENCE');

//STEP 2: Manual delineation of the tissue section (to exclude artefacts)
//This step is performed manually in the QuPath viewer, so no corresponding code is shown.

//STEP 3: Identify tissue regions (to exclude damaged areas)
//A pixel classifier ("tissueidentification") was used to detect tissue.
//Its threshold was determined by calculating average pixel intensity across a random selection of positive areas (all channels), followed by visual inspection to confirm accurate tissue detection.
//If detection was not accurate, the threshold was adjusted and re-inspected until tissue was correctly identified.
createAnnotationsFromPixelClassifier("tissueidentification", 0.0, 0.0, "INCLUDE_IGNORED")

//STEP 4: DAPI-based cell detection
//Detection parameters (including the intensity threshold) were repeatedly tested to find optimal detection settings.
//Results were visually inspected using at least 10 widely distributed cell clusters across the image to confirm accuracy.
//For all images, an intensity threshold of 120 was used, as this appeared accurate across all images.
runPlugin('qupath.imagej.detect.cells.WatershedCellDetection', '{"detectionImage":"DAPI","requestedPixelSizeMicrons":0.3,"backgroundRadiusMicrons":18.0,"backgroundByReconstruction":true,"medianRadiusMicrons":0.0,"sigmaMicrons":1.7,"minAreaMicrons":15.0,"maxAreaMicrons":100.0,"threshold":120.0,"watershedPostProcess":true,"cellExpansionMicrons":3.0,"includeNuclei":true,"smoothBoundaries":true,"makeMeasurements":true}')

//STEP 5: Classify cells for SOX2 and NeuN positivity
//The object classifier ("SOX2+NeuN_1") combines two single-measurement classifiers, one per marker, each based on:
//   - object filter: cells
//   - channel filter: marker of interest (SOX2 or NeuN)
//   - measurement: nucleus channel standard deviation
//Each single-measurement threshold was determined by calculating the average of that measurement across a random selection of at least 6 positive cells, followed by visual inspection.
//If the resulting classification was not accurate, the threshold was adjusted and re-inspected until classification appeared correct.
runObjectClassifier("SOX2+NeuN_1");

//STEP 6: Apply pipeline to remaining images in the dataset
//For each subsequent image, the same parameters were tried first for tissue detection (Step 3), cell detection (Step 4), and classification (Step 5).
//Only intensity-related thresholds for tissue detection (Step 3) and classification (Step 5) were adjusted where necessary, and only to the minimal extent needed for accurate results.
