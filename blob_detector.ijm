var synap = "synaptopodin";
var homer = "Homer1c";
var mito = "mitochondrial";
var synap_blob_radius_label = "Approximate radius of " + synap + " 'blobs'";
var synap_thresh_method_label = "Threshold method for " + synap;
var synap_blob_radius = 3;
var synap_thresh_method = "Default";
var homer_blob_radius_label = "Approximate radius of " + homer + " 'blobs'";
var homer_thresh_method_label = "Threshold method for " + homer;
var homer_blob_radius = 2;
var homer_thresh_method = "Default";
var homer_chan = 3;
var synap_chan = 1;
var mito_chan = 2;
var search_radius = 3;
var homer_chan_label = homer + " channel";
var synap_chan_label = synap + " channel";
var mito_chan_label = mito + " channel";
var all_thresh_methods = getList("threshold.methods");
var homer_synap_outputs = "Homer1c-Synaptopodin_NearestNeighbours.csv";
var homer_label = "Homer1c Object Number";
var synap_label = "Synaptopodin Object Number";
var homer_mito_outputs = "Homer1c-Mitochondrial_IntensityMeasurements.csv";
var detected_synap_outputs = "Detected_Synaptopodin_Objects.ome.tiff";
var detected_homer_objects = "Detected_Homer1c_Objects.ome.tiff";
var search_radius_label = "Specify maximum distance from " + homer + " 'blobs' to search for " + mito + " signal"

file = File.openDialog("Select input file");

setBatchMode(true);

print("Opening " + file);

run("Bio-Formats Importer", "open=[" + file + "] autoscale color_mode=Composite rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT");
getDimensions(width, height, channels, slices, frames);

Dialog.create("Title");
Dialog.addSlider(homer_chan_label, 1, channels, homer_chan);
Dialog.addSlider(synap_chan_label, 1, channels, synap_chan);
Dialog.addSlider(mito_chan_label, 1, channels, mito_chan);
Dialog.addNumber(homer_blob_radius_label, homer_blob_radius);
Dialog.addNumber(search_radius_label, search_radius);
Dialog.addNumber(synap_blob_radius_label, synap_blob_radius);
Dialog.addChoice(homer_thresh_method_label, all_thresh_methods, homer_thresh_method);
Dialog.addChoice(synap_thresh_method_label, all_thresh_methods, synap_thresh_method);
Dialog.addDirectory("Choose location to save outputs", getDirectory("cwd"));
Dialog.show();

homer_chan = Dialog.getNumber();
synap_chan = Dialog.getNumber();
mito_chan = Dialog.getNumber();
homer_blob_radius = Dialog.getNumber();
search_radius = Dialog.getNumber();
synap_blob_radius = Dialog.getNumber();
homer_thresh_method = Dialog.getChoice();
synap_thresh_method = Dialog.getChoice();
output = Dialog.getString();

input = getTitle();

run("Split Channels");

synap_title = "C" + synap_chan + "-" + input;
mito_title = "C" + mito_chan + "-" + input;
homer_title = "C" + homer_chan + "-" + input;

images = getList("image.titles");

for (i = 0; i < images.length; i++) {
	chan_number = parseInt(substring(images[i], 1, 2));
	if((chan_number != synap_chan) && (chan_number != mito_chan) && (chan_number != homer_chan)){
		close(images[i]);
	}
}

print("Finding objects in synaptopodin channel...");

selectWindow(synap_title);
synap_labelled = segmentAndLabelBlobs(synap_blob_radius, synap_thresh_method);
close(synap_title);

print("Finding objects in Homer1c channel...");

selectWindow(homer_title);
homer_labelled = segmentAndLabelBlobs(synap_blob_radius, synap_thresh_method);
close(homer_title);

print("Calculating distances between Homer1c objects and nearest Synaptopodin object...");

run("3D Distances Closest", "image_a=" + File.getNameWithoutExtension(homer_labelled) + " image_b=" + File.getNameWithoutExtension(synap_labelled) + " number=1 distance=DistCenterCenterPix distance_maximum=100");
selectWindow("ClosestDistanceCCPix");
Table.renameColumn("LabelObj", homer_label);
Table.renameColumn("O1", synap_label);
Table.renameColumn("V1", "Distance");
saveAs("Results", output + File.separator() + homer_synap_outputs);
close(homer_synap_outputs);

print("Saving outputs in " + output);

selectWindow(synap_labelled);
run("Bio-Formats Exporter", "save=[" + output + File.separator() + detected_synap_outputs + "] compression=LZW");

selectWindow(homer_labelled);
run("Bio-Formats Exporter", "save=[" + output + File.separator() + detected_homer_objects + "] compression=LZW");

print("Quantifying mitochondrial intensity within Homer1c objects...");

measureMito(mito_title, homer_labelled, search_radius);

close("*");

windows = getList("window.titles");

for (i = 0; i < windows.length; i++) {
	close(windows[i]);
}

print("All done.");

setBatchMode(false);


function segmentAndLabelBlobs(blob_radius, thresh_method){
	// Step 1: Detect blobs
	run("FeatureJ Hessian", "largest smoothing=" + blob_radius);
	hessian = getTitle();

	// Step 2: Threshold result
	setAutoThreshold(thresh_method + " stack");
	run("Convert to Mask", "method=Default background=Light");
	mask = getTitle();

	// Step 3: Label objects
	run("Connected Components Labeling", "connectivity=6 type=[16 bits]");
	output = getTitle();
	close(hessian);
	close(mask);
	
	return output;
}

function measureMito(mito, labels, search_radius){
	print("Measuring distance 0...");
	getVoxelSize(width, height, depth, unit);
	anisotropy = depth / width;
	Table.create("Results");
	measureIntensities(mito, labels, 0, "Results", "Mean");
	last_label_image = labels;
	anistropyCount = 1;
	for (i = 0; i < search_radius; i++) {
		print("Measuring distance " + (i + 1) + "...");
		selectWindow(last_label_image);
		anistropyCount = dilateLabels(i, anisotropy, anistropyCount);
		dilated_labels = getTitle();
		imageCalculator("Difference create stack", last_label_image, dilated_labels);
		difference = getTitle();
		measureIntensities(mito, difference, (i + 1), "Results", "Mean");
		temp = last_label_image;
		last_label_image = dilated_labels;
		close(difference);
		close(temp);
	}
	selectWindow("Results");
	saveAs("Results", output + File.separator() + homer_mito_outputs);
}

function measureIntensities(mito, difference, i, resultsTable, resultsColumn){
	run("Intensity Measurements 2D/3D", "input=[" + mito + "] labels=[" + difference + "] mean");
	means = Table.getColumn(resultsColumn);
	selectWindow(resultsTable);
	Table.setColumn(resultsColumn + "_" + (i + 1), means);
}

function dilateLabels(iteration, anisotropy, anistropyCount){
	inputImage = getTitle();
	if(anistropyCount >= anisotropy){
		print("Performing 3D label dilation...");
		run("Label Morphological Filters", "operation=Dilation radius=1 from_any_label");
		anistropyCount = 1;
	} else {
		print("Performing 2D label dilation...");
		run("Duplicate...", "duplicate");
		imageTitle = getTitle();
		run("Stack to Images");
		images = getList("image.titles");
		for(i = 0; i < images.length; i++){
			selectWindow(images[i]);
			getDimensions(width, height, channels, slices, frames);
			if(slices < 2){
				run("Label Morphological Filters", "operation=Dilation radius=1 from_any_label");
				close(images[i]);
			}
		}
		run("Images to Stack", " ");
		rename(inputImage + " - Dilated");
		anistropyCount++;
	}
	return anistropyCount;
}
