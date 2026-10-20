// If an apporpriate line is drawn on an embryo, this macro will 
// straighten the line and measure the average intensity at each 
// AP position along that line. 
// Output: the straightened image (for validation) 
// and a .csv of  the intensity values. 
// X0: AP position. Y0-n: genes. 

// #@ String (label="Please a suffix", description="Name field") suffix

suffix = '_'

//change me! 
// path = "/home/shannon/Dropbox/DPhil/2021/image-analysis/tbox_gene_expression/1D/data/"
// path = "C:/Users/shil5659/Dropbox/DPhil/2021/image-analysis/PSM_patterning/segmentation_onset/data/"

// path = getDir("Choose a Directory");

path = "C:/Users/admin/Dropbox/DPhil/2021/image-analysis/PSM_patterning/tbox_10ss/data/"

linewidth =50

title = getTitle();
getPixelSize(unit, pixelWidth, pixelHeight);  

// Save ROI 
roiManager("Save", path + title + suffix + ".roi");
// Save original image 
saveAs("tiff", path + title + suffix);

//straighten only works on greyscale
Stack.setDisplayMode("grayscale");
run("Straighten...", "title="+title+" line=" + linewidth + " process");
//run("Gaussian Blur...", "sigma=1 scaled stack");//run after straighten!
save(path + title + suffix + "_line.tif"); 

//make the whole image a ROI and add to ROI manager
getDimensions(width, height, channels, slices, frames);
print(channels, slices); 
for (i=1; i<=channels; i++){
	setSlice(i); 
	makeRectangle(0, 0, width, height);
	roiManager("Add");
}
roiManager("Multi Plot");
Plot.showValues(); 
saveAs("Results", path + title + suffix + ".csv");
//close results
close("Results");  
close("Roi Manager");  
close("Profiles");
//close zoom and main image
//close(title); 
rename(title);