//The ImageJ macro for nuclear envelope wrinkling index quantification was originally developed by Cosgrove et al.: Cosgrove BD, Loebel C, Driscoll TP, //Tsinman TK, Dai EN, Heo SJ, et al. Nuclear envelope wrinkling predicts mesenchymal progenitor cell mechano-response in 2D and 3D microenvironments. //Biomaterials. 2021;270: 120662. doi:10.1016/J.BIOMATERIALS.2021.120662
//We made minor changes to their macro to adapt it to our data.

// ask user to select the folder containing cropped max-projected Lamin A/C staining images
dir = getDirectory("Select A folder");
imagecount=getFileList(dir);

regFactor = 6;
//run("Z Project...", "projection=[Max Intensity]");
//title=getTitle;
//LaminAC=title;
wrinkleIndex=newArray();
title=newArray();

count=0;
for (i = 0; i < lengthOf(imagecount); i++) {
	open(dir+"/"+imagecount[i]);
	title[count]=getTitle();
	rename("image");
	run("Enhance Contrast", "saturated=0.35");
	//setMinAndMax(0, 67);
	//setMinAndMax(7, 62);
	run("Duplicate...", " ");
	run("Threshold...");
	setAutoThreshold("Default dark");
	setOption("BlackBackground", true);
	run("Convert to Mask");
	selectWindow("Threshold");
	run("Close");
	
	run("Images to Stack");
	Stack.setPosition(1,2,1);
	setTool("wand");
	waitForUser("Click nuclei center, make sure nuclei is outlined or choose different point , then press ok");
	run("Scale... ", "x=0.9 y=0.9 centered");
	Stack.setPosition(1,1,1);
	// Get max/min intensity, mean and median pixel values, std dev of intensity, and the nuclear aspect ratio. 
	run("Set Measurements...", "area mean min median standard redirect=None decimal=3");
	run("Measure");
	maxInt=getResult('Max',nResults-1);
	minInt=getResult('Min',nResults-1);
	mean=getResult('Mean',nResults-1);
	median=getResult('Median',nResults-1);
	stddev=getResult('StdDev',nResults-1);
	AR=getResult('Aspect_Ratio',nResults-1);
	IJ.deleteRows(0, nResults-1);
	getMinAndMax(min, max);
	run("Stack to Images");
	
	LaminAC="image";
	selectWindow(LaminAC);
	maxInt=((mean)-(1.5*stddev))/regFactor; 
	// Run sobdel-based edge detection 
	label=("output");
	run("FeatureJ Edges", "compute smoothing=0.95 suppress lower="+maxInt+" higher="+maxInt);
	//saveAs("PNG",dir+label+"_Cell"+i+"_rF"+regFactor+"_Edges");
	selectWindow(LaminAC);
	//selectWindow("MAX_1kpaS86 1");
	close();
	run("Images to Stack");
	
	// Select boundaries for binarized edge map and shrink slightly to remove the contrast outline of the nucleus. 
	setTool("wand");
	waitForUser("Click nuclei center, make sure nuclei is outlined or choose different point , then press ok");
	run("Scale... ", "x=0.96 y=0.96 centered");
	Stack.setPosition(1,2,1)
	// Calculate edge detection statistics 
	run("Measure");
	medianVal=median;
	stdVal=stddev;
	minVal=minInt;
	threshold=maxInt;
	meanVal=mean;
	
	wrinkleIndex[count]=getResult('Mean',nResults-1)/255*100; ; // since binarized, divided by avg if every pixel was lit (/255) in percent (*100) 
	//nuclearArea=getResult('Area');
	//run("Clear Results");
	//setResult("Wrinkling index",nResults-1,wrinkleIndex);
	close("Stack");
	count=count+1;
}
Array.show("Results",title,wrinkleIndex)
saveAs("Results",dir+"outputWI.csv");

// Output and save data 
// Array.show("Results",wrinkleIndex,nuclearArea,threshold,minVal,medianVal,meanVal,stdVal);
 //saveAs("Results",dir+output);
