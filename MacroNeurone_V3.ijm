run("Close All");
roiManager("Reset");
run("Clear Results");

//////////////////////// CHANGER LES PARAMETRES EN FONCTION DE L ACQUISITION ////////////////////////
TailleDots=5;
Lentille=1;
pixelsize=0.16/Lentille;
ChannelPattern=4;
ChannelvGLUT1=1;
ChannelPSD=3;
ChannelMAP2=2;
ProtOrder=newArray("Cy5","Cy3","GFP","DAPI","BF");
size=((TailleDots*2)/(pixelsize));
/////////////////////////////////////////////////////////////////////////////////////////////////////


//Ouverture des fichier nd + enregistrement en fichier tiff 
input=getDirectory("Sélectionner le répertoire d’entrée") ;
output_All=getDirectory("Selectionner repertoire de sortie") ;
files=getFileList(input) ;

setBatchMode(false);
for(i=0 ; i<files.length ; i++){
	run("Close All");
	run("Clear Results");
	if(endsWith(files[i], ".nd")){
	image=files[i];
	path = input + image;
    run("Bio-Formats Macro Extensions");
    Ext.setId(path);
    Ext.getCurrentFile(file);
    Ext.getSeriesCount(seriesCount);
    run("Bio-Formats Importer", "open=&path autoscale color_mode=Default view=Hyperstack stack_order=XYCZT series_");
    output=output_All+File.separator+image;
   	File.makeDirectory(output);
	saveAs("Tiff", output + File.separator + image) ;
	imageTiff=getTitle();
	ImageName_woExt=File.getNameWithoutExtension(imageTiff);
	ROImakersMassCenter(size,imageTiff,ChannelPattern,ImageName_woExt,output);
	Detection(imageTiff,ChannelMAP2,ImageName_woExt,output);
	mesurenonspec=0; 
	mesurenonspec=StackMakers(imageTiff,ProtOrder,ImageName_woExt,mesurenonspec,output);
	
	if (mesurenonspec==0) {MeasureNonSpe (ImageName_woExt,ProtOrder,ChannelPattern,ChannelvGLUT1,pixelsize,output);}
	Measure(ImageName_woExt,ProtOrder,ChannelPattern,ChannelvGLUT1,pixelsize,output);
	run("Clear Results");

	if (mesurenonspec==0) {MeasureNonSpe (ImageName_woExt,ProtOrder,ChannelPattern,ChannelPSD,pixelsize,output);}
	Measure(ImageName_woExt,ProtOrder,ChannelPattern,ChannelPSD,pixelsize,output);
	run("Clear Results");
	
	output2=output+File.separator+"Multisynapse_Measure";
   	File.makeDirectory(output2);
	Multisynapse(ImageName_woExt,ProtOrder,ChannelvGLUT1,pixelsize,output2,output);
	run("Clear Results");
	
	Multisynapse(ImageName_woExt,ProtOrder,ChannelPSD,pixelsize,output2,output);
	run("Clear Results");
	
}}
setBatchMode("exit and display");


///////////////////////////////////////         Fonction pour déterminer les ROI avec le centre de masse       ////////////////////////
function ROImakersMassCenter(size,imageTiff,ChannelPattern,ImageName_woExt,output) { 
roiManager("Reset");
selectWindow(imageTiff);
run("Duplicate...", "title=Pattern duplicate channels="+ChannelPattern);
run("Set Measurements...", "center display redirect=None decimal=3");
setAutoThreshold("Default dark");
makeRectangle(size/2.5, size/2.5, getWidth()-size/2.5, getHeight()-size/2.5);
wait(50);
run("Analyze Particles...", "size=50-Infinity display exclude add slice");
nombreDots=roiManager("Count");
roiManager("reset");


for (i = 1; i <nombreDots-1; i++) {
	xpoints=getResult("XM", i);
	ypoints=getResult("YM", i);
	setTool("rectangle");
	makeRectangle(xpoints-size/2, ypoints-size/2, size, size);
	roiManager("Add");
}
run("Clear Results");
roiManager("save", output+File.separator+"ROI_all_"+ImageName_woExt+".zip");
Selectclose("Pattern");
}

//////////////////////////////////    dETECTE des ROI où les zones sans marquage OU CORPS cellulaire    //////////


// Supprimer les zones sans marquage de MAP2 et vgLUT1 
// Combinaison  des deux channeles pour etre sur  d'avoir tout les points 
function Detection(imageTiff,ChannelMAP2,ImageName_woExt,output) { 
	selectWindow(imageTiff);
	run("Select None");
	run("Duplicate...", "title=VGLUT1 duplicate channels="+ChannelvGLUT1);
	run("Duplicate...", "title=MAP2 duplicate channels="+ChannelMAP2);
	imageCalculator("Add create", "VGLUT1","MAP2");
	selectWindow("Result of VGLUT1");
	resetThreshold();
	setAutoThreshold("Default dark ");
	run("Convert to Mask");
	setOption("BlackBackground", true);
	run("Erode");
	run("Erode");
	run("Dilate");
	run("Dilate");
	selectWindow("Result of VGLUT1");
	
// Selection des zones sans Neurones et élimination des zones correspondant au corps cellulaires  	
	for(j=0; j<roiManager("Count"); j++){
	run("Select None");
	roiManager("Select", j);
	
	if(getValue("Max")<255){
		roiManager("Rename", "PatternNonspe"+j);
		roiManager("Set Color", "blue");
	}else if(getValue("Mean")>175) {
		roiManager("Rename", "CellBody"+j);
		roiManager("Set Color", "red");
		}}	 

roiManager("save", output+File.separator+"ROI_All_"+ImageName_woExt+".zip");
jmax = roiManager("Count");
	for(j=0; j<jmax; j++){
	selectWindow(imageTiff);
	roiManager("Select", j);
	ROIcolor=Roi.getStrokeColor;
	if(ROIcolor=="red"){
		roiManager("Delete"); 
        j  = j - 1;
        jmax = jmax-1;
		}}	
Selectclose("Result of VGLUT1");
Selectclose("MAP2");
Selectclose("VGLUT1");	}


//////////////////////////////////     Decoupage sur les différents channels des zones et création d'un stack    //////////
function StackMakers(imageTiff,ProtOrder,ImageName_woExt,mesurenonspec,output) { 
 
selectWindow(imageTiff);
getDimensions(w, h, channels, slices, frames);
for (i=1; i<=channels; i++){
	selectWindow(imageTiff);
	Stack.setChannel(i);
	run("Select None");
	run("Duplicate...", "title=Channel"+i+"_Couleur");
	img=getInfo("image.title");
	nbrNonspe=0; 
	
		for(j=0; j<roiManager("Count"); j++){
	 	selectWindow(img);
		wait(50);
		roiManager("Select", j);
		ROIcolor=Roi.getStrokeColor;
		if(ROIcolor=="blue"){
		run("Duplicate...", "title=Nonspe_"+j);
		nbrNonspe=nbrNonspe+1; 
		}else {
		run("Duplicate...", "title=Spot_"+j);}}

		run("Images to Stack", "  title=Spot_ use");
		saveAs("Tiff", output+File.separator+"Stack_"+ImageName_woExt+ProtOrder[i-1]+"_c"+i);

		if (nbrNonspe>1) {
		run("Images to Stack", "  title=Nonspe_ use");
		saveAs("Tiff", output+File.separator+"StackNonspe_"+ImageName_woExt+ProtOrder[i-1]+"_c"+i);
		}else {
		mesurenonspec=1;}
		Selectclose("Channel"+i+"_Couleur");
		}return mesurenonspec;}
		


		
////////////////////////////////////////////////////////////
function MeasureNonSpe (ImageName_woExt,ProtOrder,ChannelPattern,ChannelMeasure,pixelsize,output){

selectWindow("StackNonspe_"+ImageName_woExt+ProtOrder[ChannelPattern-1]+"_c"+ChannelPattern+".tif");
getDimensions(w, h, channels, slices, frames);
NbrSlices=slices;
run("Set Measurements...", "area mean min display redirect=None decimal=3");
roiManager("reset");

for (i = 1; i <=NbrSlices; i++) {
	resetThreshold();

	//creation rgn Pattern 
	selectWindow("StackNonspe_"+ImageName_woExt+ProtOrder[ChannelPattern-1]+"_c"+ChannelPattern+".tif");
	Stack.setSlice(i);
	name=getInfo("slice.label");
	run("Select None");
	setAutoThreshold("Mean dark");
	run("Create Selection");
	roiManager("Add");
	
	if (i == 1) {
	compt=0;
	}else {
	compt=compt+1;}

	//roiManager("Select", 0);
	roiManager("Select", compt);
	roiManager("Rename", "NonSpe_"+name);

	//creation rgn PatternOUT-Background

	Treshold("Mean");	
	compt=ROIaddRename(compt,"NonSpeBack_"+name); 

	//Measure sur Pattern
	selectWindow("StackNonspe_"+ImageName_woExt+ProtOrder[ChannelMeasure-1]+"_c"+ChannelMeasure+".tif");
	Stack.setSlice(i);
	run("Select None");
	run("Set Scale...", "distance=1 known="+pixelsize+" unit=µm");
	//roiManager("Select", newArray(0,1));
	roiManager("Select", newArray(compt-1,compt));
	roiManager("Measure");}
	roiManager("save", output+File.separator+"ROI_NonSpe_"+ProtOrder[ChannelMeasure-1]+ImageName_woExt+".zip");
	}

//////////////////////////Measure vGLUT1 ou PSD95 //////////////////////////////

function Measure(ImageName_woExt,ProtOrder,ChannelPattern,ChannelMeasure,pixelsize,output) { 
// function description

selectWindow("Stack_"+ImageName_woExt+ProtOrder[ChannelPattern-1]+"_c"+ChannelPattern+".tif");
getDimensions(w, h, channels, slices, frames);
NbrSlices=slices;
run("Set Measurements...", "area mean min display redirect=None decimal=3");
roiManager("reset");

for (i = 1; i <=NbrSlices; i++) {
	resetThreshold();
	
	//creation rgn Pattern 
	selectWindow("Stack_"+ImageName_woExt+ProtOrder[ChannelPattern-1]+"_c"+ChannelPattern+".tif");
	Stack.setSlice(i);
	name=getInfo("slice.label");
	Treshold("Mean dark");
	roiManager("Add");
	
	if (i == 1) {
	compt=0;
	}else {
	compt=compt+1;}

	roiManager("Select", compt);
	roiManager("Rename", "1Pattern_"+name);

	//creation rgn PatternOUT-Background
	Treshold("Mean");
	compt=ROIaddRename (compt,"PatternBack_"+name); 
	
	//Measure sur Pattern
	selectWindow("Stack_"+ImageName_woExt+ProtOrder[ChannelPattern-1]+"_c"+ChannelPattern+".tif");
	run("Set Scale...", "distance=1 known="+pixelsize+" unit=µm");
	roiManager("Select", newArray(compt-1,compt));
	roiManager("Measure");

	//Creation rgn Vglut 
	
	selectWindow("Stack_"+ImageName_woExt+ProtOrder[ChannelMeasure-1]+"_c"+ChannelMeasure+".tif");
	Stack.setSlice(i);
	Treshold("Mean dark");	
	compt=ROIaddRename(compt,"POI_All_"+name); 
	
	//Creation rgn Vglut back
	
	Treshold("Mean");
	compt=ROIaddRename(compt,"POI_Back_"+name); 


	//Creation rgn Vglut IN 
	roiManager("Select", (compt-3));
	run("Create Mask");
	rename("MaskPatterncontrol");
	setAutoThreshold("Default");
	run("Create Selection");
	getStatistics(area, mean, min, max, std, histogram);
	areaPatterncontrol=area;
	selectWindow("Stack_"+ImageName_woExt+ProtOrder[ChannelPattern-1]+"_c"+ChannelPattern+".tif");
	roiManager("Select", (compt-1));
	run("Create Mask");
	rename("MaskVglutAllcontrol");
	imageCalculator("Subtract create","MaskPatterncontrol","MaskVglutAllcontrol");
	rename("MaskResult");
	setAutoThreshold("Default");
	run("Create Selection");
	getStatistics(area, mean, min, max, std, histogram);
	areaVlut1=area;
	
	Selectclose ("MaskPatterncontrol");
	Selectclose ("MaskVglutAllcontrol");
	Selectclose ("MaskResult");
	
	if (areaVlut1==areaPatterncontrol) {
	print("Enrichissement non measureable pour "+ProtOrder[ChannelMeasure-1]+" "+ImageName_woExt+"_"+name);
	}else{
	
	//Creation rgn Protein IN
	compt=NewROI(compt, compt-3, compt-1, "POI_IN_"+name, "AND");

	AreaProtAll= RecupStat(compt-2) ; 
	AreaProtIN= RecupStat(compt) ; 
 
	if (AreaProtAll==AreaProtIN) {
		print("Enrichissement non measureable pour "+ProtOrder[ChannelMeasure-1]+" "+ImageName_woExt+"_"+name);
	}else {

	//Creation rgn Protein OUT
	compt=NewROI(compt,compt,compt-2,"POI_OUT_"+name, "XOR");
	
	//Measure sur Pattern
	selectWindow("Stack_"+ImageName_woExt+ProtOrder[ChannelMeasure-1]+"_c"+ChannelMeasure+".tif");
	run("Set Scale...", "distance=1 known="+pixelsize+" unit=µm");
	roiManager("Select", newArray(compt-3,compt-2,compt-1,compt));
	roiManager("Measure");
}}}
roiManager("save", output+File.separator+"ROI_measure_"+ProtOrder[ChannelMeasure-1]+"_"+ImageName_woExt+".zip");
saveAs("Results", output+File.separator+"_measure_"+ProtOrder[ChannelMeasure-1]+"_"+ImageName_woExt+".csv");}


//////////////////////////////////////////////////////////////////////////////
function Multisynapse(ImageName_woExt,ProtOrder,ChannelMeasure,pixelsize,output2,output) { 

selectWindow("Stack_"+ImageName_woExt+ProtOrder[ChannelMeasure-1]+"_c"+ChannelMeasure+".tif");
resetThreshold();
getDimensions(w, h, channels, slices, frames);
NbrSlices=slices;
run("Set Measurements...", "area mean min display redirect=None decimal=3");
roiManager("reset");

roiManager("open", output+File.separator+"ROI_measure_"+ProtOrder[ChannelMeasure-1]+"_"+ImageName_woExt+".zip");
NbrROI=roiManager("count");
	for (i=0; i<NbrROI; i++) {
	roiManager("select", i);
	Nom=Roi.getName;
	Condition=startsWith(Nom, "1P");
		if (Condition!=true) {
		roiManager("delete");
		NbrROI=roiManager("count");
		i=i-1;}}
	roiManager("save", output+File.separator+"ROI_Pattern_"+ProtOrder[ChannelMeasure-1]+"_"+ImageName_woExt+".zip");
	roiManager("reset");

for (i = 1; i <=NbrSlices; i++) {
	resetThreshold();
	selectWindow("Stack_"+ImageName_woExt+ProtOrder[ChannelMeasure-1]+"_c"+ChannelMeasure+".tif");
	run("Select None");
	Stack.setSlice(i);
	name=getInfo("slice.label");
	run("Duplicate...", "use");
	rename(name+"_duplicate");
	run("Duplicate...", "use");
	rename(name+"_duplicate_filter");
	run("Gaussian Blur...", "sigma=5");
	imageCalculator("Subtract create", name+"_duplicate",name+"_duplicate_filter") ;
	selectWindow("Result of "+name+"_duplicate");
	resetThreshold();
	setAutoThreshold("Default dark");
	run("Set Scale...", "distance=1 known="+pixelsize+" unit=µm");
	roiManager("open", output+File.separator+"ROI_Pattern_"+ProtOrder[ChannelMeasure-1]+"_"+ImageName_woExt+".zip");
	run("Select None");
	roiManager("select", i-1);
	roiManager("reset");
	run("Analyze Particles...", "add");
	RoiNBR=roiManager("count");
	if (RoiNBR>=1) {
		selectWindow("Stack_"+ImageName_woExt+ProtOrder[ChannelMeasure-1]+"_c"+ChannelMeasure+".tif");
		roiManager("deselect");
		roiManager("Measure");
		roiManager("save", output2+File.separator+"ROI_Multisynapse_"+name+"_"+ProtOrder[ChannelMeasure-1]+"_"+ImageName_woExt+".zip");
		Selectclose(name+"_duplicate");
		Selectclose(name+"_duplicate_filter");
		Selectclose("Result of "+name+"_duplicate");
		roiManager("reset");
	}else{
		Selectclose(name+"_duplicate");
		Selectclose(name+"_duplicate_filter");
		Selectclose("Result of "+name+"_duplicate");
	}}
	saveAs("Results", output+File.separator+"_measureMultisynapse_"+ProtOrder[ChannelMeasure-1]+"_"+ImageName_woExt+".csv");}
	
	
//////////////////////////////////////////////////////////////////////////////

function Selectclose(Name) { 
	selectWindow(Name);
	run("Select None");
	close(); }
	
function RecupStat(ROInbr) {
	roiManager("Select", ROInbr);
	getStatistics(area, mean, min, max, std, histogram);
	return area; }
	
function Treshold(Option) {
	run("Select None");
	resetThreshold();
	setAutoThreshold(Option);
	run("Create Selection");  }
	
function ROIaddRename(compteur,ROIName) {	
	roiManager("Add");
	compteur=compteur+1;
	roiManager("Select", compteur); 
	roiManager("Rename", ROIName);
	return compteur}	
	
function NewROI(compteur,ROI1,ROI2,ROIName,type) {
	roiManager("Select", newArray(ROI1,ROI2));
	roiManager(type);
	compteurs=ROIaddRename(compteur, ROIName);
	return compteurs}
	

