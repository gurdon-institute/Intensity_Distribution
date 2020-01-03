	/******variables to adjust******/
var dotRadius = 0.65;		//the estimated radius of dots in C2
var dotThreshold = 400;		//the minimum intensity required for dots in C2
var minNucleusArea = 100;	//the minimum area required for a C1 nucleus to be included
var doFilter = true;		//do filtering to exclude nuclei with <1 dots?
var normalise = true;		//normalise intensity to fold difference?
var pixelW = 0.3229;		//the pixel width (µm) calibration to apply

//var counteq = " y = a*x+b";	//straight line
//var counteq = " y = a*exp(pow(x-b,2)/c) + d*exp(pow(x-e ,2)/f) "; //bimodal
//var counteq = " y = a + b*x + c*pow(x,2) + d*pow(x,3) + e*pow(x,4) "; //4th degree polynomial
//var counteq = " y = a + b*x + c*pow(x,2) + d*pow(x,3) "; //3rd degree polynomial
var counteq = " y = (1/a*(2*"+PI+"))*exp( -( pow(x-b,2) / (2*pow(c,2)) ) ) + d"; //Gaussian

var doteq = " y = a*x+b";	//straight line
//var doteq = " y = a*pow(b,x) ";	//exponential
//var doteq = " y = a*pow(x,2) + b "; //2nd degree polynomial
	/******variables to adjust******/

var unit = "µm";
var dotArea = pow(PI*dotRadius,2);	//dot radius ~2 pixels
var alf = newArray("a","b","c","d","e ","f");	//letters for Fit parameters
path = getDirectory("Choose a Directory");
list = getFileList(path);
run("Set Measurements...", "area centroid redirect=None decimal=5");
setOption("ShowRowNumbers", false);
if(nResults()>0){selectWindow("Results");run("Close");}
//run("Bio-Formats Macro Extensions");
setBatchMode(true);
for(i=0;i<list.length;i++){
	if(endsWith(list[i], "tiff")){
		open(path+list[i]);
		measure();
	}
	//else if(endsWith(list[i], "flex")){
	//	Ext.openImagePlus(path+list[i]);
	//	measure();
	//}
}
if(nResults()>0){analyse();}
setBatchMode("exit and display");

function analyse(){
	var binnedCounts = newArray();
	var binnedDots = newArray();
	var Xbins = newArray();
	var max = pow(2,12); //images are 12-bit
	//var nBins = floor(2*pow(max,1/3))+1;	//Rice rule
	var nBins = sqrt(max)*2;
	var binSize = max/nBins;
	
	for(a=0;a<nBins;a++){
		binnedCounts = Array.concat(binnedCounts,0);
		binnedDots = Array.concat(binnedDots,0);
		Xbins = Array.concat(Xbins,(a*binSize));
	}
	var XMin = pow(2,12);
	var XMax = -1;
	var atLeastOneDotN = 0;
	for(i=0;i<nResults();i++){
		bindex = round(getResult("C1 Mean",i)/max*nBins);
		XMin = minOf(XMin,Xbins[bindex]);
		XMax = maxOf(XMax,Xbins[bindex]);
		binnedCounts[bindex]++;
		thisN = getResult("C2 Dots",i);
		binnedDots[bindex] += thisN;
		if(thisN>=1){atLeastOneDotN++;}
	}
	print(path+" : number of nuclei with at least one dot = "+atLeastOneDotN);
	for(d=0;d<binnedDots.length;d++){
		binnedDots[d] = binnedDots[d]/binnedCounts[d];
		if(isNaN(binnedDots[d])){binnedDots[d]=0;}
	}

	if(doFilter){//filter arrays
		countF = newArray();
		dotF = newArray();
		XF = newArray();
		XMin = pow(2,12);
		XMax = -1;
		for(f=0;f<binnedDots.length;f++){
			if(binnedDots[f]>1){
				countF = Array.concat(countF,binnedCounts[f]);
				dotF = Array.concat(dotF,binnedDots[f]);
				XF = Array.concat(XF,(f*binSize));
				XMin = minOf(XMin,(f*binSize));
				XMax = maxOf(XMax,(f*binSize));
			}
		}
		binnedCounts = Array.copy(countF);	//swap in filtered versions
		binnedDots = Array.copy(dotF);
		Xbins = Array.copy(XF);
	}

	Array.getStatistics(binnedCounts, Ymin1, YMax1, Ymean1, YstdDev1);
	Array.getStatistics(binnedDots, Ymin2, YMax2, Ymean2, YstdDev2);
	YMax = maxOf(YMax1,YMax2);

	end = (XMax/binSize)+1;
	binnedCounts = Array.trim(binnedCounts,end);
	binnedDots = Array.trim(binnedDots,end);
	Xbins = Array.trim(Xbins,end);

	if(normalise){
		XMax = 0;
		for(i=0;i<Xbins.length;i++){
			Xbins[i] = (Xbins[i]/XMin);
			XMax = maxOf(XMax,Xbins[i]);
		}
		XMin = 1;
	}

Fit.doFit(counteq,Xbins,binnedCounts);
fitCounts = newArray();
fitX = newArray();
for(f=0;f<=Xbins[Xbins.length-1]+1;f+=0.1){
	//print(f+"/"+Xbins[Xbins.length-1]+" : "+Fit.f(f));
	fitX = Array.concat(fitX,f);
	fitCounts = Array.concat(fitCounts,Fit.f(f));
}
countR2 = Fit.rSquared();
for(p=0;p<Fit.nParams();p++){
	counteq = replace(counteq,alf[p],Fit.p(p));
}
counteq = replace(counteq,"\\+ \\-","-");
Fit.doFit(doteq,Xbins,binnedDots);
fitDots = newArray();
for(f=0;f<=Xbins[Xbins.length-1]+1;f+=0.1){
	//print(f+" : "+Fit.f(f));
	fitDots = Array.concat(fitDots,Fit.f(f));
}
dotR2 = Fit.rSquared();
for(p=0;p<Fit.nParams();p++){
	doteq = replace(doteq,alf[p],Fit.p(p));
}
doteq = replace(doteq,"\\+ \\-","-");
	Plot.create("Intensity Distribution : "+path,"C1 Intensity fold difference from minimum","n");
	Plot.setFrameSize(600, 600);
	Plot.setLimits(XMin,XMax,0,YMax*1.1);
	Plot.setColor("blue");
	Plot.add("x",Xbins,binnedCounts);
Plot.add("line",fitX,fitCounts);
	Plot.addText("x - Number of Nuclei",0.01,0.03);
	Plot.addText(counteq+" R^2="+countR2,0.01,0.06);
	Plot.setColor("red");
	Plot.add("circles",Xbins,binnedDots);
Plot.add("line",fitX,fitDots);
	Plot.addText("o - Mean Number of Dots per Nucleus",0.01,0.09);
	Plot.addText(doteq+" R^2="+dotR2,0.01,0.12);
	Plot.show();
}

function measure(){
	title = getTitle();
	getDimensions(W,H,C,Z,T);
	if(Z!=2){return;}	//stacks are messed up, slices=channels
	setVoxelSize(pixelW, pixelW, 1, unit);
	Stack.setSlice(1);
	run("Duplicate...", "title=nuclei");
	selectWindow("nuclei");
	run("Gaussian Blur...", "sigma=3");
	setAutoThreshold("Li dark");
	setOption("BlackBackground", true);
	run("Convert to Mask");
	run("Fill Holes");
	run("Open");
	run("Watershed");
	if(roiManager("count")>0){
		roiManager("Deselect");
		roiManager("Delete");
	}
	run("Analyze Particles...", "size="+minNucleusArea+"-Infinity display exclude add");
	selectWindow("nuclei");	close();
	for(r=0;r<roiManager("count");r++){
		res = nResults()-roiManager("count")+r;
		selectWindow(title);
		//C1 mean + coords
		Stack.setSlice(1);
		roiManager("Select", r);
		run("Add Selection...");
		getStatistics(null,C1mean);
		setResult("C1 Mean",res,C1mean);
		Roi.getBounds(roiX, roiY, roiW, roiH);
		//C2 mean
		Stack.setSlice(2);
		roiManager("Select", r);
		getStatistics(null,C2mean);
		setResult("C2 Mean",res,C2mean);
		//C2 dot count
		run("Duplicate...", "title=C2"+r);
		selectWindow("C2"+r);
		run("Select None");
		setThreshold(dotThreshold, pow(2,16));
		run("Create Selection");
		if(selectionType!=-1){getStatistics(C2area);}
		else{C2area=0;}
		setResult("C2 Area I>="+dotThreshold,res,C2area);
		selectWindow("C2"+r);	close();
		nDots = C2area/dotArea;
		setResult("C2 Dots",res,nDots);
		
	}
	updateResults();
	run("Select None");
	selectWindow(title);	close();
}