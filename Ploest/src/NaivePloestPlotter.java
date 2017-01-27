import java.awt.Color;
import java.awt.Shape;
import java.awt.geom.Rectangle2D;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.io.UnsupportedEncodingException;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.Vector;

import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartUtilities;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.axis.NumberAxis;
import org.jfree.chart.axis.ValueAxis;
import org.jfree.chart.plot.DatasetRenderingOrder;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.renderer.xy.StandardXYItemRenderer;
import org.jfree.chart.renderer.xy.XYItemRenderer;
import org.jfree.chart.renderer.xy.XYLineAndShapeRenderer;
import org.jfree.data.xy.XYDataset;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;

import dataFitters.GaussianDataFitter;
import dataFitters.PoissonDataFitter;
import jMEF.MixtureModel;
import jMEF.PVector;

public class NaivePloestPlotter {
	static boolean PLOIDY_CONTINUITY=false;//TO DO smooths the ploidy estimation by continuity filter
	static boolean AVG_PLOIDY=true;//smooths the ploidy estimation by averaging over a window of length: PLOIDY_SMOOTHER_WIDTH
	static int PLOIDY_SMOOTHER_WIDTH=49;//preferably odd number
	
	static float[] readCounts;
	static final int MAX_NB_MIXTURES=10;
	Map<String,ContigData> contigsList;
	List<String> contArrList;
	static double[] clusterMus;//final result of the MEANS (mus) of the clusters in the mixture model fitting
	static PoissonMixturePDF pmPDFResult;//final POISSON Mixture PDF
	static NaivePDF npdf;//Naive Smoothed Density Function
	PVector[] fitPoints;
	JFreeChart chart;
	static int maxX=0;
	//static int maxY=0;
	static int totalDataPoints=0;//total number of input datapoints (coverage for all windows)

	static int finalNumberOfMixtures;
	RatioFindNaive rt;//contains the ratio of each cluster to the ploidy-unit which allows computation of ploidy from contig coverage
	String writeoutStringLog="";
	
	ArrayList<ContigData>  unsolvedPloidyContigs=new ArrayList<ContigData> () ;
	ArrayList<ContigData>  removedContigs=new ArrayList<ContigData> () ;
	
	public NaivePloestPlotter(Map<String,ContigData> contList,int maxWindows, float[] rc) {
		readCounts=rc;
		contigsList=contList;
		contArrList = new ArrayList<String>(contigsList.keySet());

		
		try{
			displayScatterPlots();//Scatterplot containing the contig coverage (per slided window)		
			createFitterDataset() ;
			fitNaiveMixtureModel();	//approximate the distribution by naive smoother and infere max points				
			displayPloidyAndCoveragePlotNaive();//Plot containng both the coverage and the ploidy estimation
			
			
		}catch (Exception e){
			System.err.println("Error in PloestPlotter constructor");
		}
		
	}


	
	
	private void writeOutPloEstByFragment(PrintWriter writer , XYSeries series, String contigname ) {
		//System.out.println("-------------writeOutPloEstByFragment------------------");
		writer.println("Detailed ploidy estimation for "+contigname);
		//System.out.print("Detailed ploidy estimation for "+contigname);
		if(series.getItemCount()>0){
			Number prevPloidy=series.getMinY();//series.getY(0);
			Number prevPos=0;
			int ItemsSize=series.getItems().size();
			//System.out.println("prevPloidy "+prevPloidy+" ItemsSize:"+ItemsSize);
			for (int yv=1;yv<ItemsSize;yv++){
				//System.out.print(" "+yv+","+series.getY(yv));
				if(!series.getY(yv).equals(prevPloidy)){//segmentation point
					writer.println(" Ploidy: "+prevPloidy+" from +/- "+((int)prevPos*SamParser.windowLength/2)+" bp to +/- "+(yv*SamParser.windowLength/2)+" bp");
					//System.out.println(" Ploidy: "+prevPloidy+" from +/- "+((int)prevPos*SamParser.windowLength/2)+" bp to +/- "+(yv*SamParser.windowLength/2)+" bp");
					prevPloidy=series.getY(yv);
					prevPos=yv;
				}
		
				
			}

			//System.out.println(" Ploidy: "+prevPloidy+" from +/- "+((int)prevPos*SamParser.windowLength/2)+" bp to +/- "+(    (ItemsSize*SamParser.windowLength/2)+(SamParser.windowLength/2)  )+" bp");
			writer.println(" Ploidy: "+prevPloidy+" from +/- "+((int)prevPos*SamParser.windowLength/2)+" bp to +/- "+(    (ItemsSize*SamParser.windowLength/2)+(SamParser.windowLength/2)  )+" bp");
			writer.println();
			//System.out.println();
		}else{
			writer.println("  Not enough information to determine ploidy for "+contigname);
			writer.println();
			System.out.println("  Not enough information to determine ploidy for "+contigname);
		}
		

		//System.out.println("-------------writeOutPloEstByFragment END------------------");
		
	}




	public void fitNaiveMixtureModel(){
		System.out.println("-------------Smoothing the data by averaging values in bins------------------");
		
		npdf=new NaivePDF(readCounts);
		int k=significantMaxsInPDF(npdf);
		SamParser.barchart.BarChartWithFit(npdf,"FINALRESULT");
		rt=new RatioFindNaive(clusterMus);
		rt.writeOut();
	}
	
	
	public int indexOfMode(int[] vector){//finds the index of the most represented result in this vector
		double max=vector[0] ;
		int maxIndex=0;
		for (int ktr = 0; ktr < vector.length; ktr++) {
			if (vector[ktr] > max) {
				maxIndex=ktr;
				max=vector[ktr] ;
			}
		}
		return maxIndex;
	}

	public void displayScatterPlots() throws IOException{
		ContigData contigD;
		for (int c=0;c<contigsList.size();c++){//for each contig
			contigD=contigsList.get(contArrList.get(c));
			XYDataset data1=createPlotDataset(contigD);
			if(maxX>0){
				chart = ChartFactory.createScatterPlot(
						("Genome Coverage "+contigD.contigName), // chart title
						"Genome Position (x " +(contigD.windLength/2)+" bp)", // x axis label
						"Coverage", // y axis label
						data1, // XYDataset 
						PlotOrientation.VERTICAL,
						true, // include legend
						true, // tooltips
						false // urls
						);
				//Set range
				XYPlot xyPlot = (XYPlot) chart.getPlot();
				NumberAxis domain = (NumberAxis) xyPlot.getDomainAxis();
				domain.setRange(0.00, maxX);
				ValueAxis rangeAxis = xyPlot.getRangeAxis();	
			
				//rangeAxis.setRange(0.00, SamParser.readsDistributionMaxCoverage);
				if(contigD.maxY>0){
					rangeAxis.setRange(0.00,contigD.maxY);
				}else {
					rangeAxis.setRange(0.00,10);
					System.err.println(contigD.contigName+" doesn't have any coverage. Contig is Removed");
					
					contigsList.remove(contArrList.get(c));
					contArrList.remove(c);
					c--;
				}
				
				
				
				ChartUtilities.saveChartAsJPEG(new File(Ploest.outputFile + "//" + Ploest.projectName+ "//Contig_Coverage_Charts//Chart_Contig_"+c+".jpg"), chart, 1000, 600);
				
			}else {
				System.err.println(contigD.contigName+" is too short ("+contigD.maxLength+"bp ) for the length of the sliding window (winLength="+contigD.windLength+" bp). Contig is Removed");
				writeoutStringLog+= contigD.contigName+" is too short ("+contigD.maxLength+"bp ) for the length of the sliding window (winLength="+contigD.windLength+" bp). Contig is Removed\r\n";
			}
			
		}
		System.out.println("ScatterPlots displayed");

	}


	

public void displayPloidyAndCoveragePlotNaive()throws IOException{
		

		ContigData contigD;
		PrintWriter writer = rt.writer;
		writer.println();
		writer.println("*********************************************************");
		writer.println(writeoutStringLog);
		writer.println("*********************************************************");
		writer.println(" Ploidy estimation detailed by fragments.  Precision: +/-"+(SamParser.windowLength/2)+" bp");
		writer.println();
		for (int c=0;c<contigsList.size();c++){//for each contig
			
			contigD=contigsList.get(contArrList.get(c));
			XYPlot xyPlot = new XYPlot();
			/* SETUP SCATTER */

			// Create the scatter data, renderer, and axis
			XYDataset collection1 = createPlotDataset(contigD);
			XYItemRenderer renderer1 = new XYLineAndShapeRenderer(false, true);   // Shapes only
			ValueAxis domain1 = new NumberAxis("Genome Position (x " +(contigD.windLength/2)+" bp)");
			ValueAxis rangeAxis = new NumberAxis("Coverage");

			// Set the scatter data, renderer, and axis into plot
			xyPlot.setDataset(0, collection1);
			xyPlot.setRenderer(0, renderer1);
			xyPlot.setDomainAxis(0, domain1);
			xyPlot.setRangeAxis(0, rangeAxis);

			// Map the scatter to the first Domain and first Range
			xyPlot.mapDatasetToDomainAxis(0, 0);
			xyPlot.mapDatasetToRangeAxis(0, 0);

			/* SETUP LINE */

			// Create the line data, renderer, and axis
			XYDataset collection2 = createPloidyEstimationDatasetNaive(contigD,writer);
			if(!unsolvedPloidyContigs.contains(contigD)){
				XYItemRenderer renderer2 = new XYLineAndShapeRenderer(false , true);   // Lines only
				ValueAxis domain2 = new NumberAxis("Genome Position (x " +(contigD.windLength/2)+" bp)");
				ValueAxis range2 = new NumberAxis("Ploidy Estimation");
				range2.setUpperBound(rangeAxis.getUpperBound()/rt.bestScore.candidateUnit);
				domain2.setUpperBound(domain1.getUpperBound());
				// Set the line data, renderer, and axis into plot
				range2.setStandardTickUnits(NumberAxis.createIntegerTickUnits());//domain2.setStandardTickUnits(NumberAxis.createIntegerTickUnits());
				xyPlot.setDataset(1, collection2);
				xyPlot.setRenderer(1, renderer2);
				xyPlot.setDomainAxis(1, domain2);
				xyPlot.setRangeAxis(1, range2);
				// Map the line to the second Domain and second Range
				xyPlot.mapDatasetToDomainAxis(1, 1);
				xyPlot.mapDatasetToRangeAxis(1, 1);
				xyPlot.setDatasetRenderingOrder( DatasetRenderingOrder.FORWARD );	
				// Create the chart with the plot and a legend
				JFreeChart chart = new JFreeChart("Coverage and Ploidy Estimation :"+contigD.contigName, JFreeChart.DEFAULT_TITLE_FONT, xyPlot, true);
				String correctedContigName = contigD.contigName.replaceAll("[^a-zA-Z0-9.-]", "_");
				ChartUtilities.saveChartAsJPEG(new File(Ploest.outputFile + "//" + Ploest.projectName+ "//Ploidy_Estimation_Charts//Ploidy_Estimation_"+correctedContigName+".jpg"),chart, 1500, 900);
			}
		}
		if (unsolvedPloidyContigs.size()>0){
			int minLength=unsolvedPloidyContigs.get(0).maxLength;
			writer.println("> The following contigs could not be solved with the current window length of "+Ploest.windowLength+" bp : ");
			for (int u=0;u<unsolvedPloidyContigs.size();u++){
				writer.println("  "+unsolvedPloidyContigs.get(u).contigName);
				if (minLength>=unsolvedPloidyContigs.get(u).maxLength){
					minLength=unsolvedPloidyContigs.get(u).maxLength;
				}
			}
			writer.println(" A new estimation will be attempted with a shorter window length of "+minLength/10);
			Ploest.windowLength=minLength/10;
		}
		
		
		if (removedContigs.size()>0){
			writer.println();
			writer.println("> The following contigs were removed because of lack of data values: ");
			for (int u=0;u<removedContigs.size();u++){
				writer.println("  "+removedContigs.get(u).contigName);
			}
			writer.println(" ");
			SamParser.RUN_SECOND_ROUND=true;	
		}
		
		writer.close();
	}


	public int findIndexOfMin(double[] bicVector){
		double min=bicVector[1] ;
		int minIndex=1;
		for (int ktr = 1; ktr < bicVector.length; ktr++) {
			if ((!Double.isNaN(bicVector[ktr]))&&(bicVector[ktr]>0)&&(bicVector[ktr] < min)) {
				minIndex=ktr;
				min=bicVector[ktr] ;
			}
		}
		return minIndex;
	}

	private static XYDataset createPlotDataset(ContigData contigD) throws FileNotFoundException, UnsupportedEncodingException {
		XYSeriesCollection result = new XYSeriesCollection();
		XYSeries series = new XYSeries(" Coverage");

		//PrintWriter writer = new PrintWriter(Ploest.outputFile + "//" + Ploest.projectName+ "//plotDataSet"+maxX+".txt", "UTF-8");

		double x;
		double y;
		
		//this first loop is for the PLOESTPLOTTER
		int wInd=0;//writting index
	
		for (int i = 0; i <= (contigD.windPos.size()-1); i++) {
			if(contigD.windPos.get(i)!=null){
				x = wInd++;  			
				y = contigD.windPos.get(i);
				if(y>contigD.maxY)contigD.maxY=(int)y;

				series.add(x, y);
				//writer.println( " x:" +x + " y:"+y);
			}

		}

		result.addSeries(series);
		maxX=contigD.windPos.size();
		//writer.close();
		return result;
	}


	private  void createFitterDataset() {
	
		fitPoints=SamParser.fitPoints;

	}



	
	private  XYDataset createPloidyEstimationDatasetNaive(ContigData contigD, PrintWriter writer ) {

		XYSeriesCollection result = new XYSeriesCollection();
		XYSeries series = new XYSeries(" Ploidy Estimation 1");
		
		double [] xValues=new double[contigD.windPos.size()];
		int [] yValues=new int[contigD.windPos.size()];
		maxX=0;
		int wInd=0;
		//this loop is for ESTIMATED PLOIDY PLOTTING
		//System.out.println("ESTIMATED PLOIDY PLOTTING:");
		for (int i = 0; i < contigD.windPos.size(); i++) {

			if(contigD.windPos.get(i)!=null){
				xValues [wInd]= wInd;  			
				yValues [wInd]= getPointPloidyEstimationNaive(contigD.windPos.get(i));
				//System.out.print(" "+wInd+","+yValues [wInd]);
				if (!AVG_PLOIDY)series.add(wInd, yValues [wInd]);
				wInd++;
			}				

		}//System.out.println();
		if (--wInd>maxX){
			maxX=(int) wInd;
		}
		System.out.println("createPloidyEstimationDatasetNaive :"+contigD.contigName);
		if(AVG_PLOIDY ){		//smooth the ploidy plot by averaging the values over a window of PLOIDY_SMOOTHER_WIDTH points
			series=averagePloidyMode(contigD,series,wInd,xValues,yValues);//uses the mode over PLOIDY_SMOOTHER points
		}
		System.out.println();
		if(series.getItemCount()<5 && !removedContigs.contains(contigD)){
			unsolvedPloidyContigs.add(contigD);
			System.err.println(" CONTIG :"+contigD.contigName+" added to list for new window length ploidy estimation.");
			
		}

		result.addSeries(series);
	
		//
	
		writeOutPloEstByFragment(writer, series,contigD.contigName );//writes out the ploidy estimation detailed by fragment

		
		return result;
	}

	
	public XYSeries averagePloidyMode(ContigData contigD,XYSeries series, int wInd,double [] xValues,int [] yValues){
		String contigname=contigD.contigName;
		
		System.out.println(" AVERAGE PLOIDY MODE. wInd="+wInd+" (PLOIDY_SMOOTHER_WIDTH/2):"+(PLOIDY_SMOOTHER_WIDTH/2));
		
		//wInd == size of x data
		int currentMode=0;//the most observed ploidy value over the PLOIDY_SMOOTHER_WIDTH
		int PLOIDY_SMOOTHER_WING=PLOIDY_SMOOTHER_WIDTH/2; //length of each of the sides of the PLOIDY_SMOOTHER window before and after the position being evaluated
		int [] ploidyCounter=new int[MAX_NB_MIXTURES+1];//over the PLOIDY_SMOOTHER_WIDTH, this vector keeps track of how many times each ploidy is observed
		if (wInd > PLOIDY_SMOOTHER_WIDTH) {//we need a minimum of points to average the ploidy
			
			// solve the first  positions of the plot
			for (int v = 0; v < PLOIDY_SMOOTHER_WIDTH/2; v++) {
				if(yValues[v]<=MAX_NB_MIXTURES ){
					ploidyCounter[yValues[v]]++;
					if (ploidyCounter[currentMode]<ploidyCounter[yValues[v]]){
						currentMode=yValues[v];
					}
					if (currentMode!=0){
						series.add(xValues[v], currentMode);
						//System.out.print(" "+(int)xValues[v]+","+currentMode);
					}
				}
				
			}
			//solve the rest of the genome until the last positions-PLOIDY_SMOOTHER_WIDTH/2
			int mostRightValue;//value at the right end of the ploidy-smoother window
			int mostLeftValue=0;//value at the left end of the ploidy-smoother window
			for(int v=(PLOIDY_SMOOTHER_WIDTH/2);v<(wInd-(PLOIDY_SMOOTHER_WIDTH/2));v++){
				if(yValues[v]<=MAX_NB_MIXTURES ){
		
					mostLeftValue=yValues[v-(PLOIDY_SMOOTHER_WIDTH/2)];
					if(mostLeftValue>0 && mostLeftValue<=MAX_NB_MIXTURES ){
						ploidyCounter[mostLeftValue]--;//remove mostleft value of window
					}
		
					mostRightValue=yValues[v+(PLOIDY_SMOOTHER_WIDTH/2)];
					if (mostRightValue>0 && mostRightValue<=MAX_NB_MIXTURES )ploidyCounter[mostRightValue]++;//add mostright value of window
			
					if(mostLeftValue!=mostRightValue){//if the removed and the added are different, get the mode of the curent vector ploidyCounter
						for(int pc=0;pc<ploidyCounter.length;pc++){
							if(ploidyCounter[pc]>ploidyCounter[currentMode])currentMode=pc;
						}
					}
					if (currentMode!=0){
						series.add(xValues[v], currentMode);
						//System.out.print(" "+(int)xValues[v]+","+currentMode);
					}
				}
			}
			
			//solve the very last positions PLOIDY_SMOOTHER_WIDTH/2
			for(int v=(wInd-(PLOIDY_SMOOTHER_WIDTH/2));v<wInd;v++){
				if(yValues[v]<=MAX_NB_MIXTURES ){
					if(yValues[v]>0 && yValues[v]<=MAX_NB_MIXTURES)ploidyCounter[yValues[v]]--;
					if (ploidyCounter[currentMode]<ploidyCounter[yValues[v]]){
						for(int pc=0;pc<ploidyCounter.length;pc++){
							if(ploidyCounter[pc]>ploidyCounter[currentMode])currentMode=pc;
						}
					}
					if (currentMode!=0){
						series.add(xValues[v], currentMode);
						//System.out.print(" "+(int)xValues[v]+","+currentMode);
					}
				}
			}
			
			System.out.println("");
		}else{//not enough points, simply average over the available points
			for (int v = 0; v < wInd; v++) {
				if(yValues[v]<=MAX_NB_MIXTURES ){
					ploidyCounter[yValues[v]]++;
					if (ploidyCounter[currentMode]<ploidyCounter[yValues[v]]){
						currentMode=yValues[v];
					}
					
					if (currentMode!=0){
						series.add(xValues[v], currentMode);
						//System.out.print(" "+(int)xValues[v]+","+currentMode);
					}
				}
			}
			System.out.println();
		}
		System.out.print("END. NOW CHECK CoNTINUITY.  ");
		if(series.getItemCount()>0){
			return checkContinuity(series,20);
		}else{
			removedContigs.add(contigD);
			
			System.err.println("Error in contig :"+contigname+". This series have 0 values!!!");
			return series;
		}
		
	}
	
	private XYSeries checkContinuity(XYSeries series, int continuityLength) {
		XYSeries result=new XYSeries(" Ploidy Estimation");
		int [] yvalues=new int[series.getItemCount()];
		Number currentY=series.getY(0);//current Y value being counted
		int continousCurrents=0;//nb of contiguous observed values of the current Y value
		Number[] firstObservedValues=new Number[continuityLength];//first observed values are stored while verifying that there are at least continuityLength contiguous observations
		for (int i=0;i<series.getItemCount();i++){
			if( currentY.equals(series.getY(i) )){//if Y values is currently being observed
				continousCurrents++;
				if (continousCurrents>continuityLength){//if the continuity length is respected
					result.add(series.getX(i),currentY);//add to series
				}else if (continousCurrents==continuityLength){//the following 2 'else if' allows the first continuityLength observed values to be taken into account 
					firstObservedValues[0]=currentY;
					for (int j=0;j<continuityLength;j++){//add the stored first continuityLength observed values
						result.add((series.getX(i).intValue()-continuityLength+j+1),firstObservedValues[j]);//add to series	
					}
				}else if (continousCurrents<continuityLength){
					firstObservedValues[continousCurrents]=currentY;
				}
				
			}else{//new Y values observed
				
				continousCurrents=0;//reset counter
				currentY=series.getY(i);//reset observed Y
			}
		}
		System.out.println("END CoNTINUITY CHECK");
		return result;
	}




	public int getPointPloidyEstimationNaive(double ptCoverage){
			return(int) Math.round(ptCoverage/RatioFindNaive.candUnit);
	}


	static int getFinalNumberOfMixtures(){
		return finalNumberOfMixtures;
	}

	

	private int significantMaxsInPDF(NaivePDF naivePDF) {
		int sigMaxs = 0;//nb of significant maximums
		double threshold = SamParser.readDistributionMaxY * 0.001;// discards the values that are
													// below this threshold
		System.out.println(" SamParser.maxYh :" + SamParser.readDistributionMaxY + " Y min threshold:" + threshold);

		ArrayList<Double> yMinList = new ArrayList<Double>();
		ArrayList<Double> xMinList = new ArrayList<Double>();
		//pointers to the y values
		double left = 0;
		double mid = 0;
		double right = 0;
		//pointers to the x values
		int ind = 0;
		int lastLeftIndex = 0;
		int lastMidIndex = 0;
		int lastRightIndex = 0;

		System.out.println(" pmf.yDataPoints.length :" + naivePDF.yDataPoints.length + " Y min threshold:" + threshold+ " \nSignificant maxima in NaivePDF:");

		while (ind < naivePDF.yDataPoints.length) {

			if (naivePDF.yDataPoints[ind] != right) {
				//move and update pointers
				left = mid;
				mid = right;
				right = naivePDF.yDataPoints[ind];
				lastLeftIndex = lastMidIndex;
				lastMidIndex = lastRightIndex;
				lastRightIndex = ind;
				//System.out.println(" .     ind :" + ind + " = " + mid + "  l:" + left + " m:" + mid+ " r:" + right );

				if (right < mid && mid > left && mid > threshold) {// we count maxs only above threshold
					
					// now that we encountered a max bin, we scan the previous
					// and the following bins to find the exact maximum point in this area
					
					double maxVal = 0;// precise max Y value in the corresponding bins
					int Xindex = lastLeftIndex;//index (x value) of the maxVal

					//System.out.println(" ....    max in :" + lastRightIndex + " = " + mid + "  l:" + left + " m:" + mid+ " r:" + right + "  between leftIn:" + lastLeftIndex + " midInd:" + lastMidIndex+ " rightInd:" + lastRightIndex + " \nSCANING from lastLeftIndex:" + lastLeftIndex+ " to  :" + (lastRightIndex + naivePDF.smootherLength));

					for (int ib = lastLeftIndex; ib < (lastRightIndex + naivePDF.smootherLength); ib++) {

						if (readCounts[ib] > maxVal) {
							maxVal = readCounts[ib];
							Xindex = ib;
						}
						//System.out.println("     ib:" + ib + " Xindex:" + Xindex + " maxval:" + maxVal+ " readCounts[ib]:" + readCounts[ib]);

					}
					if(!xMinList.contains(naivePDF.xDataPoints[Xindex])){
						yMinList.add(maxVal);
						xMinList.add(naivePDF.xDataPoints[Xindex]);
						System.out.println(" ****    max in :" + naivePDF.xDataPoints[Xindex] + " = " + maxVal+ " lastLeftIndex:"+lastLeftIndex+" lastRightIndex:"+lastRightIndex+" (+ "+naivePDF.smootherLength+" = "+(lastRightIndex + naivePDF.smootherLength)+")");
						sigMaxs++;
					}else System.out.println(" REPEATED max in :" + naivePDF.xDataPoints[Xindex] + " = " + maxVal+"  ...VALUE DISCARDED...!");

					

				}
			} else {
				right = naivePDF.yDataPoints[ind];
				lastRightIndex = ind;
			}

			ind++;
		}

		System.out.println("     SIGMAXS :" + sigMaxs);
		clusterMus=new double[xMinList.size()];
		System.out.println("     CLUSTER MUS :" );
		for (int c=0;c<xMinList.size();c++){
			clusterMus[c]=xMinList.get(c);
			System.out.print("  "+clusterMus[c] );
		}
		System.out.println();

		return sigMaxs;
	}
}
