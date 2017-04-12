package abeellab.pedca;
import java.awt.Color;
import java.awt.Font;
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

public class NaivePedcaPlotter {
	
	
	
	//static variable for smoothing the ploidy estimation
	static boolean AVG_PLOIDY=true;//smooths the ploidy estimation by averaging over a window of length: PLOIDY_SMOOTHER_WIDTH
	static int PLOIDY_SMOOTHER_WIDTH=Pedca.k;//preferably odd number
	static int CONTINUITY_POINTS=20;//minimum datapoints (windows) required to evaluate the continuity of the ploidy estimation
	
	//variable for the read count distribution
	static float[] readCounts;
	static final int MAX_PLOIDY=Pedca.MAX_PLOIDY;
	Map<String,ContigData> contigsList;
	List<String> contArrList;
	static double[] clusterMus;//final result of the MEANS (mus) of the clusters in the mixture model fitting
	int sumOfReadCounts=0;
	//variables for the naive smooth
	static NaivePDF npdf;//Naive Smoothed Density Function
	
	//PVector[] fitPoints;
	int[] intFitPoints;
	JFreeChart chart;
	static int maxX;
	static int totalDataPoints;//total number of input datapoints (coverage for all windows)
	static int finalNumberOfMixtures;
	RatioFindNaive rt;//contains the ratio of each cluster to the ploidy-unit which allows computation of ploidy from contig coverage
	String writeoutStringLog="";
	
	//variable for plotting the ploidy estimation
	//variables for managing unsolved ploidy contigs and run a second round with a smaller window length
	ArrayList<ContigData>  unsolvedPloidyContigs=new ArrayList<ContigData> () ;
	ArrayList<ContigData>  removedContigs=new ArrayList<ContigData> () ;
	
	//variable for base calling upon the contigs with ploidy from the first and second cluster
    ArrayList<ContigData> continousPloidyContigsCluster1=new ArrayList<ContigData> () ; 
    ArrayList<ContigData> continousPloidyContigsCluster2=new ArrayList<ContigData> () ;
	 public ArrayList<String> continousPloidyContigsNamesCluster1;
	 public ArrayList<String> continousPloidyContigsNamesCluster2 ;
	 public int LENGTH_OF_CLUSTER_ONE_CONTIGS;
	 public int LENGTH_OF_CLUSTER_TWO_CONTIGS;	
	
	
	public NaivePedcaPlotter(Map<String,ContigData> contList,int maxWindows, float[] rc, SamParser samP) {
		
		//reset static variables
		clusterMus=null;
		npdf=null;
		maxX=0;
		totalDataPoints=0;
		rt=null;
		LENGTH_OF_CLUSTER_ONE_CONTIGS=0;
		LENGTH_OF_CLUSTER_TWO_CONTIGS=0;
		continousPloidyContigsNamesCluster1=new ArrayList<String> () ; 
		continousPloidyContigsNamesCluster2=new ArrayList<String> () ; 
		
		readCounts=rc;
		contigsList=contList;
		contArrList = new ArrayList<String>(contigsList.keySet());
		sumOfReadCounts=samP.barchart.sumOfReadCounts;
		
		try{
			createFitterDataset() ;
			fitNaiveMixtureModel();	//approximate the distribution by naive smoother and infere max points				
			displayPloidyAndCoveragePlotNaive(rt.writer);//Plot containng both the coverage and the ploidy estimation	
		}catch (Exception e){
			System.err.println("Error in PloestPlotter constructor");
		}
		
	}
	
	private  void createFitterDataset() {
		
		intFitPoints=SamParser.intFitPoints;

	}
	public void naivePloestPlotter2ndRound(Map<String,ContigData> contList,int maxWindows, float[] rc) {
		readCounts=rc;
		contigsList=contList;
		contArrList = new ArrayList<String>(contigsList.keySet());
		unsolvedPloidyContigs=new ArrayList<ContigData> () ;//reinitialize unsolvedPloidyContigs
		removedContigs=new ArrayList<ContigData> ();//reinitialize removedContigs
		try{
			rt.writeOut2ndRound();
			displayPloidyAndCoveragePlotNaive(rt.writer2ndRun);//Plot containng both the coverage and the ploidy estimation
		}catch (Exception e){
			System.err.println("Error in PloestPlotter constructor");
		}
		
	}
	
	private void writeOutPloEstByFragment(PrintWriter writer , XYSeries series, ContigData contigD ) {
		
		String contigname =contigD.contigName;
		

		Number prevPloidy=0;
		Number newPloidy=0;
		int ItemsSize=0;
		int winFrac=4;//fraction by which divide the window length to compute the real position
		int firstPloidyPos=0;
		
		//get first valid ploidy point
		if(series.getItemCount()>0){
			ItemsSize=series.getItems().size();
			prevPloidy=series.getY(firstPloidyPos);//series.getMinY();
			System.out.println(contigname+" ItemsSize:"+ItemsSize+" contig size:"+contigD.maxLength+" PEDCA ws:"+Pedca.windowLength+ " SamParser.windowLength:"+SamParser.windowLength);

			while(prevPloidy.intValue()==0 && firstPloidyPos<ItemsSize-1){
				prevPloidy=series.getY(++firstPloidyPos);
			}
			newPloidy=prevPloidy;
		}
		
		
		if(series.getItemCount()>0  && prevPloidy.intValue()!=0){
			Number prevPos=0;


			for (int yv=firstPloidyPos;yv<ItemsSize;yv++){
				newPloidy=series.getY(yv);
				if(newPloidy==null)newPloidy=0;
				if(!newPloidy.equals(prevPloidy)  ){//segmentation point
					contigD.thisContigHasContinousPloidy=false;
					if(/*prevPloidy!=null && */prevPloidy.intValue()!=0){
						
						writer.println(contigname+"\ts "+((int)prevPos*SamParser.windowLength/winFrac)+"\t"+(yv*SamParser.windowLength/winFrac)+"\t"+prevPloidy);
						prevPloidy=newPloidy;						
						prevPos=yv;
					}else{
						contigD.thisContigHasContinousPloidy=false;
						prevPloidy=newPloidy;
						prevPos=yv;
					}
				}
			}
			
			if (contigD.thisContigHasContinousPloidy && (prevPloidy.intValue()!=0) ){//coninuous contig with valid ploidy
				writer.println(contigname+"\tc "+((int)prevPos*SamParser.windowLength/winFrac)+"\t"+contigD.maxLength+"\t"+prevPloidy);
			}else{//fragmented ploidy... writeout last fragment 
				
				if(!series.getY(ItemsSize-1).equals(0)){//last frag has valid ploidy
					writer.println(contigname+"\tf "+((int)prevPos*SamParser.windowLength/winFrac)+"\t"+(ItemsSize*SamParser.windowLength/winFrac)+"\t"+prevPloidy);

				}else 	if(series.getY(ItemsSize-1).equals(0)){//last frag has invalid ploidy
						writer.println(contigname+"\tf "+((int)prevPos*SamParser.windowLength/winFrac)+"\t"+(ItemsSize*SamParser.windowLength/winFrac)+"\t"+prevPloidy);
			
				}else {
					writer.println(contigname+"\tf "+((int)prevPos*SamParser.windowLength/winFrac)+"\t"+(ItemsSize*SamParser.windowLength/winFrac)+"\t"+prevPloidy);
				}
			}
			//store all contigs with ploidy belonging to cluster 1 or 2
			if(Pedca.baseCallIsOn && contigD.thisContigHasContinousPloidy ){
				
				if( prevPloidy.intValue()==rt.bestScore.bestCNVIndexes[0]){//contigs from cluster 1
					continousPloidyContigsCluster1.add(contigD);
					continousPloidyContigsNamesCluster1.add(contigD.contigName);
					LENGTH_OF_CLUSTER_ONE_CONTIGS+=contigD.maxLength;
				}else if(prevPloidy.intValue()==rt.bestScore.bestCNVIndexes[1]){//contigs from cluster 2
					continousPloidyContigsCluster2.add(contigD);
					continousPloidyContigsNamesCluster2.add(contigD.contigName);
					LENGTH_OF_CLUSTER_TWO_CONTIGS+=contigD.maxLength;
				}
			}
		
		}else{
			writer.println(contigname+"\t-\t-\t-");
		}
	}

	//Smoothing the data by averaging values in bins
	public void fitNaiveMixtureModel(){

		npdf=new NaivePDF(readCounts);
		significantMaxsInPDF(npdf);
		SamParser.barchart.BarChartWithFit(npdf,"FINALRESULT");
		rt=new RatioFindNaive(clusterMus);
		rt.writeOut();
	}
	
	//finds the index of the most represented result in this vector
	public int indexOfMode(int[] vector){
		
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
						("Genome Coverage "+contigD.contigName+ SamParser.stringSecondRound), // chart title
						"Genome Position (x " +(SamParser.windowLength/2)+" bp)", // x axis label
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
	
				if(contigD.maxY>0){
					rangeAxis.setRange(0.00,contigD.maxY);
				}else {
					rangeAxis.setRange(0.00,10);
					System.err.println(contigD.contigName+" doesn't have any coverage. Contig is Removed");
					contigsList.remove(contArrList.get(c));
					contArrList.remove(c);
					c--;
				}
				String correctedContigName = contigD.contigName.replaceAll("[^a-zA-Z0-9.-]", "_");
				ChartUtilities.saveChartAsJPEG(new File(Pedca.outputFile + "//" + Pedca.projectName+ "//Contig_Coverage_Charts//Chart_Contig_"+correctedContigName+".jpg"), chart, 1000, 600);
			}else {
				System.err.println(contigD.contigName+" is too short ("+contigD.maxLength+"bp ) for the length of the sliding window (winLength="+contigD.windLength+" bp). Contig is Removed");
				writeoutStringLog+= contigD.contigName+" is too short ("+contigD.maxLength+"bp ) for the length of the sliding window (winLength="+contigD.windLength+" bp). Contig is Removed\r\n";
			}
		}
	}


	

public void displayPloidyAndCoveragePlotNaive( PrintWriter writ)throws IOException{	

		ContigData contigD;
		
		                                                                      
		writ.println("#**************************************************************************************");
		writ.println("#*Ploidy estimation detailed by fragments.  Precision: +/-"+(PLOIDY_SMOOTHER_WIDTH*SamParser.windowLength/(2*MAX_PLOIDY))+" bp                      *");
		writ.println("#*Only detects fragments with continuous ploidy over a sequence >"+(CONTINUITY_POINTS*(SamParser.windowLength/2) )+" bp *");
		writ.println("#***************************************************************************************");
		writ.println("#");
		writ.println("PRECISION="+(SamParser.windowLength/2));
		writ.println("#");
		writ.println("PURC="+rt.bestScore.candidateUnit);
		writ.println("#>CONTIG NAME\tFROM\tTO\tPLOIDY ESTIMATION");
		writ.println("PLOIDY=");
		
		for (int c=0;c<contigsList.size();c++){//for each contig
			
			contigD=contigsList.get(contArrList.get(c));
			
			XYPlot xyPlot = new XYPlot();

			
			/* SETUP SCATTER */
			// Create the scatter data, renderer, and axis
			XYDataset collection1 = createPlotDataset(contigD);
			XYItemRenderer renderer1 = new XYLineAndShapeRenderer(false, true);   // Shapes only
			ValueAxis domain1 = new NumberAxis("Genome Position (bp)");
			ValueAxis rangeAxis = new NumberAxis("Coverage");

			// Set the scatter data, renderer, and axis into plot
			xyPlot.setDataset(0, collection1);
			xyPlot.setRenderer(0, renderer1);
			xyPlot.setDomainAxis(0, domain1);
			xyPlot.setRangeAxis(0, rangeAxis);
			
			//if(contigD.maxY>0 && contigD.maxY>rt.candUnit*15){rangeAxis.setRange(0.00,rt.candUnit*15);	}
			//THIS OPTION ON JUST FOR PRINTING OUT THESIS (REMOVE LATER?)	
			rangeAxis.setRange(0.00,rt.candUnit*MAX_PLOIDY);
			
			// Map the scatter to the first Domain and first Range
			xyPlot.mapDatasetToDomainAxis(0, 0);
			xyPlot.mapDatasetToRangeAxis(0, 0);
			
			// Create the line data, renderer, and axis
			XYDataset collection2 = createPloidyEstimationDatasetNaive(contigD,writ);
			
			 Font fontLabel = new Font("Dialog", Font.PLAIN, 30); 
			 Font fontTitle = new Font("Dialog", Font.PLAIN, 35); 
			 Font fontLegend = new Font("Dialog", Font.PLAIN, 24); 
			 
			 xyPlot.getDomainAxis().setLabelFont(fontLabel);
			 xyPlot.getRangeAxis().setLabelFont(fontLabel);
			 xyPlot.getDomainAxis().setTickLabelFont(fontLabel);
			 xyPlot.getRangeAxis().setTickLabelFont(fontLabel);
			

			if(!unsolvedPloidyContigs.contains(contigD)){//ploidy is solved for these contigs
				XYItemRenderer renderer2 = new XYLineAndShapeRenderer(false , true);   // Lines only
				ValueAxis domain2 = new NumberAxis("Genome Position");
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
				
				domain2.setLabel("");;
				domain2.setVisible(false);
				range2.setLabelFont(fontLabel);				
				range2.setTickLabelFont(fontLabel);
				
				// Create the chart with the plot and a legend of COVERAGE AND PLOIDY		
				JFreeChart chart = new JFreeChart("Coverage and Ploidy Estimation :"+contigD.contigName, JFreeChart.DEFAULT_TITLE_FONT, xyPlot, true);
				chart.getTitle().setFont(fontTitle);
				chart.getLegend().setItemFont(fontLegend);
				
				String correctedContigName = contigD.contigName.replaceAll("[^a-zA-Z0-9.-]", "_");
				ChartUtilities.saveChartAsJPEG(new File(Pedca.outputFile + "//" + Pedca.projectName+ "//Ploidy_Estimation_Charts//Ploidy_Estimation_"+correctedContigName+"_"+SamParser.stringSecondRound+".jpg"),chart, 1500, 900);
				
			}
			
			if(unsolvedPloidyContigs.contains(contigD) && SamParser.stringSecondRound=="_2nd_Round_"){//ploidy is NOT solved for these contigs AND it is the SECOND RUN

				// Create the chart with the plot and a legend OF ONLY THE COVERAGE!!			
				JFreeChart chart = new JFreeChart("Coverage and Ploidy Estimation :"+contigD.contigName, JFreeChart.DEFAULT_TITLE_FONT, xyPlot, true);
				String correctedContigName = contigD.contigName.replaceAll("[^a-zA-Z0-9.-]", "_");
				ChartUtilities.saveChartAsJPEG(new File(Pedca.outputFile + "//" + Pedca.projectName+ "//Ploidy_Estimation_Charts//Ploidy_Estimation_"+correctedContigName+"_"+SamParser.stringSecondRound+".jpg"),chart, 1500, 900);
			}
		}
			
		if(Pedca.baseCallIsOn /*&& rt.bestScore.score>0.07*/){
			
			try {
				runBaseCallCheck();
				Pedca.baseCallIsOn =false;
			} catch (InterruptedException e) {
				System.err.println("runBaseCall error exception in NaivePloestPlotter.displayPloidyAndCoveragePlotNaive()");
				e.printStackTrace();
			}
		}
		

		//if some contigs CAN NOT BE SOLVED --> remove them
		if (removedContigs.size()>0 && unsolvedPloidyContigs.size()==0 ){ //if this is the first and only round
			printRemovedContigs( writ);
		}else if (removedContigs.size()>0 && SamParser.RUN_SECOND_ROUND ){//or if this is the second round
			printRemovedContigs( writ);
		}

		//if a second run is needed... (there'are unsolved contigs)
		if (unsolvedPloidyContigs.size()>0 && !SamParser.RUN_SECOND_ROUND){

			int minLength=unsolvedPloidyContigs.get(0).maxLength;
			writ.println("#");
			writ.println("#> The following contigs could not be solved with the current window length of "+Pedca.windowLength+" bp : ");
			writ.println("CONTIGS_TO_BE_SOLVED=");
			//compute new window length
			for (int u=0;u<unsolvedPloidyContigs.size();u++){
				writ.println("  "+unsolvedPloidyContigs.get(u).contigName+ " length: "+unsolvedPloidyContigs.get(u).maxLength+" bp");
				if (minLength>=unsolvedPloidyContigs.get(u).maxLength){
					minLength=unsolvedPloidyContigs.get(u).maxLength;
				}
			}
			
			Pedca.windowLength=(int) (minLength/(3*CONTINUITY_POINTS));
			writ.println("#>A new estimation will be attempted with a new window length of "+Pedca.windowLength);
			if(Pedca.windowLength<Pedca.MIN_WIND_LENGTH)Pedca.windowLength=Pedca.MIN_WIND_LENGTH;//if(Ploest.windowLength<8)Ploest.windowLength=8;
			SamParser.RUN_SECOND_ROUND=true;	
			//update contigs info
			for (int u=0;u<unsolvedPloidyContigs.size();u++){
				unsolvedPloidyContigs.get(u).windLength=Pedca.windowLength;
			}
		}
				
		writ.close();
	}


	private void printRemovedContigs(PrintWriter writer){
		
		writer.println("#");
		writer.println("#> The following contigs were removed because of lack of data values: ");
		writer.println("REMOVED_CONTIGS=");
		for (int u=0;u<removedContigs.size();u++){
			writer.println(removedContigs.get(u).contigName);
		}
		writer.println("#");
		
}


	private void runBaseCallCheck() throws FileNotFoundException, InterruptedException {
		
		try{
			System.out.println("runBaseCall on these contigs:");
			System.out.println("cluster 1");
			for (int i =0 ;i<continousPloidyContigsCluster1.size();i++){
				System.out.println(continousPloidyContigsCluster1.get(i).contigName);
			}
			System.out.println("cluster 2");
			for (int i =0 ;i<continousPloidyContigsCluster2.size();i++){
				System.out.println(continousPloidyContigsCluster2.get(i).contigName);
			}
			VCFManager vcfManager=new VCFManager(Pedca.vcfFile.getAbsolutePath(),this);
		
		
		Pedca.baseCallIsOn=false;
		}catch (Exception e){
			
			e.printStackTrace();
		}
			
			

	
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

		double x;
		double y;
		
		//this first loop is for the PLOESTPLOTTER
		int wInd=0;//writting index
		
		for (int i = 0; i < (contigD.windPos.size()); i++) {
			x = (wInd++ *(contigD.windLength/2));  
			if(contigD.windPos.get(i)!=null){						
				y = contigD.windPos.get(i);
			}else{			
				y = 0;
			}
			
			if(y>contigD.maxY)contigD.maxY=(int)y;
			series.add(x, y);

		}

		result.addSeries(series);
		maxX=contigD.windPos.size();
		//writer.close();
		return result;
	}






	
	private  XYDataset createPloidyEstimationDatasetNaive(ContigData contigD, PrintWriter writer ) {
		XYSeriesCollection result = new XYSeriesCollection();
		XYSeries series = new XYSeries(" Ploidy Estimation 1");
		
		double [] xValues=new double[contigD.windPos.size()];
		int [] yValues=new int[contigD.windPos.size()];
		maxX=0;
		int wInd=0;
		
		//this loop is for ESTIMATED PLOIDY PLOTTING
		for (int i = 0; i < contigD.windPos.size(); i++) {

			if(contigD.windPos.get(i)!=null){
				yValues [wInd]= (int)Math.round(contigD.windPos.get(i)/RatioFindNaive.candUnit);
			}else{
				yValues [wInd]=0;
			}
			xValues [wInd]= wInd * (contigD.windLength/2);  	

			if (!AVG_PLOIDY){	
				series.add(xValues [wInd], yValues [wInd]);
			}
			wInd++;
	}
System.out.println(contigD.contigName+"1St loop Number of values:"+contigD.windPos.size()+" contig size:"+contigD.maxLength+"  series.getItems().size():"+series.getItems().size()+" series.getItemCount():"+series.getItemCount());


		if (wInd>maxX){//--wInd
			maxX=(int) wInd;
		}

		if(AVG_PLOIDY ){		//smooth the ploidy plot by averaging the values over a window of PLOIDY_SMOOTHER_WIDTH points
			series=averagePloidyMode(contigD,series,wInd,xValues,yValues);//uses the mode over PLOIDY_SMOOTHER points
		}

		if(series.getItemCount()<10 /*&& !removedContigs.contains(contigD)*/){
			unsolvedPloidyContigs.add(contigD);			
		}

		result.addSeries(series);
System.out.println(contigD.contigName+" AV      Number of values:"+contigD.windPos.size()+" contig size:"+contigD.maxLength+"  series.getItems().size():"+series.getItems().size()+" series.getItemCount():"+series.getItemCount());

		writeOutPloEstByFragment(writer, series,contigD );//writes out the ploidy estimation detailed by fragment

		return result;
	}

	
	public XYSeries averagePloidyMode(ContigData contigD,XYSeries series, int valuesSize,double [] xValues,int [] yValues){
		String contigname=contigD.contigName;

		int currentMode=0;//the most observed ploidy value over the PLOIDY_SMOOTHER_WIDTH
		int PLOIDY_SMOOTHER_WING=PLOIDY_SMOOTHER_WIDTH/2; //length of each of the sides of the PLOIDY_SMOOTHER window before and after the position being evaluated
		int [] ploidyCounter=new int[MAX_PLOIDY+1];//over the PLOIDY_SMOOTHER_WIDTH, this vector keeps track of how many times each ploidy is observed
		int modeThresholdRate=3;//the currentMode needs to have a minimal threshold
		int nbZeroValues=0;//neglects the zero values when calculating the mode
System.out.println(contigD.contigName+" AV1           valuesSize:"+valuesSize+"  series.getItems().size():"+series.getItems().size());

		if (valuesSize > PLOIDY_SMOOTHER_WIDTH) {//we need a minimum of points to average the ploidy
			
			// solve the first  positions of the plot
			for (int v = 0; v < PLOIDY_SMOOTHER_WIDTH; v++) {
				if(yValues[v]<=MAX_PLOIDY && yValues[v]>0){
					ploidyCounter[yValues[v]]++;
					if (ploidyCounter[currentMode]<ploidyCounter[yValues[v]]){
						currentMode=yValues[v];
					}
				}else if (yValues[v]==0)nbZeroValues++;
			}
			//now store the mode over the first positions
			for (int v = 0; v < PLOIDY_SMOOTHER_WING; v++) {
				if (currentMode!=0 && ploidyCounter[currentMode]>((PLOIDY_SMOOTHER_WIDTH-(nbZeroValues/8))/modeThresholdRate)){
					series.add(xValues[v], currentMode);
				}
			}
System.out.println(contigD.contigName+" AV2           valSize-wing:"+(valuesSize-PLOIDY_SMOOTHER_WING)+"  series.getItems().size():"+series.getItems().size());
			
			//solve the core of the genome until the last positions-PLOIDY_SMOOTHER_WIDTH/2
			int mostRightValue;//value at the right end of the ploidy-smoother window
			int mostLeftValue=0;//value at the left end of the ploidy-smoother window
			
			for(int v=PLOIDY_SMOOTHER_WING;v<(valuesSize-PLOIDY_SMOOTHER_WING);v++){
				
				if(yValues[v]<=MAX_PLOIDY ){
					mostLeftValue=yValues[v-PLOIDY_SMOOTHER_WING];
					mostRightValue=yValues[v+PLOIDY_SMOOTHER_WING];
					
					if(mostLeftValue!=mostRightValue){//if the removed and the added are different, update the mode
						if(mostLeftValue>0 && mostLeftValue<=MAX_PLOIDY ){	
							ploidyCounter[mostLeftValue]--;	//remove mostleft value of window
							}else if (mostLeftValue==0)nbZeroValues--;
						
						if (mostRightValue>0 && mostRightValue<=MAX_PLOIDY ){
							ploidyCounter[mostRightValue]++;//add mostright value of window
						}else if (mostRightValue==0)nbZeroValues++;
						
						//update the mode
						for(int pc=1;pc<ploidyCounter.length;pc++){
							if(ploidyCounter[pc]>ploidyCounter[currentMode])currentMode=pc;
						}//else the mode is the same

					}
					

					if (currentMode!=0 && ploidyCounter[currentMode]>((PLOIDY_SMOOTHER_WIDTH-(nbZeroValues/8))/modeThresholdRate)){
						series.add(xValues[v], currentMode);
					}
				}
			}
System.out.println(contigD.contigName+" AV3           valuesSize:"+valuesSize+"  series.getItems().size():"+series.getItems().size());

			//solve the very last positions PLOIDY_SMOOTHER_WIDTH/2
			for(int v=(valuesSize-(PLOIDY_SMOOTHER_WIDTH/2));v<valuesSize;v++){
				if(yValues[v]<=MAX_PLOIDY ){
					mostLeftValue=yValues[v-(PLOIDY_SMOOTHER_WIDTH/2)];

					if(mostLeftValue>0 && mostLeftValue<=MAX_PLOIDY && ploidyCounter[mostLeftValue]>0 ){	
						ploidyCounter[mostLeftValue]--;	//remove mostleft value of window
					}else if (mostLeftValue==0)nbZeroValues--;

					if (ploidyCounter[currentMode]<ploidyCounter[yValues[v]]){
						for(int pc=1;pc<ploidyCounter.length;pc++){//we don't consider 0 values
							if(ploidyCounter[pc]>ploidyCounter[currentMode])currentMode=pc;
						}
					}

					if (currentMode!=0 && ploidyCounter[currentMode]>((PLOIDY_SMOOTHER_WIDTH-(nbZeroValues/8))/modeThresholdRate) ){
						series.add(xValues[v], currentMode);
					}
				}
				
			}
		}else{//not enough points, simply average over the available points
			for (int v = 0; v < valuesSize; v++) {
				if(yValues[v]<=MAX_PLOIDY ){
					ploidyCounter[yValues[v]]++;
					if (ploidyCounter[currentMode]<ploidyCounter[yValues[v]]){
						currentMode=yValues[v];
					}
					if (currentMode!=0 && ploidyCounter[currentMode]>((PLOIDY_SMOOTHER_WIDTH-(nbZeroValues/8))/modeThresholdRate)){
						series.add(xValues[v], currentMode);
					}
				}
			}
		}
System.out.println(contigD.contigName+" AV4           valuesSize:"+valuesSize+"  series.getItems().size():"+series.getItems().size());

		if(series.getItemCount()>0){
			return checkContinuity(series,CONTINUITY_POINTS,contigD);
		}else{
			removedContigs.add(contigD);
			System.err.println("Error in contig :"+contigname+". This series have 0 values!!!");
			return series;
		}

		
	}
	
	private XYSeries checkContinuity(XYSeries series, int continuityLength,ContigData contigD) {
		//check that at least continuityLength points have the same copy number estimation before deciding if a fragment has a certain ploidy		
		int ctr=0;
		XYSeries result=new XYSeries(" Ploidy Estimation");
		int [] yvalues=new int[series.getItemCount()];
		Number currentPloidY=series.getY(0);//current Y value being counted
		int continousCurrents=0;//nb of contiguous observed values of the current Y value

		for (int i=0;i<series.getItemCount();i++){
			if( currentPloidY.equals(series.getY(i)) && (series.getY(i).intValue()!=0)){//if PloidY value is currently being observed
				continousCurrents++;
				if (continousCurrents>continuityLength){//if the continuity length is respected
					result.add(series.getX(i),currentPloidY);//add to series
				}else if (continousCurrents==continuityLength){//the following 2 'else if' allows the first continuityLength observed values to be taken into account 
					for (int j=0;j<continuityLength;j++){//add the stored first continuityLength observed values
						result.add((series.getX(i).intValue()-continuityLength+j+1),currentPloidY);//add to series	
						}
				}
				
			}else /*if(series.getY(i).intValue()!=0)*/{//new Y values observed
				
				if (continousCurrents>=continuityLength){
					for (int j=0;j<continousCurrents;j++){//add the  first continousCurrents  values
						if (currentPloidY.intValue()!=0){
							result.add((series.getX(i-1).intValue()-continousCurrents+j+1),currentPloidY);//add to series	
						}
					}
				}
				continousCurrents=0;//reset counter
				currentPloidY=series.getY(i);//reset observed Y (ploidy)
			}
		}

		return result;
	}




	static int getFinalNumberOfMixtures(){
		return finalNumberOfMixtures;
	}

	

	private int significantMaxsInPDF(NaivePDF naivePDF) {
		int sigMaxs = 0;//nb of significant maximums
		//double ReadCountThreshold = Pedca.SIGNIFICANT_MIN;// discards the values that are below this threshold in the original readCount Distribution
		//int FITthreshold
		System.out.println(" Sum of all ReadCounts:"+sumOfReadCounts+" Pedca.SIGNIFICANT_MIN:"+Pedca.SIGNIFICANT_MIN+"; maxY :" + SamParser.readDistributionMaxY );

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

		while (ind < naivePDF.yDataPoints.length) {

			if (naivePDF.yDataPoints[ind] != right) {
				//move and update pointers
				left = mid;
				mid = right;
				right = naivePDF.yDataPoints[ind];
				lastLeftIndex = lastMidIndex;
				lastMidIndex = lastRightIndex;
				lastRightIndex = ind;

				if (right < mid && mid > left /*&& mid > FITthreshold*/) {// we count maxs only above threshold (deactivated here, threshold is checked below) 
					
					// now that we encountered a max bin, we scan the previous
					// and the following bins to find the exact maximum point in this area
					
					double maxVal = 0;// precise max Y value in the corresponding bins
					int Xindex = lastLeftIndex ;//index (x value) of the maxVal
					int endscaningIndex=lastRightIndex + naivePDF.smootherLength/2;
					if(endscaningIndex>readCounts.length-1)endscaningIndex=readCounts.length-1;

					for (int ib = lastLeftIndex ; ib < endscaningIndex; ib++) {
						if (readCounts[ib] > maxVal) {
							maxVal = readCounts[ib];
							Xindex = ib;
						}
					}
					
					
					if(!xMinList.contains(naivePDF.xDataPoints[Xindex]) && maxVal>Pedca.SIGNIFICANT_MIN ){// we count maxs only above threshold
						yMinList.add(maxVal);
						xMinList.add(naivePDF.xDataPoints[Xindex]);
						System.out.println(" ****    max in :" + naivePDF.xDataPoints[Xindex] + " = " + maxVal+ " lastLeftIndex:"+lastLeftIndex+" lastRightIndex:"+endscaningIndex);
						sigMaxs++;
					}else {
						System.out.print(" REPEATED max in :" + naivePDF.xDataPoints[Xindex] + " = " + maxVal+"  ...VALUE DISCARDED...!");
						if (maxVal<=Pedca.SIGNIFICANT_MIN ){
							System.out.print( " VALUE BELOW THRESHOLD");
						}System.out.println();
					
					}

					

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
