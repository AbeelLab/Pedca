import java.awt.Color;
import java.awt.Rectangle;
import java.awt.Shape;
import java.awt.geom.Rectangle2D;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;

import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartFrame;
import org.jfree.chart.ChartUtilities;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.axis.NumberAxis;
import org.jfree.chart.axis.ValueAxis;
import org.jfree.chart.plot.DatasetRenderingOrder;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.renderer.xy.StandardXYItemRenderer;
import org.jfree.chart.renderer.xy.XYItemRenderer;
import org.jfree.data.xy.XYDataset;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;

import dataFitters.*;
import jMEF.MixtureModel;
import jMEF.PVector;

public class PloestPlotter {
	static final int MAX_NB_MIXTURES=10;
	Map<String,ContigData> contigsList;
	List<String> contArrList;
	static double[] gMMweights;//result of the WEIGHTS of the gaussians in the mixture model fitting
	static double[] gMMmus;//result of the MEANS (mus) of the gaussians in the mixture model fitting
	static double[] gMMsigmas;//result of the STANDARD VARIATIONS (sigmas) of the gaussians in the mixture model fitting
	static GaussianDataFitter dfResult;//final data fit with the best fit predicted 
	PVector[] fitPoints;
	JFreeChart chart;
	static int maxX=0;
	static int maxY=0;
	static int totalDataPoints=0;//total number of input datapoints (coverage for all windows)
	GaussianMixturePDF[] gmPDF;//Y points of the gaussian fit
	static int finalNumberOfMixtures;
	RatioFind rt;//contains the ratio of the gaussians to the ploidy unit which allows computation of ploidy from contig coverage

	public PloestPlotter(Map<String,ContigData> contList,int maxWindows) {
		contigsList=contList;
		contArrList = new ArrayList<String>(contigsList.keySet());
		maxY= (int) (maxWindows*5);//Y axis range
		try{
			displayScatterPlot();
			//findTotalDataPoints();
			createFitterDataset() ;
			fitGaussianMixtureModel();
			displayPloidyEstimationScatterPlot();

		}catch (Exception e){
			System.out.println("Error in PloestPlotter constructor");
			System.err.println("Error in PloestPlotter constructor");
		}
	}
/*
	private void findTotalDataPoints() {

		for (int c=0;c<contigsList.size();c++){//for each contig
			totalDataPoints+=contigsList.get(contArrList.get(c)).windPos.length+1;
		}
		System.out.println("totalDataPoints:"+totalDataPoints);
	}
*/
	public void fitGaussianMixtureModel(){
		System.out.println("-------------Fitting Different Gaussian Mixture Models------------------");
		double aic;
		double bic;
		double[]aicsEM=new double[MAX_NB_MIXTURES+1];//each EM (Expectation Maximization) AIC value is stored in its corresponding k(number of mixtures) index
		double[]bicsEM=new double[MAX_NB_MIXTURES+1];//same for EM BIC values
		double[]aicsBSC=new double[MAX_NB_MIXTURES+1];//each BSC (Bregman Soft Clustering)AIC value is stored in its corresponding k(number of mixtures) index
		double[]bicsBSC=new double[MAX_NB_MIXTURES+1];//same for BSC BIC values
		MixtureModel[]emMMs=new MixtureModel[MAX_NB_MIXTURES+1];//contains the result of the EM fit for each of the mixtures
		MixtureModel[]bscMMs=new MixtureModel[MAX_NB_MIXTURES+1];//contains the result of the BSC fit for each of the mixtures
		GaussianDataFitter df;
		int NbOfRuns=10;
		int indofmin;
		int[] bestBSCGuess=new int[MAX_NB_MIXTURES];//best BSC guess with bic model evaluation
		int[] bestEMGuess=new int[MAX_NB_MIXTURES];//best BSC guess with EM model evaluation
		int[] correctedResults=new int[MAX_NB_MIXTURES];//best guess with bic model evaluation + min points evaluation
		gmPDF=new GaussianMixturePDF[NbOfRuns];
		for (int r=0;r<NbOfRuns;r++){
		
			for (int k=1;k<MAX_NB_MIXTURES;k++){//fit to different number k of mixtures
				df=new GaussianDataFitter (fitPoints,k );//fits a gauss mixture(by EM and BSC) to the 
				emMMs[k]=df.getEMmodel();//store EM fit values
				bscMMs[k]=df.getBSCModel();//same for BSC
				//get EM LogLikelihoods and estimate BIC and AIC values
				aic=-2*(df.getEMLogLikelihood())+(2*(k*3));//k+3 are the free parameters: k=number of models +3(weight,mu and sigma)
				bic=-0.5*(df.getEMLogLikelihood())+((k*3)*Math.log(totalDataPoints));
				aicsEM[k]=aic;
				bicsEM[k]=bic;
				//get BSC LogLikelihoods and estimate BIC and AIC values
				aic=-2*(df.getBSCLogLikelihood())+(2*(k*3));
				bic=-0.5*(df.getBSCLogLikelihood())+((k*3)*Math.log(totalDataPoints));
				aicsBSC[k]=aic;
				bicsBSC[k]=bic;	
			}
		
			indofmin=findIndexOfMin(bicsBSC);//find index of gaussian mixture with min BIC value
			bestBSCGuess[indofmin]++;//System.out.println(indofmin+" MM  BSC params:"+bscMMs[indofmin].printParams());	
			indofmin=findIndexOfMin(bicsEM);//find index of gaussian mixture with min EM value
			bestEMGuess[indofmin]++;
			gmPDF[r]=new GaussianMixturePDF(bscMMs[indofmin],0.0,(double)SamParser.readCounts.length,0.1);
			int k=significantMinsInPDF( gmPDF[r]);//finds the significant minimums (below a threshold) in the gaussian mixture PDF
												  //helps determine the real number of gaussians in the PDF and correct overfitting
			correctedResults[k]++;
		}
		//PRINT BEST GUESSES AND CORRECTED RESULTS		
		System.out.println("----------BSC and EM prediction-------- ----------");//identical results with best guess EM
		System.out.print("[");
		for (int b=0;b<bestBSCGuess.length;b++){
			System.out.print(" "+bestBSCGuess[b]);
		}
		System.out.println(" ];");
		
		System.out.println("------- best correctedResults (after minimal points detection of PDF)----");
		System.out.print("[");
		for (int b=0;b<correctedResults.length;b++){
			System.out.print(" "+correctedResults[b]);
		}
		System.out.println(" ];");
		//*/
		finalNumberOfMixtures=indexOfMode(correctedResults);
		System.out.println("-------------FINAL RESULT :"+finalNumberOfMixtures+" GAUSSIAN MODELS WITH "+correctedResults[finalNumberOfMixtures]+" % consensus---");
		
		dfResult=new GaussianDataFitter (fitPoints,finalNumberOfMixtures );//final data fit with the best number of mixture prediction
		
		gMMweights=new double[dfResult.getBSCModel().weight.length];
		gMMmus=new double[dfResult.getBSCModel().param.length];
		gMMsigmas=new double[dfResult.getBSCModel().param.length];
		for (int g=0;g<dfResult.getBSCModel().param.length;g++){
			gMMweights[g]=dfResult.getBSCModel().weight[g];
			gMMmus[g]=((jMEF.PVector)dfResult.getBSCModel().param[g]).array[0];
			gMMsigmas[g]=((jMEF.PVector)dfResult.getBSCModel().param[g]).array[1];			
		}
		SamParser.barchart.BarChartWithFit(new GaussianMixturePDF(dfResult.getBSCModel(),0.0,(double)SamParser.readCounts.length,0.1),finalNumberOfMixtures);
		
		rt=new RatioFind(gMMmus,correctedResults[finalNumberOfMixtures]);
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
	
	public void displayScatterPlot() throws IOException{

		//contArrList = new ArrayList<String>(contigsList.keySet());
		for (int c=0;c<contigsList.size();c++){//for each contig
			ContigData contigD=contigsList.get(contArrList.get(c));
			
			chart = ChartFactory.createScatterPlot(
					("Genome Coverage "+contigD.contigName), // chart title
					"Genome Position (x 1000 bp)", // x axis label
					"Coverage", // y axis label
					createPlotDataset(contigD), // XYDataset 
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
			rangeAxis.setRange(0.00, maxY);
	
			ChartUtilities.saveChartAsJPEG(new File(Ploest.outputFile + "//" + Ploest.projectName+ "//Contig_Coverage_Charts//Chart_Contig_"+contigD.contigName+".jpg"), chart, 1000, 600);

		}
	}

	public void displayPloidyEstimationScatterPlot() throws IOException{
		List<String> contArrList = new ArrayList<String>(contigsList.keySet());
		for (int c=0;c<contigsList.size();c++){//for each contig
			ContigData contigD=contigsList.get(contArrList.get(c));
			XYDataset data = createPloidyEstimationDataset(contigD);
			
			chart = ChartFactory.createScatterPlot(
					("Ploidy Estimation "+contigD.contigName), // chart title
					"Genome Position (x 1000 bp)", // x axis label
					"Coverage", // y axis label
					createPloidyEstimationDataset(contigD) , // XYDataset 
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
			rangeAxis.setRange(0.00, MAX_NB_MIXTURES);	
			
			XYItemRenderer renderer = new StandardXYItemRenderer();

			xyPlot.setDataset(0, data);
			xyPlot.setRenderer(0, renderer);
			renderer.setSeriesPaint(0, Color.blue);
		    double size = 20.0;
		    double delta = size / 2.0;
		    Shape shape1 = new Rectangle2D.Double(-delta, -delta, size, size);
			renderer.setSeriesShape(0,shape1);
			//JFreeChart overlaidChart=new JFreeChart("Ploidy Estimation", JFreeChart.DEFAULT_TITLE_FONT, xyPlot, true);

			ChartUtilities.saveChartAsJPEG(new File(Ploest.outputFile + "//" + Ploest.projectName+ "//Ploidy_Estimation_Charts//Ploidy_Estimation_Contig_"+contigD.contigName+".jpg"), chart, 1000, 600);

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

	private static XYDataset createPlotDataset(ContigData contigD) {

		XYSeriesCollection result = new XYSeriesCollection();
		XYSeries series = new XYSeries(contigD.getContigName());
		maxX=0;
		double x;
		double y;
		//this first loop is for the PLOESTPLOTTER
		for (int i = 0; i <= (contigD.windPos.length-1); i++) {
			x = i;        
			y = contigD.windPos[i];
			series.add(x, y);
			if (x>maxX)maxX=(int) x;           
		}
		result.addSeries(series);
		totalDataPoints+=maxX+1;
		return result;
	}


	private  void createFitterDataset() {
		fitPoints=new PVector[totalDataPoints];
		int ind=0;
		//contArrList = new ArrayList<String>(contigsList.keySet());
		for (int c=0;c<contigsList.size();c++){//for each contig
			ContigData contigD=contigsList.get(contArrList.get(c));
			for (int i = 0; i < (contigD.windPos.length); i++) {//for each window position
				PVector curVec=new PVector(1);
				curVec.array[0]=contigD.windPos[i];//coverage per window position
				fitPoints[ind++]=curVec;               
			}
		}
	}

	
	private  XYDataset createPloidyEstimationDataset(ContigData contigD) {
		XYSeriesCollection result = new XYSeriesCollection();
		XYSeries series = new XYSeries(contigD.getContigName());

		double x;
		double y;
		//this loop is for ESTIMATED PLOIDY PLOTTING
		for (int i = 0; i < contigD.windPos.length; i++) {
			x = i;        
			y = (double)((int)(contigD.windPos[i]/rt.bestScore.getCandidateUnit()));
			series.add(x, y);		         
		}
		result.addSeries(series);
		return result;
	}
	
	
	class MaxMinArrays{
		ArrayList<Double> maxList ;		
		ArrayList<Double> minList ;
		public MaxMinArrays(){
			maxList = new ArrayList<Double>();
			minList = new ArrayList<Double>();
		}
		public void addToMinList(Double d) {
			this.minList.add(d);
		}
		public void addToMaxList(Double d) {
			this.maxList.add(d);
		}
		public ArrayList<Double> getMaxList() {
			return maxList;
		}
		public ArrayList<Double> getMinList() {
			return minList;
		}
		public int getNumberOfSignificantMins(){//returns the number of mins in the function that are below a certain threshold
			int sigMins=0;
			double threshold=0.01;
			for (int i=0;i<minList.size();i++){
				if (minList.get(i) < threshold){
					sigMins++;
				}
			}
			return sigMins;
		}
		
	}
	
	static int getFinalNumberOfMixtures(){
		return finalNumberOfMixtures;
	}
	
	private int significantMinsInPDF(GaussianMixturePDF gmf) {
		//MaxMinArrays lists = new MaxMinArrays();
		int sigMins=0;
		double threshold=0.005;//discards the mins that are above this threshold 
		//First find Mins:
		int count = 0; // To handle special case of singleton list
		ArrayList<Double> yMinList= new  ArrayList<Double>();
		ArrayList<Double> xMinList= new  ArrayList<Double>();
		double left  = gmf.maxYvalue;
		double mid   = gmf.maxYvalue;
		double right = gmf.maxYvalue;
		int ind=0;
		while (ind<gmf.yDataPoints.length) {
			count++;
			if(gmf.yDataPoints[ind]!=right){
				left = mid;
				mid = right;
			}	
				right = gmf.yDataPoints[ind];

				if (ind!=1 && right > mid && mid < left && mid < threshold) {// we count mins only after index 2 and above threshold
					yMinList.add(mid);
					xMinList.add(gmf.xDataPoints[ind]);					
					sigMins++;
				}
			
			ind++;
		}
		//add last point as min
		yMinList.add(mid);
		xMinList.add(gmf.xDataPoints[ind-1]);					
		sigMins++;
		
		return sigMins;
	}

	
	
	
	
	
	/*
	public MixtureModel bestMixtureModel(){
		for(int i=1;i<MAX_NB_MIXTURES;i++){

		}
	}
	 */

}
