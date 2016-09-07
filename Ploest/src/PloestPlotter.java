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
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.plot.XYPlot;
import org.jfree.data.xy.XYDataset;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;

import dataFitters.*;
import jMEF.MixtureModel;
import jMEF.PVector;

public class PloestPlotter {
	static final int MAX_NB_MIXTURES=10;
	Map<String,ContigData> contigsList;

	PVector[] fitPoints;
	JFreeChart chart;
	static int maxX=0;
	static int maxY=0;
	static int totalDataPoints=0;//total number of input datapoints (coverage for all windows)
	GaussianMixturePDF[] gmPDF;//Y points of the gaussian fit


	public PloestPlotter(Map<String,ContigData> contList,int maxWindows) {
		System.out.println("maxWindows:"+maxWindows);
		contigsList=contList;
		maxY= (int) (maxWindows*5);
		try{
			displayScatterPlot();
			createFitterDataset() ;
			fitGaussianMixtureModel();

		}catch (Exception e){
			System.out.println("Error in PloestPlotter constructor");
			System.err.println("Error in PloestPlotter constructor");
		}
	}

	public void fitGaussianMixtureModel(){
		System.out.println("-------- -----fitGaussianMixtureModel-------- ----------");
		double aic;
		double bic;
		double[]aicsEM=new double[MAX_NB_MIXTURES+1];//each EM (Expectation Maximization) AIC value is stored in its corresponding k(number of mixtures) index
		double[]bicsEM=new double[MAX_NB_MIXTURES+1];//same for EM BIC values
		double[]aicsBSC=new double[MAX_NB_MIXTURES+1];//each BSC (Bregman Soft Clustering)AIC value is stored in its corresponding k(number of mixtures) index
		double[]bicsBSC=new double[MAX_NB_MIXTURES+1];//same for BSC BIC values
		MixtureModel[]emMMs=new MixtureModel[MAX_NB_MIXTURES+1];//contains the result of the EM fit for each of the mixtures
		MixtureModel[]bscMMs=new MixtureModel[MAX_NB_MIXTURES+1];//contains the result of the BSC fit for each of the mixtures
		GaussianDataFitter df;
		int NbOfRuns=100;
		int indofmin;
		int[] bestGuess=new int[MAX_NB_MIXTURES];//best guess with bic model evaluation
		int[] correctedResults=new int[MAX_NB_MIXTURES];//best guess with bic model evaluation + min points evaluation
		gmPDF=new GaussianMixturePDF[NbOfRuns];
		for (int r=0;r<NbOfRuns;r++){
			System.out.println("-------- -----fitGaussianMixtureModel-------- NbOfRuns:"+r+"----------");
			for (int k=1;k<MAX_NB_MIXTURES;k++){//fit to different number k of mixtures
				df=new GaussianDataFitter (fitPoints,k );
				
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
			/*PRINT RESULTS TO CONSOLE
			for (int k=1;k<MAX_NB_MIXTURES;k++){//printout results


				if(!Double.isNaN(aicsEM[k])){
					System.out.println(" k:"+k+"    EM aic:"+aicsEM[k]+" EM bic:"+bicsEM[k]);
				}

				if(!Double.isNaN(bicsBSC[k])){
					System.out.println(" k:"+k+"    BSC aic:"+aicsBSC[k]+" BS bic:"+bicsBSC[k]);
				}

			}
			 */
			System.out.println("--------Run:"+r+" ----------------");
			indofmin=findIndexOfMin(bicsBSC);//find index of gaussian mixture with min BIC value
			bestGuess[indofmin]++;
			System.out.println(indofmin+" MM  BSC params:"+bscMMs[indofmin].printParams());	
			gmPDF[r]=new GaussianMixturePDF(bscMMs[indofmin],0.0,(double)SamParser.readCounts.length,0.1);
	
			//significantMinsInPDF( gmPDF[r]);

			
			int k=significantMinsInPDF( gmPDF[r]);
			System.out.println("------RESULT:"+k+" mixtures ----------");
			df=new GaussianDataFitter (fitPoints,k );
			System.out.println(k+" MM   params:"+	df.getBSCModel().printParams());
			correctedResults[k]++;
			System.out.println("-------- ------------------ ----------");
			System.out.println("-------- ------------------ ----------");
			//*/
		}
		System.out.println("-------- -----best Guess-------- ----------");
		System.out.print("[");
		for (int b=0;b<bestGuess.length;b++){
			System.out.print(" "+bestGuess[b]);
		}
		System.out.println(" ];");
		System.out.println("-------- -----best correctedResults-------- ----------");
		System.out.print("[");
		for (int b=0;b<correctedResults.length;b++){
			System.out.print(" "+correctedResults[b]);
		}
		System.out.println(" ];");
	}

	public void displayScatterPlot() throws IOException{
		List<String> contArrList = new ArrayList<String>(contigsList.keySet());
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
			//System.out.println("setRange maxx="+maxX+ " maxy="+maxY);
			rangeAxis.setRange(0.00, maxY);

			// create and display a frame...
			/*
		        ChartFrame frame = new ChartFrame("Coverage for "+contigD.contigName, chart);
		        frame.pack();
		        frame.setVisible(true);	
			 */
			ChartUtilities.saveChartAsJPEG(new File(Ploest.outputFile + "//" + Ploest.projectName+ "//Chart_Contig_"+contigD.contigName+".jpg"), chart, maxX, maxY);

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
		int sumOfValues=0;
		int ind=0;
		List<String> contArrList = new ArrayList<String>(contigsList.keySet());

		for (int c=0;c<contigsList.size();c++){//for each contig
			ContigData contigD=contigsList.get(contArrList.get(c));
			for (int i = 0; i < (contigD.windPos.length); i++) {//for each window position
				PVector curVec=new PVector(1);
				curVec.array[0]=contigD.windPos[i];//coverage per window position
				sumOfValues+=contigD.windPos[i];//add to the total of all values (used by cleaning marginal values)
				fitPoints[ind++]=curVec;               
			}
		}
		/*CLEAN MARGINAL VALUES
		System.out.println("-------- -----CLEAN MARGINAL VALUES-------- ----------");
		for (int fp=0;fp<fitPoints.length;fp++){
			if ((fitPoints[fp].array[0]/sumOfValues)<0.005){
				System.out.println("fitPoint:"+fitPoints[fp].array[0]+"/"+sumOfValues);
				fitPoints[fp].array[0]=0.0;
			}
		}
		System.out.println("-------- -----CLEAN MARGINAL fin-------- ----------");
*/
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
			left = mid;
			mid = right;
			right = gmf.yDataPoints[ind];

			if (right > mid && mid < left && mid<threshold){
				yMinList.add(mid);
				xMinList.add(gmf.xDataPoints[ind]);
				sigMins++;
			}

			ind++;
		}		
		if(sigMins>3){
			System.out.println("Mins List x:"+xMinList);
			System.out.println("Mins List y:"+yMinList);
		}
		return sigMins;
	}

	
	
	
	
	
	/*
	public MixtureModel bestMixtureModel(){
		for(int i=1;i<MAX_NB_MIXTURES;i++){

		}
	}
	 */

}
