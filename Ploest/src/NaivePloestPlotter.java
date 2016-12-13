import java.awt.Color;
import java.awt.Shape;
import java.awt.geom.Rectangle2D;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
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

	static float[] readCounts;
	static final int MAX_NB_MIXTURES=10;
	Map<String,ContigData> contigsList;
	List<String> contArrList;
	static double[] clustersWeights;//final result of the WEIGHTS of the clausters in the mixture model fitting
	static double[] clusterMus;//final result of the MEANS (mus) of the clusters in the mixture model fitting
	static double[] clusterSigmas;//final result of the STANDARD VARIATIONS (sigmas) of the clausters in the mixture model fitting
	static GaussianDataFitter dfResult;//final data fit with the best fit predicted 
//	static GaussianMixturePDF gmPDFResult;//final Gaussian Mixture PDF
	static PoissonDataFitter pdfResult;//final data fit with the best POISSON fit predicted 
	static PoissonMixturePDF pmPDFResult;//final POISSON Mixture PDF
	static NaivePDF npdf;//Naive Smoothed Density Function
	PVector[] fitPoints;
	JFreeChart chart;
	static int maxX=0;
	//static int maxY=0;
	static int totalDataPoints=0;//total number of input datapoints (coverage for all windows)

	static int finalNumberOfMixtures;
	RatioFindNaive rt;//contains the ratio of the gaussians to the ploidy unit which allows computation of ploidy from contig coverage
	
	
	public NaivePloestPlotter(Map<String,ContigData> contList,int maxWindows, float[] rc) {
		readCounts=rc;
		contigsList=contList;
		contArrList = new ArrayList<String>(contigsList.keySet());

		
		try{
			displayScatterPlot();			
			createFitterDataset() ;
			fitNaiveMixtureModel();
			
			/*
			
			displayPloidyAndCoveragePlotPoisson();
			*/
	

		}catch (Exception e){
			System.err.println("Error in PloestPlotter constructor");
		}
		
	}


	
	
	public void fitNaiveMixtureModel(){
		System.out.println("-------------Fitting Different Poisson Mixture Models------------------");
		
		npdf=new NaivePDF(readCounts);
		int k=significantMaxsInPDF(npdf);
		SamParser.barchart.BarChartWithFit(npdf,finalNumberOfMixtures,"FINALRESULT");

		rt=new RatioFindNaive(clusterMus);
		rt.writeOut();
		
		
	}
	
	

	
	
	private void sortPMMs() {//sort gMMmu in ascending order, then sort their weights and sigmas in the
		//corresponding order

		Vector<Double> origCopy = new Vector<Double>(clusterMus.length);//get a copy of the original order
		for (int o=0;o<clusterMus.length;o++){
			origCopy.add(clusterMus[o]);
		}	

		java.util.Arrays.sort(clusterMus);//sort the reference vector
		int[]order=new int[origCopy.size()];
		for (int o=0;o<order.length;o++){////store the order
			order[o]=origCopy.indexOf(clusterMus[o]);
		}
	
		//now order the otherGMM
		//double[]newgMMsigmas=new double[pdfResult.getBSCModel().param.length];//newly order sigmas
		double[]newgMMweights=new double[pdfResult.getBSCModel().param.length];//newly order weights
		for (int o=0;o<order.length;o++){
			//newgMMsigmas[o]=gMMsigmas[order[o]];
			newgMMweights[o]=clustersWeights[order[o]];
		}
		//gMMsigmas=newgMMsigmas;
		clustersWeights=newgMMweights;

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
		ContigData contigD;
		for (int c=0;c<contigsList.size();c++){//for each contig
			contigD=contigsList.get(contArrList.get(c));

			XYDataset data1=createPlotDataset(contigD);
			chart = ChartFactory.createScatterPlot(
					("Genome Coverage "+contigD.contigName), // chart title
					"Genome Position (x 1000 bp)", // x axis label
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
			rangeAxis.setRange(0.00, SamParser.maxWindows);
			ChartUtilities.saveChartAsJPEG(new File(Ploest.outputFile + "//" + Ploest.projectName+ "//Contig_Coverage_Charts//Chart_Contig_"+c+".jpg"), chart, 1000, 600);
		}

	}


	public void displayPloidyAndCoveragePlot()throws IOException{
		

		ContigData contigD;

		for (int c=0;c<contigsList.size();c++){//for each contig
			System.out.println("  displayPloidyAndCoveragePlot contig:" +c);
			contigD=contigsList.get(contArrList.get(c));
			XYPlot xyPlot = new XYPlot();

			/* SETUP SCATTER */

			// Create the scatter data, renderer, and axis
			XYDataset collection1 = createPlotDataset(contigD);
			XYItemRenderer renderer1 = new XYLineAndShapeRenderer(false, true);   // Shapes only
			ValueAxis domain1 = new NumberAxis("Genome Position (x 1000 bp)");
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
			XYDataset collection2 = createPloidyEstimationDataset(contigD);
			XYItemRenderer renderer2 = new XYLineAndShapeRenderer(true, false);   // Lines only
			ValueAxis domain2 = new NumberAxis("Genome Position (x 1000 bp)");
			ValueAxis range2 = new NumberAxis("Ploidy Estimation");

			// Set the line data, renderer, and axis into plot
			domain2.setStandardTickUnits(NumberAxis.createIntegerTickUnits());
			xyPlot.setDataset(1, collection2);
			xyPlot.setRenderer(1, renderer2);
			xyPlot.setDomainAxis(1, domain2);
			xyPlot.setRangeAxis(1, range2);

			// Map the line to the second Domain and second Range
			xyPlot.mapDatasetToDomainAxis(1, 1);
			xyPlot.mapDatasetToRangeAxis(1, 1);

			xyPlot.setDatasetRenderingOrder( DatasetRenderingOrder.REVERSE );
			// Create the chart with the plot and a legend
			JFreeChart chart = new JFreeChart("Coverage and Ploidy Estimation :"+contigD.contigName, JFreeChart.DEFAULT_TITLE_FONT, xyPlot, true);
			ChartUtilities.saveChartAsJPEG(new File(Ploest.outputFile + "//" + Ploest.projectName+ "//Ploidy_Estimation_Charts//Ploidy_Estimation_Contig_"+c+".jpg"),chart, 1000, 600);

		}
	}


public void displayPloidyAndCoveragePlotPoisson()throws IOException{
		

		ContigData contigD;

		for (int c=0;c<contigsList.size();c++){//for each contig

			contigD=contigsList.get(contArrList.get(c));
			XYPlot xyPlot = new XYPlot();

			/* SETUP SCATTER */

			// Create the scatter data, renderer, and axis
			XYDataset collection1 = createPlotDataset(contigD);
			
			XYItemRenderer renderer1 = new XYLineAndShapeRenderer(false, true);   // Shapes only
			ValueAxis domain1 = new NumberAxis("Genome Position (x 1000 bp)");
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
			XYDataset collection2 = createPloidyEstimationDatasetPoisson(contigD);
			XYItemRenderer renderer2 = new XYLineAndShapeRenderer(true, false);   // Lines only
			ValueAxis domain2 = new NumberAxis("Genome Position (x 1000 bp)");
			ValueAxis range2 = new NumberAxis("Ploidy Estimation");

			// Set the line data, renderer, and axis into plot
			domain2.setStandardTickUnits(NumberAxis.createIntegerTickUnits());
			xyPlot.setDataset(1, collection2);
			xyPlot.setRenderer(1, renderer2);
			xyPlot.setDomainAxis(1, domain2);
			xyPlot.setRangeAxis(1, range2);

			// Map the line to the second Domain and second Range
			xyPlot.mapDatasetToDomainAxis(1, 1);
			xyPlot.mapDatasetToRangeAxis(1, 1);

			xyPlot.setDatasetRenderingOrder( DatasetRenderingOrder.REVERSE );
			// Create the chart with the plot and a legend
			JFreeChart chart = new JFreeChart("Coverage and Ploidy Estimation :"+contigD.contigName, JFreeChart.DEFAULT_TITLE_FONT, xyPlot, true);
			ChartUtilities.saveChartAsJPEG(new File(Ploest.outputFile + "//" + Ploest.projectName+ "//Ploidy_Estimation_Charts//Ploidy_Estimation_Contig_"+c+".jpg"),chart, 1000, 600);

		}
	}

public void displayPloidyEstimationScatterPlotPoisson() throws IOException{//not used
	//XYDataset
	for (int c=0;c<contigsList.size();c++){//for each contig
		ContigData contigD=contigsList.get(contArrList.get(c));
		XYDataset data2 = createPloidyEstimationDatasetPoisson(contigD);

		chart = ChartFactory.createScatterPlot(
				("Ploidy Estimation "+contigD.contigName), // chart title
				"Genome Position (x 1000 bp)", // x axis label
				"Coverage", // y axis label
				data2 , // XYDataset 
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
		xyPlot.setRenderer(0, renderer);
		renderer.setSeriesPaint(0, Color.blue);
		double size = 20.0;
		double delta = size / 2.0;
		Shape shape1 = new Rectangle2D.Double(-delta, -delta, size, size);
		renderer.setSeriesShape(0,shape1);

		ChartUtilities.saveChartAsJPEG(new File(Ploest.outputFile + "//" + Ploest.projectName+ "//Ploidy_Estimation_Charts//Ploidy_Estimation_Contig_"+contigD.contigName+".jpg"),chart, 1000, 600);

	}
}


	public void displayPloidyEstimationScatterPlot() throws IOException{//not used
		//XYDataset
		for (int c=0;c<contigsList.size();c++){//for each contig
			ContigData contigD=contigsList.get(contArrList.get(c));
			XYDataset data2 = createPloidyEstimationDataset(contigD);

			chart = ChartFactory.createScatterPlot(
					("Ploidy Estimation "+contigD.contigName), // chart title
					"Genome Position (x 1000 bp)", // x axis label
					"Coverage", // y axis label
					data2 , // XYDataset 
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
			xyPlot.setRenderer(0, renderer);
			renderer.setSeriesPaint(0, Color.blue);
			double size = 20.0;
			double delta = size / 2.0;
			Shape shape1 = new Rectangle2D.Double(-delta, -delta, size, size);
			renderer.setSeriesShape(0,shape1);

			ChartUtilities.saveChartAsJPEG(new File(Ploest.outputFile + "//" + Ploest.projectName+ "//Ploidy_Estimation_Charts//Ploidy_Estimation_Contig_"+contigD.contigName+".jpg"),chart, 1000, 600);

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

		//PrintWriter writer = new PrintWriter(Ploest.outputFile + "//" + Ploest.projectName+ "//plotDataSet"+maxX+".txt", "UTF-8");

		double x;
		double y;
		//this first loop is for the PLOESTPLOTTER
		int wInd=0;//writting index
	
		for (int i = 0; i <= (contigD.windPos.size()-1); i++) {
			if(contigD.windPos.get(i)!=null){
				x = wInd++;  			
				y = contigD.windPos.get(i);
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


	private  XYDataset createPloidyEstimationDataset(ContigData contigD) {
		System.out.println("createPloidyEstimationDataset beg");
		XYSeriesCollection result = new XYSeriesCollection();
		XYSeries series = new XYSeries(" Ploidy Estimation");

		double x;
		double y;
		maxX=0;
		int wInd=0;
		//this loop is for ESTIMATED PLOIDY PLOTTING
		for (int i = 0; i < contigD.windPos.size(); i++) {

			if(contigD.windPos.get(i)!=null){
				x = wInd++;  			
				y = getPointPloidyEstimation(contigD.windPos.get(i));
				series.add(x, y);
				if (x>maxX)maxX=(int) x; 
				//writer.println( " x:" +x + " y:"+y);
			}				

		}
		System.out.println();
		result.addSeries(series);
		System.out.println("createPloidyEstimationDataset end");
		return result;
	}
	
	private  XYDataset createPloidyEstimationDatasetPoisson(ContigData contigD) {
		XYSeriesCollection result = new XYSeriesCollection();
		XYSeries series = new XYSeries(" Ploidy Estimation");

		double x;
		double y;
		maxX=0;
		int wInd=0;
		//this loop is for ESTIMATED PLOIDY PLOTTING
		for (int i = 0; i < contigD.windPos.size(); i++) {

			if(contigD.windPos.get(i)!=null){
				x = wInd++;  			
				y = getPointPloidyEstimationPoisson(contigD.windPos.get(i));
				series.add(x, y);
				if (x>maxX)maxX=(int) x; 
				//writer.println( " x:" +x + " y:"+y);
			}				

		}
		System.out.println();
		result.addSeries(series);

		return result;
	}
	
	public double getPointPloidyEstimation(double ptCoverage){//computes the probability of the data point 
		//belonging to all gaussians and returns the best option

		double [] pdfVals=new double[clusterMus.length];
		for (int mm=0;mm<clusterMus.length;mm++){//for all mixtures computes the probability pdfVal of the data point belonging to it
			pdfVals[mm]=GaussianMixturePDF.pdf(ptCoverage,clusterMus[mm],clusterSigmas[mm]);
		}
		double max=pdfVals[0];
		int maxInd=0;
		for (int mm=0;mm<pdfVals.length;mm++){//search for the highest pdfVal
			//System.out.print(" //mm:"+mm +" pdfVal:"+pdfVals[mm]);
			if(pdfVals[mm]>max){
				max=pdfVals[mm];
				maxInd=mm;
			}
		}
		//System.out.println();
		return (double)rt.bestScore.bestCNVIndexes[maxInd];
	}

	public double getPointPloidyEstimationPoisson(double ptCoverage){//computes the probability of the data point 
		//belonging to all gaussians and returns the best option

		double [] pdfVals=new double[clusterMus.length];
		for (int mm=0;mm<clusterMus.length;mm++){//for all mixtures computes the probability pdfVal of the data point belonging to it
			pdfVals[mm]=PoissonMixturePDF.pdf(ptCoverage,clusterMus[mm]);
		}
		double max=pdfVals[0];
		int maxInd=0;
		for (int mm=0;mm<pdfVals.length;mm++){//search for the highest pdfVal
			//System.out.print(" //mm:"+mm +" pdfVal:"+pdfVals[mm]);
			if(pdfVals[mm]>max){
				max=pdfVals[mm];
				maxInd=mm;
			}
		}
		//System.out.println();
		return (double)rt.bestScore.bestCNVIndexes[maxInd];
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

	

	private int significantMaxsInPDF(NaivePDF naivePDF) {
		int sigMaxs = 0;//nb of significant maximums
		double threshold = SamParser.maxY * 0.001;// discards the values that are
													// below this threshold
		System.out.println(" SamParser.maxYh :" + SamParser.maxY + " Y min threshold:" + threshold);

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

					yMinList.add(maxVal);
					xMinList.add(naivePDF.xDataPoints[Xindex]);
					System.out.println(" ****    max in :" + naivePDF.xDataPoints[Xindex] + " = " + maxVal);
					sigMaxs++;

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
