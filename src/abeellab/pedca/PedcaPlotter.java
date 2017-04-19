package abeellab.pedca;
import java.awt.Color;
import java.awt.Rectangle;
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
import org.jfree.chart.ChartFrame;
import org.jfree.chart.ChartPanel;
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

import dataFitters.*;
import jMEF.MixtureModel;
import jMEF.PVector;

//******************************************************************
//
//
//DEPRECATED.
//(Was used to fit the read distribution with a Gaussian or a Poisson mixture)
//
//HAS BEEN REPLACED BY NaivePedcaPlotter
//
//
//
//******************************************************************

public class PedcaPlotter {
	static final int MAX_NB_MIXTURES=10;
	//static final double SIGMA_FACTOR=2;//accepted factor for standard variation criterion of a point belonging to a gaussian  
	Map<String,ContigData> contigsList;
	List<String> contArrList;
	static double[] gMMweights;//final result of the WEIGHTS of the gaussians in the mixture model fitting
	static double[] gMMmus;//final result of the MEANS (mus) of the gaussians in the mixture model fitting
	static double[] gMMsigmas;//final result of the STANDARD VARIATIONS (sigmas) of the gaussians in the mixture model fitting
	static GaussianDataFitter dfResult;//final data fit with the best fit predicted 
	static GaussianMixturePDF gmPDFResult;//final Gaussian Mixture PDF
	static PoissonDataFitter pdfResult;//final data fit with the best POISSON fit predicted 
	static PoissonMixturePDF pmPDFResult;//final POISSON Mixture PDF
	//PVector[] fitPoints;
	PVector[] intFitPoints;
	JFreeChart chart;
	static int maxX=0;
	//static int maxY=0;
	static int totalDataPoints=0;//total number of input datapoints (coverage for all windows)
	GaussianMixturePDF[] gmPDF;//Y points of the gaussian fit
	PoissonMixturePDF[] pmPDF;//Y points of the poisson fitpmPDF
	static int finalNumberOfMixtures;
	RatioFind rt;//contains the ratio of the gaussians to the ploidy unit which allows computation of ploidy from contig coverage

	public PedcaPlotter(Map<String,ContigData> contList,int maxWindows) {
		contigsList=contList;
		contArrList = new ArrayList<String>(contigsList.keySet());

		
		
		
		try{
			displayScatterPlot();
			createFitterDataset() ;
		//gaussianDensityEstimation();//TEST NEW METHOD
			
			
			fitPoissonMixtureModel();//fitGaussianMixtureModel();
			displayPloidyAndCoveragePlotPoisson();//displayPloidyAndCoveragePlot();
			
			//displayPloidyEstimationScatterPlot();

		}catch (Exception e){
			System.err.println("Error in PloestPlotter constructor");
		}
		
	}

	
	public void gaussianDensityEstimation(){
		System.out.println("-------------http://grepcode.com/file/repo1.maven.org/maven2/nz.ac.waikato.cms.weka/weka-stable/3.6.7/weka/estimators/KernelEstimator.java------------------");
		
	
		
	}
	
	
	public void fitPoissonMixtureModel(){
		System.out.println("-------------Fitting Different Poisso Mixture Models------------------");
		//System.out.println("-------------fitGaussianMixtureModel  MAX NB OF MIXTURE SET TO 90 !!!!------------------");
		double aic;
		double bic;
		double[]aicsEM=new double[MAX_NB_MIXTURES+1];//each EM (Expectation Maximization) AIC value is stored in its corresponding k(number of mixtures) index
		double[]bicsEM=new double[MAX_NB_MIXTURES+1];//same for EM BIC values
		double[]aicsBSC=new double[MAX_NB_MIXTURES+1];//each BSC (Bregman Soft Clustering)AIC value is stored in its corresponding k(number of mixtures) index
		double[]bicsBSC=new double[MAX_NB_MIXTURES+1];//same for BSC BIC values
		MixtureModel[]emMMs=new MixtureModel[MAX_NB_MIXTURES+1];//contains the result of the EM fit for each of the mixtures
		MixtureModel[]bscMMs=new MixtureModel[MAX_NB_MIXTURES+1];//contains the result of the BSC fit for each of the mixtures
		GaussianDataFitter gaussDF;
		PoissonDataFitter poissDF;
		int NbOfRuns=Pedca.nbOfRuns;
		int indofmin;
		int[] bestBSCGuess=new int[MAX_NB_MIXTURES];//best BSC guess with bic model evaluation
		int[] bestEMGuess=new int[MAX_NB_MIXTURES];//best BSC guess with EM model evaluation
		int[] correctedResults=new int[MAX_NB_MIXTURES];//best guess with bic model evaluation + min points evaluation

	
		pmPDF=new PoissonMixturePDF[NbOfRuns];
		if(Pedca.forceK==0){
			NbOfRuns=1;//DELETE WHEN DONE TESTING !!!!!!!
			for (int r=0;r<NbOfRuns;r++){
				System.out.println("run --"+r);
				for (int k=1;k<MAX_NB_MIXTURES;k++){//fit to different number k of mixtures
					System.out.println("run --"+r+" fitting "+k+" poissons/"+MAX_NB_MIXTURES);
					poissDF=new PoissonDataFitter (intFitPoints,k );

					emMMs[k]=poissDF.getEMmodel();//store EM fit values
					bscMMs[k]=poissDF.getBSCModel();//same for BSC
					
					//get EM LogLikelihoods and estimate BIC and AIC values

					aic=-2*(poissDF.getEMLogLikelihood())+(2*(2*k));//2*(2k) are the free parameters: k=number of models +2(weight, and lambda), the whle multiplied again by 2 as established by the AIC formula
					bic=-0.5*(poissDF.getEMLogLikelihood())+((2*k)*Math.log(SamParser.totalDataPoints));
				
					aicsEM[k]=aic;
					bicsEM[k]=bic;
					//get BSC LogLikelihoods and estimate BIC and AIC values
					aic=-2*(poissDF.getBSCLogLikelihood())+(2*(2*k));
					bic=-0.5*(poissDF.getBSCLogLikelihood())+((2*k)*Math.log(SamParser.totalDataPoints));
					aicsBSC[k]=aic;
					bicsBSC[k]=bic;	
					
					PoissonMixturePDF pmPDFtemp=new PoissonMixturePDF(poissDF.getEMmodel(),0.0,(double)SamParser.readCounts.length);
	//	SamParser.barchart.BarChartWithFit(pmPDFtemp,k,"tempFit");
				}

				System.out.println("  --");
				for (int k=0;k<aicsEM.length;k++){
					System.out.print("MIXTURE OF POISSONS "+k+"                \t");
				}System.out.println();
				for (int k=0;k<aicsEM.length;k++){
					System.out.print("             aicEM :"+ aicsEM[k]+"\t");
				}System.out.println();
				for (int k=0;k<aicsEM.length;k++){
					System.out.print("             bicEM:"+bicsEM[k]+"\t");
				}System.out.println();
				for (int k=0;k<aicsEM.length;k++){
					System.out.print("             AICBSC:"+aicsBSC[k]+"\t");
				}System.out.println();
				for (int k=0;k<aicsEM.length;k++){
					System.out.print("             BicBSC:"+bicsBSC[k]+"\t");
				}System.out.println();

				indofmin=findIndexOfMin(bicsEM);//find index of poisson mixture with min EM value
				System.out.println("   bestEMGuess --"+bestEMGuess[indofmin]);
				bestEMGuess[indofmin]++;
		
				System.out.println("     fitPOISSONMixtureModel indofmin --"+indofmin+ "mixt model params:"+emMMs[indofmin].printParams());
				
				pmPDF[r]=new PoissonMixturePDF(emMMs[indofmin],0.0,(double)SamParser.readCounts.length);
				int k=significantMinsInPDF( pmPDF[r]);//finds the significant minimums (below a threshold) in the gaussian mixture PDF
				//helps determine the real number of gaussians in the PDF and correct overfitting
				System.out.println("fitPoissonMixtureModel significantMinsInPDF ="+k);
				//SamParser.barchart.BarChartWithFit(gmPDF[r],(100+r));

				correctedResults[k]++;
			}
			//PRINT BEST GUESSES AND CORRECTED RESULTS		
			System.out.println("---------- EM prediction-------- ----------");//identical results with best guess EM
			System.out.print("[");
			for (int b=0;b<bestEMGuess.length;b++){
				System.out.print(" "+bestEMGuess[b]);
			}
			System.out.println(" ];");

			System.out.println("------- best correctedResults (after minimal points detection of PDF)----");
			System.out.print("[");
			for (int b=0;b<correctedResults.length;b++){
				System.out.print(" "+correctedResults[b]);
			}
			System.out.println(" ];");
			finalNumberOfMixtures=indexOfMode(correctedResults);
			System.out.println("-------------FINAL RESULT :"+finalNumberOfMixtures+" POISSON MODELS WITH "+(100*correctedResults[finalNumberOfMixtures]/NbOfRuns)+" % consensus---");
		}else{//force input number of mixtures
			finalNumberOfMixtures=Pedca.forceK;//
			//finalNumberOfMixtures=8;//RESET VALUE DELETING THIS
			//System.out.println("-------------FINAL RESULT, GAUSSIAN MODELS WITH NB OF MIXTURE SET TO 90 !!!! (instead of :"+finalNumberOfMixtures+" reset values in PloestPlotter) forced mixtures");
			System.out.println("-------------FINAL RESULT, POISSON MODELS WITH :"+finalNumberOfMixtures+"  forced mixtures");

		}

		//print final gaussian fitted result with final number of mixtures
		pdfResult=new PoissonDataFitter (intFitPoints,finalNumberOfMixtures );//final data fit with the best number of mixture prediction
		System.out.print("PoissonDataFitter for result done");
		gMMweights=new double[pdfResult.getBSCModel().weight.length];
		gMMmus=new double[pdfResult.getBSCModel().param.length];
	
		for (int g=0;g<pdfResult.getBSCModel().param.length;g++){

			gMMweights[g]=pdfResult.getBSCModel().weight[g];
			gMMmus[g]=((jMEF.PVector)pdfResult.getBSCModel().param[g]).array[0];
			
		}
		
		sortPMMs();//sort the result vectors (gMMmus,gMMweights and gMMsigmas) following ascending order of gMMmus
		
		pmPDFResult=new PoissonMixturePDF(pdfResult.getBSCModel(),0.0,(double)SamParser.readCounts.length);
		System.out.println("Ready to print FINALRESULT");

		rt=new RatioFind(gMMmus,100*correctedResults[finalNumberOfMixtures]/NbOfRuns);
		
		if(rt!=null)rt.writeOutPoisson();

	}
	
	public void fitGaussianMixtureModel(){
		System.out.println("-------------Fitting Different Gaussian Mixture Models------------------");
		//System.out.println("-------------fitGaussianMixtureModel  MAX NB OF MIXTURE SET TO 90 !!!!------------------");
		double aic;
		double bic;
		double[]aicsEM=new double[MAX_NB_MIXTURES+1];//each EM (Expectation Maximization) AIC value is stored in its corresponding k(number of mixtures) index
		double[]bicsEM=new double[MAX_NB_MIXTURES+1];//same for EM BIC values
		double[]aicsBSC=new double[MAX_NB_MIXTURES+1];//each BSC (Bregman Soft Clustering)AIC value is stored in its corresponding k(number of mixtures) index
		double[]bicsBSC=new double[MAX_NB_MIXTURES+1];//same for BSC BIC values
		MixtureModel[]emMMs=new MixtureModel[MAX_NB_MIXTURES+1];//contains the result of the EM fit for each of the mixtures
		MixtureModel[]bscMMs=new MixtureModel[MAX_NB_MIXTURES+1];//contains the result of the BSC fit for each of the mixtures
		GaussianDataFitter df;
		int NbOfRuns=Pedca.nbOfRuns;
		int indofmin;
		int[] bestBSCGuess=new int[MAX_NB_MIXTURES];//best BSC guess with bic model evaluation
		int[] bestEMGuess=new int[MAX_NB_MIXTURES];//best BSC guess with EM model evaluation
		int[] correctedResults=new int[MAX_NB_MIXTURES];//best guess with bic model evaluation + min points evaluation

		gmPDF =new GaussianMixturePDF[NbOfRuns];
		pmPDF=new PoissonMixturePDF[NbOfRuns];
		if(Pedca.forceK==0){
			for (int r=0;r<NbOfRuns;r++){
				System.out.println("run --"+r);
				for (int k=0;k<MAX_NB_MIXTURES;k++){//fit to different number k of mixtures
					df=new GaussianDataFitter (intFitPoints,k );//fits a gauss mixture(by EM and BSC) to the 
					emMMs[k]=df.getEMmodel();//store EM fit values
					bscMMs[k]=df.getBSCModel();//same for BSC
					//get EM LogLikelihoods and estimate BIC and AIC values
					aic=-2*(df.getEMLogLikelihood())+(2*(k*3));//k+3 are the free parameters: k=number of models +3(weight,mu and sigma)
					bic=-0.5*(df.getEMLogLikelihood())+((k*3)*Math.log(SamParser.totalDataPoints));
					aicsEM[k]=aic;
					bicsEM[k]=bic;
					//get BSC LogLikelihoods and estimate BIC and AIC values
					aic=-2*(df.getBSCLogLikelihood())+(2*(k*3));
					bic=-0.5*(df.getBSCLogLikelihood())+((k*3)*Math.log(SamParser.totalDataPoints));
					aicsBSC[k]=aic;
					bicsBSC[k]=bic;	
				}


				indofmin=findIndexOfMin(bicsEM);//find index of gaussian mixture with min EM value
				bestEMGuess[indofmin]++;
				indofmin=findIndexOfMin(bicsBSC);//find index of gaussian mixture with min BIC value			
				bestBSCGuess[indofmin]++;//System.out.println(indofmin+" MM  BSC params:"+bscMMs[indofmin].printParams());
				//bscMMs[indofmin].printParams()
				//System.out.println("fitGaussianMixtureModel indofmin --"+indofmin+ "mixt model params:"+bscMMs[indofmin].printParams());
				//bscMMs[indofmin]
				gmPDF[r]=new GaussianMixturePDF(bscMMs[indofmin],0.0,(double)SamParser.readCounts.length);
				int k=significantMinsInPDF( gmPDF[r]);//finds the significant minimums (below a threshold) in the gaussian mixture PDF
				//helps determine the real number of gaussians in the PDF and correct overfitting
				System.out.println("fitGaussianMixtureModel significantMinsInPDF ="+k);
				//SamParser.barchart.BarChartWithFit(gmPDF[r],(100+r));

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
			finalNumberOfMixtures=indexOfMode(correctedResults);
			System.out.println("-------------FINAL RESULT :"+finalNumberOfMixtures+" GAUSSIAN MODELS WITH "+(100*correctedResults[finalNumberOfMixtures]/NbOfRuns)+" % consensus---");
		}else{//force input number of mixtures
			finalNumberOfMixtures=Pedca.forceK;//
			//finalNumberOfMixtures=8;//RESET VALUE DELETING THIS
			//System.out.println("-------------FINAL RESULT, GAUSSIAN MODELS WITH NB OF MIXTURE SET TO 90 !!!! (instead of :"+finalNumberOfMixtures+" reset values in PloestPlotter) forced mixtures");
			System.out.println("-------------FINAL RESULT, GAUSSIAN MODELS WITH :"+finalNumberOfMixtures+"  forced mixtures");

		}

		//print final gaussian fitted result with final number of mixtures
		dfResult=new GaussianDataFitter (intFitPoints,finalNumberOfMixtures );//final data fit with the best number of mixture prediction

		gMMweights=new double[dfResult.getBSCModel().weight.length];
		gMMmus=new double[dfResult.getBSCModel().param.length];
		gMMsigmas=new double[dfResult.getBSCModel().param.length];
		for (int g=0;g<dfResult.getBSCModel().param.length;g++){
			gMMweights[g]=dfResult.getBSCModel().weight[g];
			gMMmus[g]=((jMEF.PVector)dfResult.getBSCModel().param[g]).array[0];
			gMMsigmas[g]=((jMEF.PVector)dfResult.getBSCModel().param[g]).array[1];			
		}
		sortGMMs();//sort the result vectors (gMMmus,gMMweights and gMMsigmas) following ascending order of gMMmus
		gmPDFResult=new GaussianMixturePDF(dfResult.getBSCModel(),0.0,(double)SamParser.readCounts.length);
		SamParser.barchart.BarChartWithFit(gmPDFResult,finalNumberOfMixtures);

		rt=new RatioFind(gMMmus,100*correctedResults[finalNumberOfMixtures]/NbOfRuns);
		if(rt!=null)rt.writeOut();
	}

	private void sortGMMs() {//sort gMMmu in ascending order, then sort their weights and sigmas in the
		//corresponding order

		Vector<Double> origCopy = new Vector<Double>(gMMmus.length);//get a copy of the original order
		for (int o=0;o<gMMmus.length;o++){
			origCopy.add(gMMmus[o]);
		}	

		java.util.Arrays.sort(gMMmus);//sort the reference vector
		int[]order=new int[origCopy.size()];
		for (int o=0;o<order.length;o++){////store the order
			order[o]=origCopy.indexOf(gMMmus[o]);
		}
		//now order the otherGMM
		double[]newgMMsigmas=new double[dfResult.getBSCModel().param.length];//newly order sigmas
		double[]newgMMweights=new double[dfResult.getBSCModel().param.length];//newly order weights
		for (int o=0;o<order.length;o++){
			newgMMsigmas[o]=gMMsigmas[order[o]];
			newgMMweights[o]=gMMweights[order[o]];
		}
		gMMsigmas=newgMMsigmas;
		gMMweights=newgMMweights;
	}
	
	
	private void sortPMMs() {//sort gMMmu in ascending order, then sort their weights and sigmas in the
		//corresponding order

		Vector<Double> origCopy = new Vector<Double>(gMMmus.length);//get a copy of the original order
		for (int o=0;o<gMMmus.length;o++){
			origCopy.add(gMMmus[o]);
		}	

		java.util.Arrays.sort(gMMmus);//sort the reference vector
		int[]order=new int[origCopy.size()];
		for (int o=0;o<order.length;o++){////store the order
			order[o]=origCopy.indexOf(gMMmus[o]);
		}
	
		//now order the otherGMM
		double[]newgMMweights=new double[pdfResult.getBSCModel().param.length];//newly order weights
		for (int o=0;o<order.length;o++){
			newgMMweights[o]=gMMweights[order[o]];
		}
		gMMweights=newgMMweights;

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
			rangeAxis.setRange(0.00, SamParser.readsDistributionMaxCoverage);
			ChartUtilities.saveChartAsJPEG(new File(Pedca.outputFile + "//" + Pedca.projectName+ "//Contig_Coverage_Charts//Chart_Contig_"+c+".jpg"), chart, 1000, 600);
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
			ChartUtilities.saveChartAsJPEG(new File(Pedca.outputFile + "//" + Pedca.projectName+ "//Ploidy_Estimation_Charts//Ploidy_Estimation_Contig_"+c+".jpg"),chart, 1000, 600);

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
			ChartUtilities.saveChartAsJPEG(new File(Pedca.outputFile + "//" + Pedca.projectName+ "//Ploidy_Estimation_Charts//Ploidy_Estimation_Contig_"+c+".jpg"),chart, 1000, 600);

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

		ChartUtilities.saveChartAsJPEG(new File(Pedca.outputFile + "//" + Pedca.projectName+ "//Ploidy_Estimation_Charts//Ploidy_Estimation_Contig_"+contigD.contigName+".jpg"),chart, 1000, 600);

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

			ChartUtilities.saveChartAsJPEG(new File(Pedca.outputFile + "//" + Pedca.projectName+ "//Ploidy_Estimation_Charts//Ploidy_Estimation_Contig_"+contigD.contigName+".jpg"),chart, 1000, 600);

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

		intFitPoints=SamParser.fitPoints;

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

		double [] pdfVals=new double[gMMmus.length];
		for (int mm=0;mm<gMMmus.length;mm++){//for all mixtures computes the probability pdfVal of the data point belonging to it
			pdfVals[mm]=GaussianMixturePDF.pdf(ptCoverage,gMMmus[mm],gMMsigmas[mm]);
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

		double [] pdfVals=new double[gMMmus.length];
		for (int mm=0;mm<gMMmus.length;mm++){//for all mixtures computes the probability pdfVal of the data point belonging to it
			pdfVals[mm]=PoissonMixturePDF.pdf(ptCoverage,gMMmus[mm]);
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

	private int significantMinsInPDF(GaussianMixturePDF gmf) {
		//MaxMinArrays lists = new MaxMinArrays();
		int sigMins=0;
		double threshold=Pedca.SIGNIFICANT_MIN;//discards the mins that are above this threshold 
		//First find Mins:
		//int count = 0; // To handle special case of singleton list
		ArrayList<Double> yMinList= new  ArrayList<Double>();
		ArrayList<Double> xMinList= new  ArrayList<Double>();
		double left  = gmf.maxYvalue;
		double mid   = gmf.maxYvalue;
		double right = gmf.maxYvalue;
		int ind=0;
		while (ind<gmf.yDataPoints.length) {
			//count++;
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

	private int significantMinsInPDF(PoissonMixturePDF pmf) {
		//MaxMinArrays lists = new MaxMinArrays();
		int sigMins=0;
		double threshold=Pedca.SIGNIFICANT_MIN;//discards the mins that are above this threshold 
		//First find Mins:
		//int count = 0; // To handle special case of singleton list
		ArrayList<Double> yMinList= new  ArrayList<Double>();
		ArrayList<Double> xMinList= new  ArrayList<Double>();
		double left  = pmf.maxYvalue;
		double mid   = pmf.maxYvalue;
		double right = pmf.maxYvalue;
		int ind=0;
		while (ind<pmf.yDataPoints.length) {
			//count++;
			if(pmf.yDataPoints[ind]!=right){
				left = mid;
				mid = right;
			}	
			right = pmf.yDataPoints[ind];

			if (ind!=1 && right > mid && mid < left && mid < threshold) {// we count mins only after index 2 and above threshold
				yMinList.add(mid);
				xMinList.add(pmf.xDataPoints[ind]);					
				sigMins++;
			}

			ind++;
		}
		//add last point as min
		yMinList.add(mid);
		xMinList.add(pmf.xDataPoints[ind-1]);					
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
