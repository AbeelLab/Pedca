package abeellab.pedca;
import java.awt.Font;
import java.io.File;
import java.io.IOException;

import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartUtilities;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.axis.AxisSpace;
import org.jfree.chart.axis.CategoryAxis;
import org.jfree.chart.axis.CategoryLabelPositions;
import org.jfree.chart.axis.NumberAxis;
import org.jfree.chart.axis.NumberTickUnit;
import org.jfree.chart.axis.ValueAxis;
import org.jfree.chart.labels.StandardXYToolTipGenerator;
import org.jfree.chart.plot.DatasetRenderingOrder;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.renderer.xy.StandardXYItemRenderer;
import org.jfree.chart.renderer.xy.XYAreaRenderer;
import org.jfree.chart.renderer.xy.XYItemRenderer;
import org.jfree.chart.renderer.xy.XYLineAndShapeRenderer;
import org.jfree.data.category.DefaultCategoryDataset;
import org.jfree.data.time.TimeSeries;
import org.jfree.data.time.TimeSeriesCollection;
import org.jfree.data.xy.XYDataset;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;

import jMEF.MixtureModel;
import jMEF.PVector;

public class BarChart {
	JFreeChart histChart;
	JFreeChart overlaidChart;
	DefaultCategoryDataset histDataset;
	float [] normReadCounts;
	int maxX;
	float maxY=0;
	int NB_OF_BASECALL_BiNS=20;
	Font fontLabels = new Font("Dialog", Font.PLAIN, 25); 
	Font fontTitle = new Font("Dialog", Font.PLAIN, 30); 
	Font fontLegend = new Font("Dialog", Font.PLAIN, 20); 
	int sumOfReadCounts=0;
	//Constructor for Read Count distribution
	public BarChart (int [] readCounts) {
		normReadCounts=normalize(readCounts);
		// Create a simple Bar chart
		histDataset = new DefaultCategoryDataset();
		for (int r=0;r<normReadCounts.length;r++){			
			if(normReadCounts[r]>maxY)maxY=normReadCounts[r];
			if(r >maxX)maxX=r;
		}
		int domainTick=maxX/20;

		for (int r=0;r<normReadCounts.length;r++){			
			histDataset.setValue(normReadCounts[r], "#Contigs",new DomainCategory(r,domainTick));
		}

		// Create a simple Bar chart
		histChart = ChartFactory.createBarChart("Read Count Distribution. wl: "+Pedca.windowLength+" "+SamParser.stringSecondRound, "Reads Count", "%Contigs", histDataset,
				PlotOrientation.VERTICAL, false, true, false);

		histChart.getCategoryPlot().getRangeAxis().setAutoRange(true);
		histChart.getCategoryPlot().getDomainAxis().setCategoryLabelPositions(CategoryLabelPositions.UP_45);
		histChart.getCategoryPlot().getDomainAxis().setLabelFont(fontLabels);
		histChart.getCategoryPlot().getRangeAxis().setLabelFont(fontLabels);
		histChart.getCategoryPlot().getDomainAxis().setTickLabelFont(fontLabels);
		histChart.getCategoryPlot().getRangeAxis().setTickLabelFont(fontLabels);
		histChart.getTitle().setFont(fontTitle);
	
		
		try {
			ChartUtilities.saveChartAsJPEG(new File(Pedca.outputFile + "//" + Pedca.projectName+ "//readsDistribution"+SamParser.stringSecondRound+".jpg"), histChart, 1500, 900);
		} catch (IOException e) {
			System.err.println("Problem occurred creating chart.");
		}
		System.out.println("BarChar printed "+SamParser.stringSecondRound+". maxY="+maxY);
	}
	
	
	//baseCall chart constructor
	public BarChart (double [] baseCalls,int cluster, int maxY) {
		 
		fontLabels = new Font("Dialog", Font.PLAIN, 30); 
		fontTitle = new Font("Dialog", Font.PLAIN, 35); 
		fontLegend = new Font("Dialog", Font.PLAIN, 25); 
		 
		int[] bins=new int [NB_OF_BASECALL_BiNS];
		int bin;
		Double product;
		
		for (int i=0;i<baseCalls.length;i++){//fill the bins
			
			product=baseCalls[i]*NB_OF_BASECALL_BiNS;
			product=(product-(product%1));//(gets the right bin where this value goes)
			bin=product.intValue();
			bins[bin]++;
		}

		// Create a simple Bar chart
		Integer perc;
		histDataset = new DefaultCategoryDataset();
		for (int r=0;r<bins.length;r++){	
			perc=(r*5);
			histDataset.setValue(bins[r], "Base Call %",perc);
		}
		
		histChart = ChartFactory.createBarChart("BaseCall Distribution. Cluster nb:"+cluster, "Base Call %", "Number of occurrences ", histDataset,
				PlotOrientation.VERTICAL, false, true, false);
		
		histChart.getCategoryPlot().getRangeAxis().setRange(0.00, maxY);
        histChart.getCategoryPlot().getDomainAxis().setLabelFont(fontLabels);
        histChart.getCategoryPlot().getRangeAxis().setLabelFont(fontLabels);
        histChart.getCategoryPlot().getDomainAxis().setTickLabelFont(fontLabels);
        histChart.getCategoryPlot().getRangeAxis().setTickLabelFont(fontLabels);
  
        histChart.getTitle().setFont(fontTitle);
        
		
		try {
			ChartUtilities.saveChartAsJPEG(new File(Pedca.outputFile + "//BaseCall//BaseCallHistogramCluster_"+cluster+".jpg"), histChart, 1500, 900);
		} catch (IOException e) {
			System.err.println("Problem occurred creating base call chart.");
		}
		System.out.println("BarChar BaseCall cluster "+cluster+" printed ");
	}

	//@SuppressWarnings("deprecation")
	public void BarChartWithFit (GaussianMixturePDF gaussFit,int r) {
		
		if (gaussFit==null){
			System.out.println("BarChartWithFit gaussFit==null r:"+r);
		}else{
			
			final XYDataset data1 = createHistDataset(gaussFit) ;//histogram of readCounts
			final XYItemRenderer renderer1 = new StandardXYItemRenderer();

			final NumberAxis domainAxis = new NumberAxis("ReadsCounts");
			//domainAxis.setTickMarkPosition(DateTickMarkPosition.MIDDLE);
			final ValueAxis rangeAxis = new NumberAxis("%Contigs");
			final XYPlot plot = new XYPlot  (data1, domainAxis, rangeAxis, renderer1);


			// add a second dataset and renderer...
			final XYDataset data2 = createFitCurveDataset(gaussFit);
			final XYItemRenderer renderer2 = new StandardXYItemRenderer();

			plot.setDataset(1, data2);
			plot.setRenderer(1, renderer2);

			plot.setDatasetRenderingOrder(DatasetRenderingOrder.REVERSE);

			// return a new chart containing the overlaid plot...
			overlaidChart=new JFreeChart("Gauss Mixture Model Fit of Reads Distribution", JFreeChart.DEFAULT_TITLE_FONT, plot, true);


			try {
				ChartUtilities.saveChartAsJPEG(new File(Pedca.outputFile + "//" + Pedca.projectName+ "//readsDistributionGaussianFitted"+r+".jpg"), overlaidChart, 1000, 600);
			} catch (IOException e) {
				System.err.println("Problem occurred creating chart.");
			}
		}
	}

	
	

public void BarChartWithFit (NaivePDF naivePDF, String title) {

	 fontLabels = new Font("Dialog", Font.PLAIN, 20); 
	 fontTitle = new Font("Dialog", Font.PLAIN, 25); 
	 fontLegend = new Font("Dialog", Font.PLAIN, 18); 
	
	
	if (naivePDF==null){
		System.out.println("BarChartWithFit naive Fit==null");
	}else{
		
		final XYDataset data1 = createHistDataset(naivePDF) ;//histogram of readCounts
		final XYItemRenderer renderer1 = new StandardXYItemRenderer();
		final NumberAxis domainAxis = new NumberAxis("ReadsCounts");
		final ValueAxis rangeAxis = new NumberAxis("%Contigs");
		final XYPlot plot = new XYPlot  (data1, domainAxis, rangeAxis, renderer1);
		
		// add a second dataset and renderer...
		final XYDataset data2 = createFitCurveDataset(naivePDF);
		final XYItemRenderer renderer2 = new StandardXYItemRenderer();

		plot.setDataset(1, data2);
		plot.setRenderer(1, renderer2);
	
		plot.getDomainAxis().setLabelFont(fontLabels);
		plot.getRangeAxis().setLabelFont(fontLabels);
		plot.getDomainAxis().setTickLabelFont(fontLabels);
		plot.getRangeAxis().setTickLabelFont(fontLabels);
		
				
		plot.setDatasetRenderingOrder(DatasetRenderingOrder.FORWARD);
		// return a new chart containing the overlaid plot...
		overlaidChart=new JFreeChart("Read Distribution Fit. wl:"+Pedca.windowLength+" bp", JFreeChart.DEFAULT_TITLE_FONT, plot, true);
		overlaidChart.getTitle().setFont(fontTitle);
		overlaidChart.getLegend().setItemFont(fontLegend);
		try {
			ChartUtilities.saveChartAsJPEG(new File(Pedca.outputFile + "//" + Pedca.projectName+ "//readsDistributionFitted"+title+".jpg"), overlaidChart, 1000, 600);
		} catch (IOException e) {
			System.err.println("Problem occurred creating chart.");
		}
		System.out.println("BarChar of naive PDF printed");
	}
}
	

private XYDataset createHistDataset( GaussianMixturePDF gaussFit) {
	XYSeriesCollection result = new XYSeriesCollection();
	XYSeries series = new XYSeries("Reads Counts");
	//maxX=0;

	double ind=(gaussFit.beg-0.5);

	for (int i = 0; i<normReadCounts.length; i++) {
		while (ind<i+0.5){
			//System.out.println("|ind:"+ind+" i:"+i+" x:"+normReadCounts[i]+" y:"+normReadCounts[i]);
			series.add(ind, normReadCounts[i]);
			ind+=gaussFit.step;
		}

	}
	result.addSeries(series);

	return result;
}

	private XYDataset createHistDataset( NaivePDF naiveFit) {
		XYSeriesCollection result = new XYSeriesCollection();
		XYSeries series = new XYSeries("Reads Counts");
		

		double ind=(naiveFit.beg-0.5);

		for (int i = 0; i<normReadCounts.length; i++) {
			while (ind<i+0.5){
				//System.out.println("|ind:"+ind+" i:"+i+" x:"+normReadCounts[i]+" y:"+normReadCounts[i]+" naiveFit.step:"+naiveFit.step);
				series.add(ind, normReadCounts[i]);
				ind+=naiveFit.step;
				if(normReadCounts[i]>naiveFit.maxYHISTOGRAMvalue)naiveFit.maxYHISTOGRAMvalue=normReadCounts[i];
			}

		}
		result.addSeries(series);
		return result;
	}
	private XYDataset createHistDataset( PoissonMixturePDF gaussFit) {
		XYSeriesCollection result = new XYSeriesCollection();
		XYSeries series = new XYSeries("Reads Counts");
		//maxX=0;

		double ind=(gaussFit.beg-0.5);

		for (int i = 0; i<normReadCounts.length; i++) {
			while (ind<i+0.5){
				//System.out.println("|ind:"+ind+" i:"+i+" x:"+normReadCounts[i]+" y:"+normReadCounts[i]);
				series.add(ind, normReadCounts[i]);
				ind+=gaussFit.step;
			}

		}
		result.addSeries(series);

		return result;
	}

	private XYDataset createFitCurveDataset(NaivePDF naiveFit) {
		XYSeriesCollection result = new XYSeriesCollection();
		XYSeries series = new XYSeries("PDF Naive Smoother Fit");
		
		//correction for overlapping the 2 datasets
		int xOffSet=naiveFit.smootherWing/2;//to correct for the bins width
		naiveFit.yRatioCorrection=naiveFit.maxYHISTOGRAMvalue/naiveFit.maxYFITvalue;//to correct the y data normalization ratio
		//System.out.println("createFitCurveDataset naiveFit.maxy="+naiveFit.maxYFITvalue+" ");
		//add to series
		for (int i = 0; i<naiveFit.yDataPoints.length; i++) {
			series.add(naiveFit.xDataPoints[i]+xOffSet, naiveFit.yDataPoints[i]*naiveFit.yRatioCorrection);
		}
		
		result.addSeries(series);//and send
		
		return result;
	}
	
	private XYDataset createFitCurveDataset(GaussianMixturePDF gaussFit) {
		XYSeriesCollection result = new XYSeriesCollection();
		XYSeries series = new XYSeries("Gaussian Mixture Fit");
		
		for (int i = 0; i<gaussFit.yDataPoints.length; i++) {

			series.add(gaussFit.xDataPoints[i], gaussFit.yDataPoints[i]);
			//if (x>maxX)maxX=(int) x;           
		}
		
		result.addSeries(series);//and send
		
		return result;
	}
	
	private XYDataset createFitCurveDataset(PoissonMixturePDF gaussFit) {
		XYSeriesCollection result = new XYSeriesCollection();
		XYSeries series = new XYSeries("Poisson Mixture Fit");
		System.out.print("createFitCurveDataset [ ");
		for (int i = 0; i<gaussFit.yDataPoints.length; i++) {
			series.add(gaussFit.xDataPoints[i], gaussFit.yDataPoints[i]);
			System.out.print(" "+gaussFit.xDataPoints[i]+" "+gaussFit.yDataPoints[i]+";");          
		}System.out.println(" ]");
		
		result.addSeries(series);//and send
		
		return result;
	}

	private float[] normalize(int[] readCounts) {
		float[] normalizedReadCounts=new float[readCounts.length];
		
		for (int i=0;i<readCounts.length;i++){
			sumOfReadCounts+=readCounts[i];
		}
		//System.out.println("SUM:"+sum);
		for (int i=0;i<readCounts.length;i++){
			normalizedReadCounts[i]=100*(float)readCounts[i]/sumOfReadCounts;
		}
		return normalizedReadCounts;
	}


	class DomainCategory implements Comparable<DomainCategory> {
		   Integer value;
		   String stValue;

		   DomainCategory(int val,int tick) {
		      value = val; 
		      stValue = val%tick==0?""+val:"";}

		   public int compareTo(DomainCategory key) { return value.compareTo(key.value); }

		   public String toString() { return stValue; }
		}
}