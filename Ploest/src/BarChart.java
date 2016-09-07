import java.io.File;
import java.io.IOException;

import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartUtilities;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.axis.NumberAxis;
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

public class BarChart {
	JFreeChart histChart;
	JFreeChart overlaidChart;
	DefaultCategoryDataset histDataset;
	float [] normReadCounts;
	int maxX;

	public BarChart (int [] readCounts) {
		normReadCounts=normalize(readCounts);
		// Create a simple Bar chart

		histDataset = new DefaultCategoryDataset();
		for (int r=0;r<normReadCounts.length;r++){			
			histDataset.setValue(normReadCounts[r], "#Contigs",(Integer)r);
		}

		histChart = ChartFactory.createBarChart("Read Count Distribution", "Reads Count", "%Contigs", histDataset,
				PlotOrientation.VERTICAL, false, true, false);
		try {
			ChartUtilities.saveChartAsJPEG(new File(Ploest.outputFile + "//" + Ploest.projectName+ "//readsDistribution.jpg"), histChart, 1000, 600);
		} catch (IOException e) {
			System.err.println("Problem occurred creating chart.");
		}
		//System.out.println("BarChar printed");
	}

	//@SuppressWarnings("deprecation")
	public void BarChartWithFit (GaussianMixturePDF gaussFit,int r) {
		if (gaussFit==null)System.out.println("BarChartWithFit gaussFit==null");
		final XYDataset data1 = createHistDataset(gaussFit) ;//histogram of readCounts
		final XYItemRenderer renderer1 = new StandardXYItemRenderer();

		final NumberAxis domainAxis = new NumberAxis("ReadsCounts");
		//domainAxis.setTickMarkPosition(DateTickMarkPosition.MIDDLE);
		final ValueAxis rangeAxis = new NumberAxis("%Contigs");
		final XYPlot plot = new XYPlot  (data1, domainAxis, rangeAxis, renderer1);


		// add a second dataset and renderer...
		final XYDataset data2 = createCurveDataset(gaussFit);
		final XYItemRenderer renderer2 = new StandardXYItemRenderer();

		plot.setDataset(1, data2);
		plot.setRenderer(1, renderer2);

		plot.setDatasetRenderingOrder(DatasetRenderingOrder.FORWARD);

		// return a new chart containing the overlaid plot...
		overlaidChart=new JFreeChart("Overlaid Plot Example", JFreeChart.DEFAULT_TITLE_FONT, plot, true);


		try {
			ChartUtilities.saveChartAsJPEG(new File(Ploest.outputFile + "//" + Ploest.projectName+ "//readsDistributionGaussianFitted"+r+".jpg"), overlaidChart, 1000, 600);
		} catch (IOException e) {
			System.err.println("Problem occurred creating chart.");
		}
		//System.out.println("BarChar printed");
	}

	private XYDataset createHistDataset( GaussianMixturePDF gaussFit) {
		XYSeriesCollection result = new XYSeriesCollection();
		XYSeries series = new XYSeries("Reads Counts");
		//maxX=0;
		
		if (gaussFit==null)System.out.println("gaussFit==null");
		
		double ind=(gaussFit.beg-0.5);
		
		//System.out.println(" createHistDataset ; xDataPoints.length:"+gaussFit.xDataPoints.length+" normReadCounts.length:"+normReadCounts.length);
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

	private XYDataset createCurveDataset(GaussianMixturePDF gaussFit) {
		XYSeriesCollection result = new XYSeriesCollection();
		XYSeries series = new XYSeries("Gaussian Mixture Fit");
		//maxX=0;


		for (int i = 0; i<gaussFit.yDataPoints.length; i++) {
	
			series.add(gaussFit.xDataPoints[i], gaussFit.yDataPoints[i]);
			//if (x>maxX)maxX=(int) x;           
		}
		result.addSeries(series);

		return result;
	}

	private float[] normalize(int[] readCounts) {
		float[] normalizedReadCounts=new float[readCounts.length];
		int sum=0;
		for (int i=0;i<readCounts.length;i++){
			sum+=readCounts[i];
		}
		//System.out.println("SUM:"+sum);
		for (int i=0;i<readCounts.length;i++){
			normalizedReadCounts[i]=(float)readCounts[i]/sum;
		}
		return normalizedReadCounts;
	}


}