import java.io.File;
import java.io.IOException;

import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartUtilities;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.data.category.DefaultCategoryDataset;

public class BarChart {
	JFreeChart chart;
	DefaultCategoryDataset dataset;

	public BarChart (int [] readCounts) {
		float [] normReadCounts=normalize(readCounts);
		// Create a simple Bar chart
		
		dataset = new DefaultCategoryDataset();
		for (int r=0;r<normReadCounts.length;r++){			
			dataset.setValue(normReadCounts[r], "#Contigs",(Integer)r);
		}
		
		chart = ChartFactory.createBarChart("Read Count Distribution", "Reads Count", "%Contigs", dataset,
				PlotOrientation.VERTICAL, false, true, false);
		try {
			ChartUtilities.saveChartAsJPEG(new File(Ploest.outputFile + "//" + Ploest.projectName+ "//readsDistribution.jpg"), chart, 1000, 600);
		} catch (IOException e) {
			System.err.println("Problem occurred creating chart.");
		}
		//System.out.println("BarChar printed");
	}

	private float[] normalize(int[] readCounts) {
		float[] normalizedReadCounts=new float[readCounts.length];
		int sum=0;
		for (int i=0;i<readCounts.length;i++){
			sum+=readCounts[i];
		}
		System.out.println("SUM:"+sum);
		for (int i=0;i<readCounts.length;i++){
			normalizedReadCounts[i]=(float)readCounts[i]/sum;
		}
		return normalizedReadCounts;
	}
}