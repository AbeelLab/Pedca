import java.io.File;
import java.io.IOException;

import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartUtilities;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.data.category.DefaultCategoryDataset;

public class BarChart {

	public BarChart (int [] readCounts) {
		// Create a simple Bar chart
		DefaultCategoryDataset dataset = new DefaultCategoryDataset();
		for (int r=0;r<readCounts.length;r++){
			
			dataset.setValue(readCounts[r], "#Contigs",(Integer)r);
		}
		
		JFreeChart chart = ChartFactory.createBarChart("Read Count Distribution", "Reads Count", "#Contigs", dataset,
				PlotOrientation.VERTICAL, false, true, false);
		try {
			ChartUtilities.saveChartAsJPEG(new File(Ploest.outputFile + "//" + Ploest.projectName+ "//readsDistribution.jpg"), chart, 5000, 3000);
		} catch (IOException e) {
			System.err.println("Problem occurred creating chart.");
		}
		//System.out.println("BarChar printed");
	}
}