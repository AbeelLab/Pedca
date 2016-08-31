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

import DataFitter.DattaFiter;
import jMEF.PVector;

public class PloestPlotter {
	Map<String,ContigData> contigsList;
	PVector[] fitPoints;
	static int maxX=0;
	static int maxY=0;
	static int totalDataPoints=0;//total number of input datapoints (coverage for all windows)
	
	public PloestPlotter(Map<String,ContigData> contList,int maxWindows) {
		System.out.println("maxWindows:"+maxWindows);
		contigsList=contList;
		maxY=350;
		try{
		displayScatterPlot();
		createFitterDataset() ;
		DattaFiter df;
		for (int k=1;k<10;k++){//fit to different number k of mixtures
			df=new DattaFiter (fitPoints,k );
			System.out.println(" k:"+k+" loglikelihood:"+df.getGaussLogLikelihood());
		}
	
		}catch (Exception e){
        	
        }
	}
	
	
	public void displayScatterPlot() throws IOException{
		List<String> contArrList = new ArrayList<String>(contigsList.keySet());
		for (int c=0;c<contigsList.size();c++){//for each contig
			ContigData contigD=contigsList.get(contArrList.get(c));
			JFreeChart chart = ChartFactory.createScatterPlot(
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
   // System.out.print("maxX:"+maxX);

    totalDataPoints+=maxX+1;
    //System.out.println("cds totalDataPoints:"+totalDataPoints);
    return result;
	}


	private  void createFitterDataset() {
		fitPoints=new PVector[totalDataPoints];
		int ind=0;
		List<String> contArrList = new ArrayList<String>(contigsList.keySet());
	
		for (int c=0;c<contigsList.size();c++){//for each contig
			ContigData contigD=contigsList.get(contArrList.get(c));
			for (int i = 0; i < (contigD.windPos.length); i++) {//for each window position
				PVector curVec=new PVector(1);
	            curVec.array[0]=contigD.windPos[i];
	            fitPoints[ind++]=curVec;               
		    }
		}

	}



}
