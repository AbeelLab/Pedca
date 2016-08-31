package DataFitter;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.Scanner;
import java.util.Vector;

import Tools.ExpectationMaximization1D;
import Tools.KMeans;
import jMEF.BregmanSoftClustering;
import jMEF.MixtureModel;
import jMEF.PVector;
import jMEF.UnivariateGaussian;

public class DattaFiter {

	double gaussLogLikelihood;//logLikelihood of this gauss mix model
	/**
	 * Main function.
	 * @param args
	 */
	public DattaFiter (PVector[] points,int n ) {//fit the datapoints to a mixture of n gaussians 

		// Display
		String title = "";
		title += "+---------------------------------------------+\n";
		title += "| EM Gauss Fitter | (with K-means verification)\n";
		title += "+---------------------------------------------+\n";
		System.out.print(title);

		// Variables

		Vector<PVector>[] clusters = KMeans.run(points, n);

	
		// Classical EM
		MixtureModel mmc;
		mmc = ExpectationMaximization1D.initialize(clusters);
		mmc = ExpectationMaximization1D.run(points, mmc);
		gaussLogLikelihood=mmc.getLogLikelihod();
		//System.out.println("Mixure model estimated using classical EM \n" + mmc + "\n");
		
		
		// Bregman soft clustering
		MixtureModel mmef;
		mmef = BregmanSoftClustering.initialize(clusters, new UnivariateGaussian());
		mmef = BregmanSoftClustering.run(points, mmef);
		//System.out.println("Mixure model estimated using Bregman soft clustering \n" + mmef + "\n");

	}

	public double getGaussLogLikelihood(){
		return gaussLogLikelihood;
	}
	
}
