package dataFitters;

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

public class GaussianDattaFiter {

	double emLogLikelihood;//EM logLikelihood of this gauss mix model
	double bscLogLikelihood;//BSC logLikelihood of this gauss mix model
	MixtureModel mmc;// Expectation Maximization Mixture Model
	MixtureModel mmef;//Bregman soft clustering Mixture Model
	/**
	 * Main function.
	 * @param args
	 */
	public GaussianDattaFiter (PVector[] points,int n ) {//fit the datapoints to a mixture of n gaussians 

		// Display
		String title = "";
		title += "+---------------------------------------------+\n";
		title += "| EM Gauss Fitter | (with K-means verification)\n";
		title += "+---------------------------------------------+\n";
		//System.out.print(title);

		// Variables

		Vector<PVector>[] clusters = KMeans.run(points, n);

	
		// Classical EM
		
		mmc = ExpectationMaximization1D.initialize(clusters);
		mmc = ExpectationMaximization1D.run(points, mmc);
		emLogLikelihood=mmc.getEMLogLikelihod();
		//System.out.println("Mixure model estimated using classical EM \n" + mmc + "\n");
		
		
		// Bregman soft clustering
		//MixtureModel mmef;
		mmef = BregmanSoftClustering.initialize(clusters, new UnivariateGaussian());
		mmef = BregmanSoftClustering.run(points, mmef);
		bscLogLikelihood=mmef.getBSCLogLikelihod();
		//System.out.println("Mixure model estimated using Bregman soft clustering \n" + mmef + "\n");

	}

	public double getEMLogLikelihood(){
		return emLogLikelihood;
	}
	
	public double getBSCLogLikelihood(){
		return bscLogLikelihood;
	}
	
	public MixtureModel getEMmodel(){
		return mmc;
	}
	
	public MixtureModel getBSCModel(){
		return mmef;
	}
}
