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
import jMEF.Poisson;


public class PoissonDataFitter {

	double emLogLikelihood;//EM logLikelihood of this poisson mix model
	double bscLogLikelihood;//BSC logLikelihood of this poisson mix model
	MixtureModel mmopc;// Expectation Maximization Mixture Model
	MixtureModel mopbsc;;//Poisson Bregman soft clustering Mixture Model
	/**
	 * Main function.
	 * @param args
	 */
	public PoissonDataFitter (PVector[] points,int n ) {//fit the datapoints to a mixture of n poissons 

		// Display
		String title = "";
		title += "+---------------------------------------------+\n";
		title += "| EM Poisson Fitter / Bergman Soft Clustering  \n";
		title += "+---------------------------------------------+\n";
		System.out.print(title);

		// Variables
		Vector<PVector>[] clusters = KMeans.run(points, n);

		// Classical EM
		mmopc = ExpectationMaximization1D.initialize(clusters);
		mmopc = ExpectationMaximization1D.run(points, mmopc);
		emLogLikelihood=mmopc.getEMLogLikelihod();

		// Bregman soft clustering
		mopbsc = BregmanSoftClustering.initialize(clusters, new Poisson());
		mopbsc = BregmanSoftClustering.run(points, mopbsc);
		bscLogLikelihood=mopbsc.getBSCLogLikelihod();
		//System.out.println("Mixure model estimated using Bregman soft clustering \n" + mmef + "\n");
		System.out.println("  Poisson Fitter end");
	}

	public double getEMLogLikelihood(){
		return emLogLikelihood;
	}
	
	public MixtureModel getEMmodel(){
		return mmopc;
	}
	
	public double getBSCLogLikelihood(){
		return bscLogLikelihood;
	}
	
	
	public MixtureModel getBSCModel(){
		return mopbsc;
	}
}
