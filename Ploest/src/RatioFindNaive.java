import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.io.UnsupportedEncodingException;
import java.math.RoundingMode;
import java.text.DecimalFormat;
import java.util.Stack;

public class RatioFindNaive
{
	static int MAX_NB_MIXTURES=10;
	double[] ds;
	DecimalFormat df = new DecimalFormat("#.##");
	CNVscore[] scores;
	CNVscore bestScore=null;
	
	int consensus;//% of consensus in corrected results (certainty of this prediction)


	
	
	public RatioFindNaive(double[] ds){
		
		df.setRoundingMode(RoundingMode.CEILING);
	
		
		this.ds=ds;
	
		scores=new CNVscore[(MAX_NB_MIXTURES+1-ds.length)];//vector with all the scores for each of the posible ratio solutions

		double candUnit=0.0;
		double product=00;
		double[] productsVect;
		for (int i=(MAX_NB_MIXTURES+1-ds.length);i>0;i--){//TODO: I'would rewrite this section to start evaluating the candidate unit 
															//from the bigger cluster, instead of the smallest, which will give a higher 
															//accuracy estimation of candUnit
			candUnit=ds[0]/i;
			
			productsVect=new double[MAX_NB_MIXTURES];
			for (int j=i+1;j<MAX_NB_MIXTURES;j++){//for all candidate Alternative Nb Of Mixtures higher than candUnit
				product=candUnit*j;//get product ( theoretical mean value) = candUnit * candAlternativeNbOfMixture
				productsVect[j]=product;//store all products for all candAlternativeNbOfMixture
			}
			scores[i-1]=nextNearestMeans(productsVect,i);
		}
		bestScore=findMinScore(scores);
	
	}
	
	
	private CNVscore nextNearestMeans(double[] productsVect, int ratio) {
		CNVscore result = null;
		int[] resultInds=new int[productsVect.length];
		double[]resultDists=new double[productsVect.length];
		double minDist;
		int minInd=0;
		double candDist;
		for (int i=0;i<productsVect.length;i++){//for all possible products (candUnit* other nb of mixtures higher than THIS ratio)
			if (productsVect[i]!=0.0){//start with THIS value
				minDist=productsVect[i]*MAX_NB_MIXTURES;
				for (int nd=1;nd<ds.length;nd++){
					candDist=round(100*Math.abs(ds[nd]-productsVect[i])/SamParser.maxWindows,3);//get distance to all next means values (divided by maxWindows to normalize)
					if (candDist<minDist){
						minDist=candDist;
						minInd=nd;

						resultDists[i]=minDist;
						resultInds[i]=minInd;
					}
				}

			}

		}
		result=new CNVscore(resultInds, resultDists, ratio);
		return result;	
	}
	
	
	public static double round(double value, int places) {
	    if (places < 0) throw new IllegalArgumentException();

	    long factor = (long) Math.pow(10, places);
	    value = value * factor;
	    long tmp = Math.round(value);
	    return (double) tmp / factor;
	}

	public class CNVscore{
		boolean respectsMaxNbOfMixtures=true;//validates this score if CN of secondary ds (d1,d2,etc...) are > 0
		double [] bestMinDistances=new double[ds.length];
		int [] bestCNVIndexes=new int[ds.length];
		double score=0.0;
		final double candidateUnit;

		private CNVscore(int [] cnvi,double [] md,int ratio){
			//initialize vectors
			bestMinDistances=new double[ds.length];
			bestCNVIndexes=new int[ds.length];
			for(int i=0;i<ds.length;i++){
				bestMinDistances[i]=findMax(md);//initialize with max distance
			}

			bestCNVIndexes[0]=ratio;
			bestMinDistances[0]=0.0;
			candidateUnit=ds[0]/ratio;
			computeScore(cnvi,md);
			/*
		System.out.print(" cnv[");
		for(int i=0;i<md.length;i++){
			System.out.print(cnvi[i]+"   ");
		}
		System.out.println("]");
		System.out.print("  md[");
		for(int i=0;i<md.length;i++){
			System.out.print(df.format(md[i])+" ");
		}
		System.out.println("]");
			 */

		}

		public double getCandidateUnit(){
			return candidateUnit;
		}
		
		private double findMax(double[] md) {
			double max=0.0;
			for(int i=0;i<md.length;i++){
				if(md[i]>max)max=md[i];
			}
			return max;
		}

		private void computeScore(int [] cnvi,double [] md) {

			//find for each remaining gaussian, the minimum index
			for(int d=1;d<ds.length;d++){//for each of the remaining gaussians d
				//System.out.println("== bestMinDistances length:"+bestMinDistances.length+" d:"+d+" md.length:"+md.length );
				for (int cnv=0;cnv<cnvi.length;cnv++){//find its minimum distance 
					if(cnvi[cnv]==d && md[cnv]<bestMinDistances[d]){//if it's the right cnv and is a new min
						bestMinDistances[d]=md[cnv];//change min
						bestCNVIndexes[d]=cnv;//store its index
					}
				}

			}
			//compute score
			for(int d=0;d<ds.length;d++){
				if (d>0 ){
					score+=bestMinDistances[d];
					if (bestCNVIndexes[d]==0)respectsMaxNbOfMixtures=false;
				}
			}
			
		}
		
		
		public void printScore(){
			System.out.println("========================");
			for(int d=0;d<ds.length;d++){
				System.out.println("gaussian cluster number:"+d+" Copy number estimation:"+bestCNVIndexes[d]+" Distance error:"+bestMinDistances[d]);
			}
			System.out.println();
			System.out.println("Estimation score:"+score);
			System.out.println("Maximum Nb Of Mixtures respected = "+respectsMaxNbOfMixtures);
			System.out.println("========================");
		}
	}

	private CNVscore findMinScore(CNVscore[] scs) {
		CNVscore result=null;
		int indMin=scs.length+1;//initialization with a false control value
		
		double min = ds[ds.length-1];
		boolean hasValidScore=false;
		for(int i=0;i<scs.length;i++){
			if(scs[i].respectsMaxNbOfMixtures){
				if (!hasValidScore){//is first valid score
					hasValidScore=true;
					indMin=i;//initialize indMin with a real value;
					min=scs[i].score;//and min
				}
				if (scs[i].score<min){
					min=scs[i].score;
					indMin=i;
				}
			}
		}
		if(indMin!=scs.length+1){
			result=scs[indMin];
			System.out.println("%%%%%%%%% BEST SCORE %%%%%%%%%%%%");
			result.printScore();
		}else{
			System.out.println("No CN mixture was able to satisfy the constraints. Result == null");
		}
		return result;

	}
	
	public void writeOut() {
		PrintWriter writer=null ;
		try {
			writer = new PrintWriter(Ploest.outputFile + "//" + Ploest.projectName+ "//"+Ploest.projectName+"PloidyEstimation.txt", "UTF-8");
			writer.println(">FINAL NUMBER OF CLUSTERS: "+NaivePloestPlotter.clusterMus.length );
			writer.println("\n");
	
			for (int g=0;g<NaivePloestPlotter.clusterMus.length;g++){
				writer.println("cluster center ["+g+"]:\t"+NaivePloestPlotter.clusterMus[g]);		
			}
			writer.println("\n");
			
			writer.println("> PLOIDY ESTIMATION :"+ Ploest.projectName);
			writer.println("> CLUSTER_NUMBER \tCOPY_NUMBER_ESTIMATION \t DISTANCE_ERROR (% MAX reads counts) ");
			for(int d=0;d<ds.length;d++){
				writer.println("\t\t"+d+" \t\t\t"+bestScore.bestCNVIndexes[d]+" \t\t\t"+bestScore.bestMinDistances[d]);
			}
	
			writer.println("\n");
			writer.println(">Estimation distance score: "+bestScore.score);
			//writer.println(">Estimation consensus: "+consensus+" %");
			writer.println(">Maximum Nb Of Mixtures respected = "+bestScore.respectsMaxNbOfMixtures);
			writer.close();
		} catch (FileNotFoundException | UnsupportedEncodingException e) {
			if(writer!=null)writer.close();
			// TODO Auto-generated catch block
			e.printStackTrace();
		}		
	}
	
	public void writeOutPoisson() {
		try {
			PrintWriter writer = new PrintWriter(Ploest.outputFile + "//" + Ploest.projectName+ "//"+Ploest.projectName+"PloidyEstimation.txt", "UTF-8");
			
			
			writer.println(">FINAL NUMBER OF MIXTURES: "+PloestPlotter.finalNumberOfMixtures+" poissons\n" );
			writer.println("\n");
			
			for (int g=0;g<PloestPlotter.pdfResult.getBSCModel().param.length;g++){
				writer.println("weight["+g+"]:"+PloestPlotter.gMMweights[g]+" \tmu["+g+"]"+PloestPlotter.gMMmus[g]);		
			}
				
			writer.println("\n");
			
			writer.println("> PLOIDY ESTIMATION :"+ Ploest.projectName);
			writer.println("> POISSON_CLUSTER_NUMBER \tCOPY_NUMBER_ESTIMATION \t DISTANCE_ERROR ");
			for(int d=0;d<ds.length;d++){
				writer.println("\t\t"+d+" \t\t\t"+bestScore.bestCNVIndexes[d]+" \t\t\t"+bestScore.bestMinDistances[d]);
			}
			writer.println("\n");
			writer.println("Estimation distance score: "+bestScore.score);
			writer.println("Estimation consensus: "+consensus+" %");//100*correctedResults[finalNumberOfMixtures]/NbOfRuns
			writer.println("Maximum Nb Of Mixtures respected = "+bestScore.respectsMaxNbOfMixtures);
		
			writer.close();
		} catch (FileNotFoundException | UnsupportedEncodingException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

		
	}


}
