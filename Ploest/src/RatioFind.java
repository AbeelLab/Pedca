import java.math.RoundingMode;
import java.text.DecimalFormat;
import java.util.Stack;

public class RatioFind
{
	static int MAX_NB_MIXTURES=10;
	double[] ds;
	DecimalFormat df = new DecimalFormat("#.##");
	CNVscore[] scores;

	public RatioFind(double[] ds){
		df.setRoundingMode(RoundingMode.CEILING);
		this.ds=ds;
		java.util.Arrays.sort(ds);
		scores=new CNVscore[(MAX_NB_MIXTURES+1-ds.length)];

		double candUnit=0.0;
		double product=00;
		double[] productsVect;
		for (int i=(MAX_NB_MIXTURES+1-ds.length);i>0;i--){
			candUnit=ds[0]/i;
			productsVect=new double[MAX_NB_MIXTURES];

			for (int j=i+1;j<MAX_NB_MIXTURES;j++){//for all candidate Alternative Nb Of Mixtures higher than candUnit
				product=candUnit*j;//get product ( theoretical mean value) = candUnit * candAlternativeNbOfMixture
				productsVect[j]=product;//store all products for all candAlternativeNbOfMixture
			}
			scores[i-1]=nextNearestMeans(productsVect,i);	  
		}
		findMinScore(scores);

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
					candDist=Math.abs(ds[nd]-productsVect[i]);//get distance to all next means values 
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

	private class CNVscore{
		boolean respectsMaxNbOfMixtures=true;//validates this score if CN of secondary ds (d1,d2,etc...) are > 0
		int ratio;
		double [] bestMinDistances=new double[ds.length];
		int [] bestCNVIndexes=new int[ds.length];
		double score=0.0;

		private CNVscore(int [] cnvi,double [] md,int ratio){
			//initialize vectors
			bestMinDistances=new double[ds.length];
			bestCNVIndexes=new int[ds.length];
			for(int i=0;i<ds.length;i++){
				bestMinDistances[i]=findMax(md);//initialize with max distance
			}

			bestCNVIndexes[0]=ratio;
			bestMinDistances[0]=0.0;

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
				System.out.println("d:"+d+" cnv:"+bestCNVIndexes[d]+" dist:"+bestMinDistances[d]);
			}
			System.out.println();
			System.out.println("score:"+score);
			System.out.println("Respects Max Nb Of Mixtures = "+respectsMaxNbOfMixtures);
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


}