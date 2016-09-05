import jMEF.MixtureModel;
import jMEF.PVector;;

/******************************************************************************
 *  Compilation:  javac Gaussian.java
 *  Execution:    java Gaussian x mu sigma
 *
 *  Function to compute the Gaussian pdf (probability density function)
 *  and the Gaussian cdf (cumulative density function)
 *
 *  % java Gaussian 820 1019 209
 *  0.17050966869132111
 *
 *  % java Gaussian 1500 1019 209
 *  0.9893164837383883
 *
 *  % java Gaussian 1500 1025 231
 *  0.9801220907365489
 *
 *  The approximation is accurate to absolute error less than 8 * 10^(-16).
 *  Reference: Evaluating the Normal Distribution by George Marsaglia.
 *  http://www.jstatsoft.org/v11/a04/paper
 *
 ******************************************************************************/

public class GaussianMixturePDF {
	double [] yDataPoints;
	double [] mus;
	double [] sigmas;
    
	public GaussianMixturePDF(MixtureModel mm,double beg,double end, double step){
		mus=new double[mm.size];
		sigmas=new double[mm.size];
		yDataPoints=new double[(int) ((end-beg)/step)];
		for (int p=0;p<mm.size;p++){//for each param extract mu and sigma
			mus[p]=((jMEF.PVector)mm.param[p]).array[0];	
			sigmas[p]=((jMEF.PVector)mm.param[p]).array[1];	
		}
		double currentY;
		double dp=beg;
		int dpInd=0;
		while (dp<end){//for each datapoint
			currentY=0.0;
			
			for(int mixtElem=0;mixtElem<mm.size;mixtElem++){//for each mixture member
				currentY+=mm.weight[mixtElem]*pdf(dp,mus[mixtElem],sigmas[mixtElem]);
			}
			yDataPoints[dpInd++]=currentY;
			dp+=step;
		}
		/*
		System.out.println();
		for (int ii=0;ii<yDataPoints.length;ii++){
			System.out.print(yDataPoints[ii]+" ");
		}
		System.out.println();
		System.out.println(" PDF y points ; length:"+yDataPoints.length+ " "+mm.size+" mixt elems");
		*/
	}
	
	
	// return pdf(x) = standard Gaussian pdf
    public static double pdf(double x) {
        return Math.exp(-x*x / 2) / Math.sqrt(2 * Math.PI);
    }

    // return pdf(x, mu, signma) = Gaussian pdf with mean mu and stddev sigma
    public static double pdf(double x, double mu, double sigma) {
        return phi((x - mu) / sigma) / sigma;
    }

    // return phi(x) = standard Gaussian pdf
    public static double phi(double x) {
        return pdf(x);
    }



    // return phi(x, mu, signma) = Gaussian pdf with mean mu and stddev sigma

    public static double phi(double x, double mu, double sigma) {
        return pdf(x, mu, sigma);
    }



 

}
