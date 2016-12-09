


import jMEF.MixtureModel;
import jMEF.PVector;;



public class PoissonMixturePDF {
	double [] yDataPoints;
	double [] xDataPoints;
	double [] lambdas;

	double beg;//bgining of x axis
	double end;//end of x axis
	double step;//step on x axis
    double maxYvalue=0.0;
    MixtureModel mixtMod;
    
	public PoissonMixturePDF(MixtureModel mm,double beg,double end){
		step=((end-beg)/200);//    /500);
		mixtMod=mm;
		System.out.println("PoissonMixturePDF beg:"+beg+" end:"+end+" step:"+step+" xdatapoint size:"+((int) ((end-beg)/step))+" mixturemodel:"+mm.toString());
		lambdas=new double[mm.size];
		//sigmas=new double[mm.size];
		this.beg=beg;
		this.end=end;
	
		yDataPoints=new double[(int) ((end-beg)/step)];
		xDataPoints=new double[yDataPoints.length];
		
	
		for (int p=0;p<mm.size;p++){//for each param extract mu and sigma
			lambdas[p]=((jMEF.PVector)mm.param[p]).array[0];	
			//sigmas[p]=((jMEF.PVector)mm.param[p]).array[1];	
		}
		double currentY;
		
		double dp=beg;//datapoint
		int dpInd=0;
System.out.print("        DATAPOINTS: [");
		while (dp<(end-step)){//for each datapoint
			currentY=0.0;
			
			for(int mixtElem=0;mixtElem<mm.size;mixtElem++){//for each mixture member
				currentY+=mm.weight[mixtElem]*pdf(dp,lambdas[mixtElem]);
				if( dp<3)System.out.print(" w"+mixtElem+":"+mm.weight[mixtElem]);
			}
			//if(dp>30 && dp<35)
			
			xDataPoints[dpInd]=dp;
			yDataPoints[dpInd++]=currentY;
			if (currentY>maxYvalue)maxYvalue=currentY;
			
			System.out.print("("+dp+" ,"+currentY+" ); ");
			
			dp+=step;
		}
		System.out.println();
		System.out.println("PoissonMixturePDF END");

	}
	
	
	// return pdf(x) = standard Gaussian pdf
	/*
    public static double pdf(double x) {
        return Math.exp(-x*x / 2) / Math.sqrt(2 * Math.PI);
    }
	 */
	public static double factorial ( double input )
	{
	  double x, fact = 1;
	  for ( x = input; x > 1; x--)
	     fact *= x;

	  return fact;

	}
    // return pdf(x, mu) = Poisson pdf with mean and stddev lambda
    public static double pdf(double x, double lambda) {
    	
    	double denominator=factorial(x);
    	
    	double numerator=Math.pow(lambda,x)*Math.exp(-lambda);
    	
    	if( x<3)System.out.print(" {x:"+x+" l:"+lambda+" n:"+numerator+"/ d:"+denominator+"="+(numerator / denominator)+"} ");
        
    	return (numerator / denominator);
    }






 

}

