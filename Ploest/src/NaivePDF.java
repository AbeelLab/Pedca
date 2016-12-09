


import jMEF.MixtureModel;
import jMEF.PVector;;



public class NaivePDF {
	double [] yDataPoints;
	double [] xDataPoints;


	int[] readCounts;

	double beg;//bgining of x axis
	double end;//end of x axis
	double step;//step on x axis
    double maxYvalue=0.0;
    MixtureModel mixtMod;
    
	public NaivePDF(int[] rc){
		
		System.out.print("naivepdf NON SOFT rc =[");
		for (int p=0;p<rc.length;p++){//for each param extract mu and sigma

			System.out.print(p+" "+rc[p]+" ;");
		}
		System.out.println("]");
		
		int smootherLength=35;//must be odd number
		int smootherWing=smootherLength/2;//part of the smootherLength before or after the pointer
		readCounts=rc;
		System.out.println("NaivePDF smootherLength:"+smootherLength+" smootherWing:"+smootherWing+ " rc size"+rc.length+" [");
	

	
		yDataPoints=new double[readCounts.length];
		xDataPoints=new double[readCounts.length];
	
		for (int p=0;p<readCounts.length;p++){//for each point
			double sum=0;
			int substract=0;
			for(int cp=(p-smootherWing);cp<(p+smootherWing);cp++){//average over the smoother window
				
				if(cp>0 && cp<readCounts.length){
					sum=readCounts[cp];
				}else {
					substract++;	
					System.out.print(" cp"+cp+" sum:"+sum+ " subst:"+substract+" ");
				}
			}
			
			for(int cp=(p-smootherWing);cp<(p+smootherWing);cp++){//
				
				if(cp>0 && cp<readCounts.length){
					xDataPoints[cp]=(double)cp;
					yDataPoints[cp]=sum/(smootherLength-substract);
					
				}
			}
			p=p+smootherWing;
			
			//System.out.print(xDataPoints[p]+" "+yDataPoints[p]+" ;");
		}
		System.out.println();
		
		
		System.out.print("naivepdf SOFT readcounts =[");
		for (int p=0;p<readCounts.length;p++){//for each param extract mu and sigma
			System.out.print(xDataPoints[p]+" "+yDataPoints[p]+" ;");
		}
		System.out.println("]");
		/*
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
		 */
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

