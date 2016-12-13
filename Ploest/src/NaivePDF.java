


import jMEF.MixtureModel;
import jMEF.PVector;;



public class NaivePDF {
	double [] yDataPoints;
	double [] xDataPoints;


	float[] readCounts;

	double beg;//bgining of x axis
	double end;//end of x axis
	double step=1;//step on x axis
    double maxYvalue=0.0;
    double maxXvalue=0.0;
    MixtureModel mixtMod;
    static int smootherLength;
    static int smootherWing;
    
	public NaivePDF(float[] rc){
		
		System.out.print("naivepdf NON SOFT rc =[");
		for (int p=0;p<rc.length;p++){//for each param extract mu and sigma

			System.out.print(p+" "+rc[p]+" ;");
		}
		System.out.println("]");
		readCounts=rc;
		maxXvalue=readCounts.length;
		
		smootherLength=(int) (maxXvalue/(3*NaivePloestPlotter.MAX_NB_MIXTURES));//the window of our smoother must be able to discretize over at least 10 different clusters (MAX_NB_MIXTURES). The minimum length should be 2X. We go for a safer 3X
		if((smootherLength & 1) == 0   ){
			System.out.println("smoother Length = "+(smootherLength+1)+" instead of "+smootherLength);
			smootherLength++;//must be odd number
		}
		smootherWing=smootherLength/2;//part of the smootherLength before or after the pointer
		
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
					if(yDataPoints[cp]>maxYvalue)maxYvalue=yDataPoints[cp];
					
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
		

	
	}
	
	

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

