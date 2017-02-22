


import jMEF.MixtureModel;
import jMEF.PVector;;



public class NaivePDF {
	double [] yDataPoints;
	double [] xDataPoints;


	float[] readCounts;//normalized read counts

	double beg;//bgining of x axis
	double end;//end of x axis
	double step=1;//step on x axis
    double maxYFITvalue=0.0;
    double maxYHISTOGRAMvalue=0.0;
    double maxXvalue=0.0;
    int peakYvalueIndex=0;//index f the peak value (therefore in both histogram and fit). Used to find the correct ratio of both plots
    MixtureModel mixtMod;
    static int smootherLength;
    static int smootherWing;
    
	public NaivePDF(float[] rc){
		/*
		System.out.print("naivepdf NON SOFT rc =[");
		for (int p=0;p<rc.length;p++){//for each param extract mu and sigma
			System.out.print(p+" "+rc[p]+" ;");
		}System.out.println("]");
		*/
		
		readCounts=rc;
		maxXvalue=readCounts.length;
		
		smootherLength=(int) (maxXvalue/(Math.round(3*NaivePloestPlotter.MAX_NB_MIXTURES)));//the window of our smoother must be able to discretize over at least 10 different clusters (MAX_NB_MIXTURES). The minimum length should be 2X. We go for a safer 3X
		if((smootherLength & 1) == 0   ){//if even number
			smootherLength++;//must be odd number
		}
		smootherWing=smootherLength/2;//part of the smootherLength before or after the pointer
		
		yDataPoints=new double[readCounts.length];
		xDataPoints=new double[readCounts.length];
		

		for (int p=0;p<readCounts.length;p++){//for each point
			double sum=0;
			int substract=0;//substract these bins (for begining and end of genome)
			
			for(int cp=(p-smootherWing);cp<(p+smootherWing);cp++){//sum over the smoother window
				
				if(cp>0 && cp<readCounts.length){
					sum+=readCounts[cp];
				}else {
					substract++;	//but substract these points if the wind goes out of the existing range
				}
			}
			double average=sum/(smootherLength-substract);
			
			for(int cp=(p-smootherWing);cp<(p+smootherWing);cp++){//set the average value over the smoothing window	
				if(cp>0 && cp<readCounts.length){
					xDataPoints[cp]=(double)cp;
					yDataPoints[cp]=average;//remove points that were outside of range
					if(yDataPoints[cp]>maxYFITvalue){
						maxYFITvalue=yDataPoints[cp];
					}
				}
			}
			p=p+smootherWing-1;//move p to next bin start point

		}
	
		/*
		System.out.println(" smootherLength:"+smootherLength);
		System.out.print("naivepdf SOFT readcounts =[");
		for (int p=0;p<readCounts.length;p++){//for each param extract mu and sigma
			System.out.print(xDataPoints[p]+" "+yDataPoints[p]+" ;");
		}System.out.println("]");
		*/
	}
	
	

	public static double factorial ( double input )
	{
	  double x, fact = 1;
	  for ( x = input; x > 1; x--)
	     fact *= x;

	  return fact;

	}
    // return distance to the center mu of this cluster 
    public static double pdf(double x, double mu) {
    	//if(x<150)System.out.println("x:"+x+" m:"+mu+" d:"+Math.abs(x-mu));
    	return Math.abs(x-mu);
    }






 

}

