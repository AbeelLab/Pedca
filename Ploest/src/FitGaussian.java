import org.jfree.data.xy.XYDataset;

public class FitGaussian {
	double[][]  data;
	static int MAXNBMODELS=8;
	MultivariateNormalMixtureExpectationMaximization mixture;
	public FitGaussian(double[][]  d){
		data=d;
		System.out.println("data.length "+ data.length);
		System.out.println("data[1].length "+ data[1].length );
		mixture=new MultivariateNormalMixtureExpectationMaximization(data);
		
		for(int nbModels=2;nbModels<MAXNBMODELS;nbModels++){
			mixture.fit(mixture.estimate(data,nbModels));
			System.out.println("loglikelihood of "+ nbModels+" ="+mixture.getLogLikelihood());
			//mixture.getFittedModel();
		}
		
		
		
	}
	
}
