import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.UnsupportedEncodingException;
import java.util.ArrayList;

import dataFitters.GaussianDataFitter;
public class Ploest {
/*	
	static String projectName="CBS_Novogene_";
	static String inputFile ="C://Users//Mel//Documents//BIOINFORMATICS//DELFT_Research//Data//CBS_Novogene//sorted_CBS1483Novogene.bam";
	static String outputFile ="C://Users//Mel//Documents//BIOINFORMATICS//DELFT_Research//Data//CBS_Novogene";
*/	

	static String projectName="PloEst200";
	static String inputFile ="C://Users//Mel//Documents//BIOINFORMATICS//DELFT_Research//Data//SimulatedReads//200bp.sam";
	static String outputFile ="C://Users//Mel//Documents//BIOINFORMATICS//DELFT_Research//Data//SimulatedReads//PloEst";
	static File fin=new File("C://Users//Mel//Documents//BIOINFORMATICS//DELFT_Research//Data//SimulatedReads//samFiles");
	static int windowLength=500;
	static File currentFolder;
	static int COV_RATE=10;
	static double SIGNIFICANT_MIN=0.01;//threshold to calculate the minimal points in the gaussian mixture plot to be considered significant (in PloestPlotter.significantMin)
	static int forceK=0;//optional: force this nb Of Gaussian Mixtures
	static int nbOfRuns=100;//nb of runs of
	
	
	
	public static void main(String[] args) {
		int argsIndex[]=new int[10];  //0:-p,projectName 
								     //1:-i,inputFile
								     //2:-s,SIGNIFICANT_MIN
								     //3:-w,windowLength
								     //4:-c,COV_RATE
								     //5:-r,numberOfRuns
								     //6:-o,outputFile
									 //7:-k,force nbOfMixtures
									 //8:-m,multiple files is set
	
		if (args.length > 0) {
			//-help
			if ((args[0]).equals("-h") || (args[0]).equals("-help") || (args[0]).equals("help") || (args[0]).equals("help")){
				printHelp();
			}else {		
				for(int i=0;i<args.length;i++){
				String arg=args[i];
					switch(arg){
					case "-p": argsIndex[0]=i+1; break;
					case "-i": argsIndex[1]=i+1; break;
					case "-s": argsIndex[2]=i+1; break;
					case "-w": argsIndex[3]=i+1; break;
					case "-c": argsIndex[4]=i+1; break;
					case "-r": argsIndex[5]=i+1; break;
					case "-o": argsIndex[6]=i+1; break;
					case "-k": argsIndex[7]=i+1; break;
					case "-m":argsIndex[8]=i+1; break;
					default: break;
					}					
				}
			
				try {
					projectName=args[argsIndex[0]];						
					outputFile = args[argsIndex[6]];
					
					if(argsIndex[2]!=0)	SIGNIFICANT_MIN=Double.parseDouble(args[argsIndex[2]]);
					if(argsIndex[3]!=0)	windowLength=Integer.parseInt(args[argsIndex[3]]);
					if(argsIndex[4]!=0)COV_RATE=Integer.parseInt(args[argsIndex[4]]);//bigger, more detail . Default 10
					if(argsIndex[5]!=0)nbOfRuns=Integer.parseInt(args[argsIndex[5]]);
					if(argsIndex[7]!=0)forceK=Integer.parseInt(args[argsIndex[7]]);
					if (argsIndex[8]==0 ){
						inputFile=args[ argsIndex[1]];
					}else{
						fin=new File(args[argsIndex[8]]);
					}
				} catch (Error e) {
					
					System.out.println("Could not read input arguments");
					printHelp();
					e.printStackTrace();
				} 
			

			}
		}
		currentFolder =	 new File(outputFile + "//" + projectName+ "//Contig_Coverage_Charts//");
		currentFolder.mkdirs();
		currentFolder =	 new File(outputFile + "//" + projectName+ "//Ploidy_Estimation_Charts//");
		currentFolder.mkdirs();

		SamParser bp=null;
		try {
			if (args.length > 0){
				if(argsIndex[8]==0 ){
					bp=new SamParser(inputFile,outputFile);
				}else{
					bp=new SamParser(fin,outputFile);
				}
			}
			bp=new SamParser(inputFile,outputFile);
		} catch (FileNotFoundException e) {
			System.out.println("FileNotFoundException");
			e.printStackTrace();
		} catch (UnsupportedEncodingException e) {
			System.out.println("UnsupportedEncodingException");
			e.printStackTrace();
		} catch (IOException e) {
			System.out.println("IOException with fin samparser");
			e.printStackTrace();
		}


	}
	//0:-p,projectName 
    //1:-i,inputFile
    //2:-s,SIGNIFICANT_MIN
    //3:-w,windowLength
    //4:-c,COV_RATE
    //5:-r,numberOfRuns
    //6:-o,outputFile
	 //7:-k,force nbOfMixtures

	public static void printHelp() {
		System.out.println("\n +++++++++++++++++++++++++++++++++++++++++++++++++++++++");
		System.out.println("\nPloEst -help:");
		System.out.println(
				"\nUSAGE:    java –jar Ploest.jar  -p <projectname> -i < inputFile> -c <coverage sampling rate> -w <windowlength> -s <significant min threshold> -k <forced number of feautures>  -r <NbOfRuns> -o <outputFolder> ");
		System.out.println("\n(projectname)         – (String)  Prefix used to generate the results file.");
		System.out.println(
				"-i (inputfile)          - (String)   Pathway (absolute if necessary) to the input file containing the original sequence or contig list.");
		System.out.println("                                  Must be in fasta or .txt file");
		System.out.println("-m (multiple libraries)   - (optional) if this option is selected, the input will be the path to a .txt file containing all sam/bam files to be used instead of input");

		System.out.println(
				"-o (outputFolder)       - (String)   Pathway (absolute if necessary) to the autogenerated output folder that will contain the results('./' if the folder points to the current directory)");
		System.out.println(
				"-w (windows length)       - (int)      Length of the sliding window, to measure the coverage inside  contig. Default 500 bp");
		System.out.println(
				"-c (coverage rate)        - (int)    Rate to average the coverage measured by the sliding window");

		System.out.println("-k (k mixtures)        - (optional int)      Runs PloEst forcing k number of gaussians mixtures to fit the reads distribution");
		System.out.println("-s (significant min)   - (optional double) Threshold to consider a minimum point in the gausian fit plot to be sigificant. Default 0.01");

		System.out.println(
				"-r (Number of runs)   - (int)  Runs Gaussian mixture data fitter this number of times before averaging the results. Default 100");
		
		
		System.exit(1);

		//*/
	}

}
