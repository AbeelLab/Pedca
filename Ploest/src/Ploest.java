import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.UnsupportedEncodingException;
import java.util.ArrayList;

import dataFitters.GaussianDataFitter;
public class Ploest {
/*	BASECLEAR
	-p PloestBaseClear2000 -i C:\Users\Mel\Documents\BIOINFORMATICS\DELFT_Research\Data\PastorianusCBS_1483\BAM\sorted_CBS1483Pastorianus.bam
-o C:\Users\Mel\Documents\BIOINFORMATICS\DELFT_Research\Data\PastorianusCBS_1483 -v C:\Users\Mel\Documents\BIOINFORMATICS\DELFT_Research\Data\PastorianusCBS_1483\VCF\PilonPastorianusCBS.vcf
*/	
/*	SIMULATEDDATASET

	-p PloestSim500All -i C:\Users\Mel\Documents\BIOINFORMATICS\DELFT_Research\Data\Simulated\BAMs\AllSimLibs.bam
-o C:\Users\Mel\Documents\BIOINFORMATICS\DELFT_Research\Data\Simulated\ -v C:\Users\Mel\Documents\BIOINFORMATICS\DELFT_Research\Data\Simulated\PhasingAndVCFs\AllChroms\Pilon200_500Sim.vcf


	static String inputFile ="C://Users//Mel//Documents//BIOINFORMATICS//DELFT_Research//Data//CBS_Novogene//sorted_CBS1483Novogene.bam";
	static String outputFile ="C://Users//Mel//Documents//BIOINFORMATICS//DELFT_Research//Data//SimulatedReads//CBS_Novogene";
	static File fin=new File("C://Users//Mel//Documents//BIOINFORMATICS//DELFT_Research//Data//SimulatedReads//samFiles");
*/	
	static String projectName="PloEstCBS1483_wlenght2000";	
	static String inputFile ="C://Users//Mel//Documents//BIOINFORMATICS//DELFT_Research//Data//CBS_1483//sorted_CBS1483.bam";
	static String outputFile ="C://Users//Mel//Documents//BIOINFORMATICS//DELFT_Research//Data//SimulatedReads//CBS_1483";
	static File fin=new File("C://Users//Mel//Documents//BIOINFORMATICS//DELFT_Research//Data//CBS_1483//sorted_CBS1483.bam");

	static File vcfFile;
	static boolean baseCallIsOn=false;
	static int forceK=0;//deprecated (for gaussian and poisson fitting)
	static int windowLength=500;
	static int k=49;//default 49 , mode smoother window
	static File currentFolder;
	static int COV_RATE=100;
	static double SIGNIFICANT_MIN=0.01;//threshold to calculate the minimal points in the gaussian mixture plot to be considered significant (in PloestPlotter.significantMin)
	static int nbOfRuns=100;//nb of runs of
	
	
	
	public static void main(String[] args) {
		long startTime = System.currentTimeMillis();
		//printHelp();
		int argsIndex[]=new int[10];  //0:-p,projectName 
								     //1:-i,inputFile (bam)
								     //2:-s,SIGNIFICANT_MIN
								     //3:-w,windowLength
								     //4:-c,COV_RATE
								     //5:-r,numberOfRuns//deprecated
								     //6:-o,outputFile
									 //7:-k,mode average window
									 //8:-m,multiple files is set
									 //9:-v, vcf input file
							
	
		if (args.length > 0) {
			System.out.print(" Ploest "); 
			for(int i=0;i<args.length;i++){
				System.out.print(" " +args[i]+" "+args[++i]);
			} 
			System.out.println();
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
					case "-v":argsIndex[9]=i+1; break;
					default: break;
					}					
				}
			
				try {
										
					outputFile = args[argsIndex[6]];
					
					if(argsIndex[2]!=0)	SIGNIFICANT_MIN=Double.parseDouble(args[argsIndex[2]]);
					if(argsIndex[3]!=0)	windowLength=Integer.parseInt(args[argsIndex[3]]);
					projectName=args[argsIndex[0]]+windowLength;	
					if(argsIndex[4]!=0)COV_RATE=Integer.parseInt(args[argsIndex[4]]);//bigger, more detail . Default 10
					if(argsIndex[5]!=0)nbOfRuns=Integer.parseInt(args[argsIndex[5]]);
					if(argsIndex[7]!=0)k=Integer.parseInt(args[argsIndex[7]]);
					if (argsIndex[8]==0 ){
						inputFile=args[ argsIndex[1]];
					}else{
						fin=new File(args[argsIndex[8]]);
					}
					if(argsIndex[9]!=0){
						vcfFile=new File(args[argsIndex[9]]);
						baseCallIsOn=true;
					}
				} catch (Error e) {
					
					System.out.println("Could not read input arguments");
					printHelp();
					e.printStackTrace();
				} 
		

			}
		}
		//currentFolder =	 new File(outputFile + "//" + projectName+ "//Contig_Coverage_Charts//");
		//currentFolder.mkdirs();
		currentFolder =	 new File(outputFile + "//" + projectName+ "//Ploidy_Estimation_Charts//");
		currentFolder.mkdirs();

		SamParser bp=null;
		//COV_RATE=windowLength/5;
		System.out.println("COV_RATE="+COV_RATE);
		try {
			if (args.length > 0){
				if(argsIndex[8]==0 ){
					bp=new SamParser(inputFile,outputFile);
				
				}else{
					bp=new SamParser(fin,outputFile);

				}
			}else bp=new SamParser(inputFile,outputFile);// IS THIS RIGHT??? (AFTER THE ELSE)
			
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

		long endTime   = System.currentTimeMillis();
		long totalTime = endTime - startTime;
		System.out.println("TOTAL TIME : "+totalTime/1000);
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
