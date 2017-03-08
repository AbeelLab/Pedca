import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.UnsupportedEncodingException;
import java.util.ArrayList;

import dataFitters.GaussianDataFitter;

public class Ploest {
	/*
	 *
	 * SIMULATEDDATASET
	 * 
	 * -p PloestSim500All -i
	 * C:\Users\Mel\Documents\BIOINFORMATICS\DELFT_Research\Data\Simulated\BAMs\
	 * AllSimLibs.bam -o
	 * C:\Users\Mel\Documents\BIOINFORMATICS\DELFT_Research\Data\Simulated\ -v
	 * C:\Users\Mel\Documents\BIOINFORMATICS\DELFT_Research\Data\Simulated\
	 * PhasingAndVCFs\AllChroms\Pilon200_500Sim.vcf
	 * 
	 * 
	 *  BASECLEAR -p PloestBaseClear2000 -i
	 * C:\Users\Mel\Documents\BIOINFORMATICS\DELFT_Research\Data\
	 * PastorianusCBS_1483\BAM\sorted_CBS1483Pastorianus.bam -o
	 * C:\Users\Mel\Documents\BIOINFORMATICS\DELFT_Research\Data\
	 * PastorianusCBS_1483 -v
	 * C:\Users\Mel\Documents\BIOINFORMATICS\DELFT_Research\Data\
	 * PastorianusCBS_1483\VCF\PilonPastorianusCBS.vcf
	 *
	 * args[0] = "-p";
				args[1] = "PloestBaseClear" ;	
				args[2] = "-w";
				args[3] = winLengths[wlInd];
				args[4] = "-i";
				args[5] = "C://Users//Mel//Documents//BIOINFORMATICS//DELFT_Research//Data//PastorianusCBS_1483//BAM//sorted_CBS1483Pastorianus.bam";
				args[6] = "-o";
				args[7] = "C://Users//Mel//Documents//BIOINFORMATICS//DELFT_Research//Data//PastorianusCBS_1483";
				args[8] = "-v";
				args[9] = "C://Users//Mel//Documents//BIOINFORMATICS//DELFT_Research//Data//PastorianusCBS_1483//VCF//PilonPastorianusCBS.vcf";
			
	 * 
	 * NOVOGENE
	
	Novogene to CBS1483_59contigs
	
				SIGNIFICANT_MIN=0.03;
				args[0] = "-p";
				args[1] = "PloestNovogene" ;////args[1] = "NovogeneVdB2anc" ;
				//args[1] = "PedcaNovogene2ancestors" ;	
				
				args[2] = "-w";
				args[3] = winLengths[wlInd];
				args[4] = "-i";
				//args[5] = "C://Users//Mel//Documents//BIOINFORMATICS//DELFT_Research//Data//CBS_Novogene//BAM//sorted_nov2anc.bam";
				//args[5] = "C://Users//Mel//Documents//BIOINFORMATICS//DELFT_Research//Data//CBS_Novogene//BAM//sorted_NovoVdBNewMap2anc.bam";
				args[5] = "C://Users//Mel//Documents//BIOINFORMATICS//DELFT_Research//Data//CBS_Novogene//BAM//sorted_CBSNovogene.bam";

				args[6] = "-o";
				//args[7] = "C://Users//Mel//Documents//BIOINFORMATICS//DELFT_Research//Data//CBS_Novogene//PedcaNovogene2ancestors";
				args[7] = "C://Users//Mel//Documents//BIOINFORMATICS//DELFT_Research//Data//CBS_Novogene//PedcaNovogene";
				args[8] = "-v";
				args[9] = "C://Users//Mel//Documents//BIOINFORMATICS//DELFT_Research//Data//CBS_Novogene//Pilon//PilonCBSNovogene.vcf";
	 * 
	Novogene to bothAncestors(eubayanus+cerevisiae)
	
				sigMinRate=0.1;
				args[0] = "-p";
				//args[1] = "PloestNovogene" ;////args[1] = "NovogeneVdB2anc" ;
				args[1] = "PedcaNovogene2ancestors" ;	
				args[2] = "-w";
				args[3] = winLengths[wlInd];
				args[4] = "-i";
				args[5] = "C://Users//Mel//Documents//BIOINFORMATICS//DELFT_Research//Data//CBS_Novogene//BAM//sorted_nov2anc.bam";
				//args[5] = "C://Users//Mel//Documents//BIOINFORMATICS//DELFT_Research//Data//CBS_Novogene//BAM//sorted_NovoVdBNewMap2anc.bam";

				args[6] = "-o";
				args[7] = "C://Users//Mel//Documents//BIOINFORMATICS//DELFT_Research//Data//CBS_Novogene";
	
	 * 
	 * BASECLEAR
	 * 
	 * args[0] = "-p";
				args[1] = "PloestBaseClear" ;	
				args[2] = "-w";
				args[3] = winLengths[wlInd];
				args[4] = "-i";
				args[5] = "C://Users//Mel//Documents//BIOINFORMATICS//DELFT_Research//Data//CBS_Baseclear//BAM//sorted_CBS1483Pastorianus.bam";
				args[6] = "-o";
				args[7] = "C://Users//Mel//Documents//BIOINFORMATICS//DELFT_Research//Data//CBS_Baseclear";
				args[8] = "-v";
				args[9] = "C://Users//Mel//Documents//BIOINFORMATICS//DELFT_Research//Data//CBS_Baseclear//VCF//PilonPastorianusCBS.vcf";
							
	 args[0] = "-p";
				args[1] = "PedcaBaseClear59bowtiePseudo" ;	
				args[2] = "-w";
				args[3] = winLengths[wlInd];
				args[4] = "-i";
				//args[5] = "C://Users//Mel//Documents//BIOINFORMATICS//DELFT_Research//Data//CBS_Baseclear//BAM//sorted_BaseVdBNewMap2anc.bam";
				args[5] = "C://Users//Mel//Documents//BIOINFORMATICS//DELFT_Research//Data//CBS_Baseclear//BAM//sorted_CBS_bowtie_pseudo_Baseclear.bam";
				args[6] = "-o";
				args[7] = "C://Users//Mel//Documents//BIOINFORMATICS//DELFT_Research//Data//CBS_Baseclear//CBSBowtiePseudo";
				//args[8] = "-v";
				//args[9] = "C://Users//Mel//Documents//BIOINFORMATICS//DELFT_Research//Data//CBS_Baseclear//VCF//PilonPastorianusCBS.vcf";
	
	 * 
	 * CRUZI
	 * 
	  			args[0] = "-p";
				args[1] = "PedcaTCruzi_CLB_PacBio_Assembly" ;	
				args[2] = "-w";
				args[3] = winLengths[wlInd];
				args[4] = "-i";
				//args[5] = "//tudelft.net//staff-bulk//ewi//insy//DBL//mcarbajo//Genomes//Cruzi//Bams//sorted_CruziEsmeraldoLikeIllumina_full.bam";
				//args[5] = "C://Users//Mel//Documents//BIOINFORMATICS//DELFT_Research//Data//TCruzi//Bams//sorted_CruziNon_EsmeraldoLikeIllumina_full.bam";
				args[5] = "C://Users//Mel//Documents//BIOINFORMATICS//DELFT_Research//Data//TCruzi//Bams//sorted_CruziPacBioAssemblyCLB1.bam";

				args[6] = "-o";
				args[7] = "C://Users//Mel//Documents//BIOINFORMATICS//DELFT_Research//Data//TCruzi//PEDCA_CLB_PacBioAssembly";
				//args[8] = "-v";
				//args[9] = "C://Users//Mel//Documents//BIOINFORMATICS//DELFT_Research//Data//TCruzi//Pilon//PilonCruziEsmLike.vcf";

	 * 
	 */
	static String projectName = "";
	static String inputFile = "";
	static String outputFile = "";
	//static File fin = new File("C://Users//Mel//Documents//BIOINFORMATICS//DELFT_Research//Data//CBS_1483//sorted_CBS1483.bam");

	static File vcfFile;
	static boolean baseCallIsOn = false;
	static int forceK = 0;// deprecated (for gaussian and poisson fitting)
	static int windowLength = 500;
	static int k = 49;// default 49 , mode smoother window
	
	static File currentFolder;
	static int COV_RATE = 100;
	static double SIGNIFICANT_MIN = 0.05;//default =0.05 threshold to calculate the minimal
											// points in the gaussian mixture
											// plot to be considered significant
											// (in PloestPlotter.significantMin)
	static int nbOfRuns = 100;// nb of runs of

	public static void main(String[] args) {
		args=new String[10];
		long startTimeGeneral = System.currentTimeMillis();
		String[] winLengths = {"1000"};//"400","500","750","1000","2000"};//"50","100","250","400","500","750","1000","1500","2000", "3000","5000", "6500", "9000", "10000", "15000","20000", "30000", "40000", "50000","75000"};		
		for (int wlInd = 0; wlInd < winLengths.length; wlInd++) {
			try {
				
				SIGNIFICANT_MIN=0.05;
				
				args[0] = "-p";
				args[1] = "PedcaBaseclear" ;	
				args[2] = "-w";
				args[3] = winLengths[wlInd];
				args[4] = "-i";
				args[5] = "C://Users//Mel//Documents//BIOINFORMATICS//DELFT_Research//Data//CBS_Baseclear//BAM//sorted_BaseVdBNewMap2anc.bam";
				args[6] = "-o";
				args[7] = "C://Users//Mel//Documents//BIOINFORMATICS//DELFT_Research//Data//CBS_Baseclear//VdBNewMap2anc";
				args[8] = "-v";
				args[9] = "C://Users//Mel//Documents//BIOINFORMATICS//DELFT_Research//Data//CBS_Baseclear//Pilon//PilonBaseclearVdBNew2anc.vcf";
			      
			
				
				
				
				
				runPloest(args);
				

			} catch (Exception e) {
				System.err.println(args[1] + args[3]+" could not be run");
				e.printStackTrace();
			}

		}
		long endTimeGeneral = System.currentTimeMillis();
		long totalTimeGeneral = endTimeGeneral - startTimeGeneral;
		System.out.println("TOTAL TIME GENERAL ("+winLengths.length+" run/s): " + totalTimeGeneral / 1000);
	}


	public static void runPloest(String[] args) {
		
		long startTime = System.currentTimeMillis();
		baseCallIsOn = false;
		// printHelp();
		int argsIndex[] = new int[10]; // 0:-p,projectName
										// 1:-i,inputFile (bam)
										// 2:-s,SIGNIFICANT_MIN
										// 3:-w,windowLength
										// 4:-c,COV_RATE
										// 5:-r,numberOfRuns//deprecated
										// 6:-o,outputFile
										// 7:-k,mode average window
										// 8:-m,multiple files is set
										// 9:-v, vcf input file

		if (args.length > 0) {
			System.out.print(" Ploest ");
			for (int i = 0; i < args.length; i++) {
				System.out.print(" " + args[i] + " " + args[++i]);
			}
			System.out.println();
			// -help
			if ((args[0]).equals("-h") || (args[0]).equals("-help") || (args[0]).equals("help")
					|| (args[0]).equals("help")) {
				printHelp();
			} else {

				for (int i = 0; i < args.length; i++) {
					String arg = args[i];
				
					switch (arg) {
					case "-p":
						argsIndex[0] = i + 1;
						break;
					case "-i":
						argsIndex[1] = i + 1;
						break;
					case "-s":
						argsIndex[2] = i + 1;
						break;
					case "-w":
						argsIndex[3] = i + 1;
						break;
					case "-c":
						argsIndex[4] = i + 1;
						break;
					case "-r":
						argsIndex[5] = i + 1;
						break;
					case "-o":
						argsIndex[6] = i + 1;
						break;
					case "-k":
						argsIndex[7] = i + 1;
						break;
					case "-m":
						argsIndex[8] = i + 1;
						break;
					case "-v":
						argsIndex[9] = i + 1;
						//System.out.println("case -v stored not supposed to happen:");
						break;
					default:
						break;
					}
				}
				try {

					outputFile = args[argsIndex[6]];
					System.out.println("outputFile :"+outputFile);
					if (argsIndex[2] != 0)
						SIGNIFICANT_MIN = Double.parseDouble(args[argsIndex[2]]);
					if (argsIndex[3] != 0)
						windowLength = Integer.parseInt(args[argsIndex[3]]);
					projectName = args[argsIndex[0]] + windowLength;
					if (argsIndex[4] != 0)
						COV_RATE = Integer.parseInt(args[argsIndex[4]]);// bigger, more detail. Default= 10
					if (argsIndex[5] != 0)
						nbOfRuns = Integer.parseInt(args[argsIndex[5]]);
					if (argsIndex[7] != 0)
						k = Integer.parseInt(args[argsIndex[7]]);
					if (argsIndex[8] == 0) {
						inputFile = args[argsIndex[1]];
					} else {
						System.err.println("argsIndex[8] :"+argsIndex[8] + " inputFile = "+inputFile);
						//TO DO
						//implement multiple input files
						//fin = new File(args[argsIndex[8]]);
					}
					
					
					if (argsIndex[9] != 0) {
						System.out.println("vcf arg9!=0");
						vcfFile = new File(args[argsIndex[9]]);
						baseCallIsOn = true;
					}else baseCallIsOn=false;
				} catch (Error e) {

					System.out.println("Could not read input arguments");
					printHelp();
					e.printStackTrace();
				}

			}
		}

		currentFolder = new File(outputFile + "//" + projectName + "//Ploidy_Estimation_Charts//");
		currentFolder.mkdirs();
		
		SamParser bp = null;
		// COV_RATE=windowLength/5;
		System.out.println("COV_RATE=" + COV_RATE);
		try {
			if (args.length > 0) {
				if (argsIndex[8] == 0) {//one single input file
					bp = new SamParser(inputFile, outputFile);

				} else {
					System.out.println("The multiple file options is not yet implemeted. Try Plpoest again running one file at a time " );
					//TO DO
				}
			} else  bp = new SamParser(inputFile, outputFile);// IS THIS RIGHT??? (AFTER THE ELSE)

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

		long endTime = System.currentTimeMillis();
		long totalTime = endTime - startTime;
		System.out.println("TOTAL TIME : ["+args[3]+"] :"+ totalTime / 1000);
	}
	
	
	
	
	public static void printHelp() {
		System.out.println("\n +++++++++++++++++++++++++++++++++++++++++++++++++++++++");
		System.out.println("\nPloEst -help:");
		System.out.println(
				"\nUSAGE:    java –jar Ploest.jar  -p <projectname> -i < input sam/bam File> -o <outputFolder> -c <coverage sampling rate> -w <windowlength> -s <significant min threshold> -k <free>  -r <free>  \n");
		System.out.println("-p (projectname)         – (String)  Prefix used to generate the results file.");
		System.out.println(
				"-i (inputfile)          - (String)   Pathway (absolute if necessary) to the input file containing the original sequence or contig list.");
		System.out.println("                                  Must be in fasta or .txt file");
		System.out.println(
				"-m (free parameter)   - (optional) ");

		System.out.println(
				"-o (outputFolder)       - (String)   Pathway (absolute if necessary) to the autogenerated output folder that will contain the results('./' if the folder points to the current directory)");
		System.out.println(
				"-w (windows length)       - (int)      Length of the sliding window, to measure the coverage inside  contig. Default 500 bp");
		System.out.println(
				"-c (coverage rate)        - (int)    Rate to average the coverage measured by the sliding window");

		System.out.println(
				"-k (mode smoother window)        - (optional int)     number of points over which the ploidy estimation is smoothed. The mode over k numbers of windows is used to average the values of the bin. Default=49 ");
		System.out.println(
				"-s (significant min)   - (optional double) Threshold to consider a minimum point in the fit plot to be sigificant. Default 0.01");

		System.out.println(
				"-r (free parameter)   - (int)  ");

		System.exit(1);

		// */
	}

}
