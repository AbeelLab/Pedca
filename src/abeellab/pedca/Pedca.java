package abeellab.pedca;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.UnsupportedEncodingException;
import java.util.ArrayList;
import java.util.List;

import dataFitters.GaussianDataFitter;

public class Pedca {
	/*
	 *
	 * SIMULATEDDATASET
	 * 
	 * 
	 * -p PloestSim  -i C:\Users\Mel\Documents\BIOINFORMATICS\DELFT_Research\Data\Simulated\BAMs\AllSimLibs.bam -o C:\Users\Mel\Documents\BIOINFORMATICS\DELFT_Research\Data\Simulated -v  C:\Users\Mel\Documents\BIOINFORMATICS\DELFT_Research\Data\Simulated\PhasingAndVCFs\AllChroms\Pilon200_500Sim.vcf -w 3000
	 * 
	 * 
	 *  BASECLEAR -p PedcaBaseClear59bowtiePseudo 
	 *  -i  C:\Users\Mel\Documents\BIOINFORMATICS\DELFT_Research\Data\CBS_Baseclear\BAM\sorted_CBS_bowtie_pseudo_Baseclear.bam 
	 *  -o C:\Users\Mel\Documents\BIOINFORMATICS\DELFT_Research\Data\CBS_Baseclear\CBSBowtiePseudo -m 
	 *  -v  C:\Users\Mel\Documents\BIOINFORMATICS\DELFT_Research\Data\
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
	-p PedcaNovogene59bowtiePseudo   -i  C:\Users\Mel\Documents\BIOINFORMATICS\DELFT_Research\Data\CBS_Novogene\BAM\sorted_CBS_bowtie_pseudo_Novogene.bam  -o C:\Users\Mel\Documents\BIOINFORMATICS\DELFT_Research\Data\CBS_Novogene\CBSBowtiePseudo
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
	 * -p PedcaBaseClear59bowtiePseudo -w  10000 -c 15 -b 3.5 -i C:\Users\Mel\Documents\BIOINFORMATICS\DELFT_Research\Data\CBS_Baseclear\BAM\sorted_CBS_bowtie_pseudo_Baseclear.bam
-o C:\Users\Mel\Documents\BIOINFORMATICS\DELFT_Research\Data\CBS_Baseclear\CBSBowtiePseudo  -v C:\Users\Mel\Documents\BIOINFORMATICS\DELFT_Research\Data\CBS_Baseclear\Pilon\PilonBaseclearVdBNew2anc.vcf
	 * 			args[0] = "-p";
				args[1] = "PedcaBaseclear" ;	
				args[2] = "-w";
				args[3] = winLengths[wlInd];
				args[4] = "-i";
				args[5] = "C://Users//Mel//Documents//BIOINFORMATICS//DELFT_Research//Data//CBS_Baseclear//BAM//sorted_BaseVdBNewMap2anc.bam";
				args[6] = "-o";
				args[7] = "C://Users//Mel//Documents//BIOINFORMATICS//DELFT_Research//Data//CBS_Baseclear//VdBNewMap2anc";
				args[8] = "-v";
				args[9] = "C://Users//Mel//Documents//BIOINFORMATICS//DELFT_Research//Data//CBS_Baseclear//Pilon//PilonBaseclearVdBNew2anc.vcf";
			      							
	 			args = new String[10];
		
		args[0] = "-p";
		args[1] = "PedcaBaseClear59bowtiePseudo" ;	
		args[2] = "-w";
		args[3] = "500";
		args[4] = "-i";
		args[5] = "C://Users//Mel//Documents//BIOINFORMATICS//DELFT_Research//Data//CBS_Baseclear//BAM//sorted_CBS_bowtie_pseudo_Baseclear.bam";
		args[6] = "-o";
		args[7] = "C://Users//Mel//Documents//BIOINFORMATICS//DELFT_Research//Data//CBS_Baseclear//CBSBowtiePseudo";
		args[8] = "-v";
		args[9] = "C://Users//Mel//Documents//BIOINFORMATICS//DELFT_Research//Data//CBS_Baseclear//Pilon//PilonPastorianusCBS.vcf";

	 * 
	 * CRUZI
	 * 
	 * -p PEDCA_CLB_PacBioAssembly -i C:\Users\Mel\Documents\BIOINFORMATICS\DELFT_Research\Data\TCruzi\Bams\sorted_CruziPacBioAssemblyCLB1.bam  -o C:\Users\Mel\Documents\BIOINFORMATICS\DELFT_Research\Data\TCruzi\PEDCA_CLB_PacBioAssembly -w 7700 -b 2.5 -d 0.90 -s 0.2 
	 * 
	  			SIGNIFICANT_MIN=0.05;
				
				args[0] = "-p";
				args[1] = "PedcaTCruzi_CLB_PacBio_Assembly" ;	
				args[2] = "-w";
				args[3] = winLengths[wlInd];
				args[4] = "-i";
				args[5] = "C:\Users\Mel\Documents\BIOINFORMATICS\DELFT_Research\Data\TCruzi\Bams\sorted_CruziPacBioAssemblyCLB1.bam";

				args[6] = "-o";
				args[7] = "C:\Users\Mel\Documents\BIOINFORMATICS\DELFT_Research\Data\TCruzi\PEDCA_CLB_PacBioAssembly";
				
				args[8] = "-v";
				args[9] = "C:\Users\Mel\Documents\BIOINFORMATICS\DELFT_Research\Data\TCruzi\Pilon\PilonCruziPacBio.vcf";

	 * 
	 */
	static String projectName = "";
	static String inputFile = "";
	static String outputFile = "";

	static File vcfFile;
	static boolean baseCallIsOn = false;
	static int forceK = 0;// deprecated (for gaussian and poisson fitting)
	static int windowLength = 500;
	static int MIN_WIND_LENGTH=16;
	static int k = 49;// default 49 , mode smoother window
	
	static File currentFolder;
	static int COV_RATE = 100;
	static double SIGNIFICANT_MIN = 0.1;//default =0.1 threshold to calculate the minimal
											// points in the gaussian mixture
											// plot to be considered significant
											// (in PloestPlotter.significantMin)
	static int nbOfRuns = 100;// nb of runs of
	static boolean MULTIRUN=false;
	static double BIN_FACTOR=2.5;
	static double USED_DATA=0.95;//percentage of coverage data used. Top highest values will be rejected
	static String currentProjectName="";
	public static void main(String[] args) {
		/*
		args = new String[10];
		
		args[0] = "-p";
		args[1] = "PedcaBaseClear59bowtiePseudo" ;	
		args[2] = "-w";
		args[3] = "500";
		args[4] = "-i";
		args[5] = "C://Users//Mel//Documents//BIOINFORMATICS//DELFT_Research//Data//CBS_Baseclear//BAM//sorted_CBS_bowtie_pseudo_Baseclear.bam";
		args[6] = "-o";
		args[7] = "C://Users//Mel//Documents//BIOINFORMATICS//DELFT_Research//Data//CBS_Baseclear//CBSBowtiePseudo";
		args[8] = "-v";
		args[9] = "C://Users//Mel//Documents//BIOINFORMATICS//DELFT_Research//Data//CBS_Baseclear//Pilon//PilonPastorianusCBS.vcf";

		*/
		
		runPloest(args);		
	}


	public static void runPloest(String[] args) {
		
		long startTime = System.currentTimeMillis();
		baseCallIsOn = false;
		// printHelp();
		int argsIndex[] = new int[13]; // 0:-p,projectName
										// 1:-i,inputFile (bam)
										// 2:-s,SIGNIFICANT_MIN
										// 3:-w,windowLength
										// 4:-c,COV_RATE
										// 5:-r,numberOfRuns//deprecated
										// 6:-o,outputFile
										// 7:-k,mode average window
										// 8:-m,automatic multiple windows
										// 9:-v, vcf input file
										//10:-b, number of fit bins
										//11:-d, used coverage data
		if (args.length > 0) {
			for (int i = 0; i < args.length; i++) {
				System.out.print(" " + args[i] );
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
						argsIndex[8] = i ;
						MULTIRUN=true;
						break;
					case "-v":
						argsIndex[9] = i + 1;
						break;
					case "-b":
						argsIndex[10] = i + 1;
						break;
					case "-d":
						argsIndex[11] = i + 1;
						break;
					default:
						break;
					}
				}
				if(MULTIRUN){
					MULTIRUN=false;
					List <String> newArgslist=new ArrayList <String>();
		
					System.out.println(" MULTIRUN  ");
				
					for (int n=0;n<args.length;n++){
						if(!args[n].equals("-m") && !args [n].equals("-w")){
							newArgslist.add(args[n]);
						}
						if(args [n].equals("-w")){
							n++;
						}
					}
					
					String [] newArgs=new String [newArgslist.size()];
					for (int n=0;n<newArgslist.size();n++){
						newArgs[n]=newArgslist.get(n);
					}
					
					
					runMultipleWindows( newArgs);
					return ;
					
				}
				try {
					// 0:-p,projectName
					// 1:-i,inputFile (bam)
					// 2:-s,SIGNIFICANT_MIN
					// 3:-w,windowLength
					// 4:-c,COV_RATE
					// 5:-r,numberOfRuns//deprecated
					// 6:-o,outputFile
					// 7:-k,mode average window
					// 8:-m,automatic multiple windows
					// 9:-v, vcf input file
					//10:-b, number of fit bins
					//11:-d, used coverage data
					
					outputFile = args[argsIndex[6]];
					System.out.println("outputFile :"+outputFile);
					if (argsIndex[2] != 0)
						SIGNIFICANT_MIN = Double.parseDouble(args[argsIndex[2]]);
					System.out.println("SIGNIFICANT_MIN="+SIGNIFICANT_MIN);
					if (argsIndex[3] != 0){
						windowLength = Integer.parseInt(args[argsIndex[3]]);
						if (windowLength<MIN_WIND_LENGTH)windowLength=MIN_WIND_LENGTH;
					}
						
					currentProjectName = args[argsIndex[0]] + windowLength;
					projectName = currentProjectName;
					if (argsIndex[4] != 0){
						COV_RATE = Integer.parseInt(args[argsIndex[4]]);// bigger, more detail. Default= 100
					}
					if (argsIndex[5] != 0)
						nbOfRuns = Integer.parseInt(args[argsIndex[5]]);
					if (argsIndex[7] != 0)
						k = Integer.parseInt(args[argsIndex[7]]);
					if (argsIndex[8] == 0) {
						inputFile = args[argsIndex[1]];
					} else {
						System.err.println("argsIndex[8] :"+argsIndex[8] + " inputFile = "+inputFile);

					}
					
					
					if (argsIndex[9] != 0) {
						System.out.println("ALLELE FREQUENCIES ANALYSIS CALLED");
						vcfFile = new File(args[argsIndex[9]]);
						baseCallIsOn = true;
					}else baseCallIsOn=false;
					
					if (argsIndex[10] != 0 && Double.parseDouble(args[argsIndex[10]])<4.0 && Double.parseDouble(args[argsIndex[10]])>1.0){
						BIN_FACTOR = Double.parseDouble(args[argsIndex[10]])/10;
					}else if (argsIndex[10] != 0){
						BIN_FACTOR = Double.parseDouble(args[argsIndex[10]])/10;
						if (BIN_FACTOR>4.0 || BIN_FACTOR <1.0){
							System.err.println(" BIN_FACTOR value non accepted, must be >1.0 and <4.0. Will be run with default value "+BIN_FACTOR );
							BIN_FACTOR=2.5;
						}
					}
					if (argsIndex[11] != 0){
						USED_DATA = Double.parseDouble(args[argsIndex[11]]);
					}	
					
				} catch (Error | Exception e ) {

					System.err.println("Could not read input arguments");
					printHelp();
					e.printStackTrace();
				}

			}
		}

		currentFolder = new File(outputFile + "//" + projectName + "//Ploidy_Estimation_Charts//");
		currentFolder.mkdirs();
		
		SamParser bp = null;
		// COV_RATE=windowLength/5;
		System.out.println("WIND_LENGTH=" + windowLength);
		System.out.println("COV_RATE=" + COV_RATE);
		System.out.println("MIN_WIND_LENGTH=" + MIN_WIND_LENGTH);
		System.out.println("MODE SMOOTHER WINDOW=" + k);
		System.out.println("NUMBER_OF_FIT_BINS=" + BIN_FACTOR*10);
		System.out.println("USED_DATA=" + USED_DATA * 100 + "%");
		
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
		System.out.println("TOTAL TIME : ["+currentProjectName+"] :"+ totalTime / 1000);
	}
	
	
	public static void runMultipleWindows(String[] args){
		String[] newArgs=new String[args.length+2];
		long startTimeGeneral = System.currentTimeMillis();
		String[] winLengths = {"500","750","1000","2000","3000"};//"50","100","250","400","500","750","1000","1500","2000", "3000","5000", "6500", "9000", "10000", "15000","20000", "30000", "40000", "50000","75000"};		
		for (int wlInd = 0; wlInd < winLengths.length; wlInd++) {
			try {
				
				for (int i=0;i<args.length;i++){
					newArgs[i]=args[i];
				}
				newArgs[args.length]= "-w";
				newArgs[args.length+1]= winLengths[wlInd];
				
				System.out.println( " ARGS for wl="+winLengths[wlInd]);
				for (int i=0;i<args.length;i++){
					System.out.print(newArgs[i] + " ");
				}
				System.out.println( " ");
				
				runPloest(newArgs);
				

			} catch (Exception e) {
				System.err.println(args[1] + args[3]+" could not be run");
				e.printStackTrace();
			}

		}
		long endTimeGeneral = System.currentTimeMillis();
		long totalTimeGeneral = endTimeGeneral - startTimeGeneral;
		System.out.println("TOTAL TIME GENERAL ("+winLengths.length+" run/s): " + totalTimeGeneral / 1000);
	}
	
	
	
	
	public static void printHelp() {
		System.out.println("\n +++++++++++++++++++++++++++++++++++++++++++++++++++++++");
		System.out.println("\nPedca -help:");
		System.out.println(
				"\nUSAGE:    java -jar Pedca.jar  -p <project name> -i < input sam/bam File> -o <output Folder> <<OPTIONAL PARAMETERS>> \n");
		System.out.println("REQUIRED PARAMETERS:");
		System.out.println("-p (project name)          - (String)  Prefix used to generate the results file.");
		System.out.println(
				"-i (input file)            - (Pathway to .bam/.sam file)   Pathway (absolute if necessary) to the input file containing the alignment file. Must be a .bam or .sam file");
		System.out.println(	"-o (output Folder)         - (String)   Pathway (absolute if necessary) to the auto-generated output folder that will contain the results('./' if the folder points to the current directory)");

		System.out.println("\nOPTIONAL PARAMETERS:");
		System.out.println(
				"-m (multi Run)            - (no parameters) Runs a preselected set of default window lengths {500,750,1000,2000,3000}");
		System.out.println(
				"-w (windows length)       - (int)      Length of the sliding window, to measure the coverage inside  contig. Default 500 bp");
		System.out.println(
				"-c (coverage rate)        - (int)    Rate factor for the coverage sampling in the Read count distribution. Default 100. The smaller it is, the less bins are sampled");

		System.out.println(
				"-k (mode smoother window) - (int)     Number of points over which the ploidy estimation is smoothed. The mode over k numbers of windows is used to average the values of the bin. Default=49 ");
		System.out.println(
				"-s (significant min)      - (double) Threshold to consider a cluster peak in the read count to be significant. Default 0.1");
		System.out.println(
				"-b (fitter bin factor)      - (double)  Affects the number of bins used to FIT the read count distribution. Default 2.5; Recommended between min=2.0 and max=4.0");
		System.out.println(
				"-d (coverage data to use) - (double) Fraction of coverage data that to be used in the read count distribution to infer the different ploidies and their ratio. Default 0.97");

		System.out.println(
				"-v (allele frequencies)   - (Pathway to .vcf file)  Pathway (absolute if necessary) to the file containing the variant calling. Must be a .vcf file");

		System.exit(1);

		// */
	}

}
