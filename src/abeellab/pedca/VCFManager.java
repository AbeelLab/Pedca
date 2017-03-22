package abeellab.pedca;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Scanner;

public class VCFManager {
	
	List<String> humanLines = new ArrayList<String>();
	ArrayList<ArrayList<Double>> vcfMatrix1 = new ArrayList<ArrayList<Double>>();
	ArrayList<ArrayList<Double>> vcfMatrix2 = new ArrayList<ArrayList<Double>>();
	String outputFileRoot = Ploest.outputFile +  "\\BaseCall";//= Ploest.outputFile + "\\" + Ploest.projectName + "\\BaseCall";
	String endFix = "";
	String outputMatrixFile1 = "Matrix1stCluster.vcf";
	String outputMatrixFile2 = "Matrix2ndCluster.vcf";
	final static int depthThreshold=10;//coverage threshold to consider a variation valid
	NaivePloestPlotter ploestPlotter;
	
	// constructor from vcf file
	public VCFManager(String inputVcf_file,NaivePloestPlotter n) throws FileNotFoundException, InterruptedException {
		System.out.println("cluster 1 size:"+n.LENGTH_OF_CLUSTER_ONE_CONTIGS+" cluster 2 size:"+n.LENGTH_OF_CLUSTER_TWO_CONTIGS);
		ploestPlotter=n;
		File pilonOutFolder = new File(outputFileRoot);
		pilonOutFolder.mkdirs();
		vcfExtractor(inputVcf_file);
	}

	public void vcfExtractor(String inputFile) throws FileNotFoundException, InterruptedException {
		System.out.println("--------------------     vcfExtractor    ----------------------"+inputFile);
		System.out.println("   ContigsNamesCluster1" + ploestPlotter.continousPloidyContigsNamesCluster1.size()+" NaivePloestPlotter.continousPloidyContigsNamesCluster2:"+ploestPlotter.continousPloidyContigsNamesCluster2.size());
		String currentChromosome = "";
		double chromNumber = 0;
		int nbOfVarsCluster1 = 0;
		int nbOfVarsCluster2 = 0;
		
		// solve the paths
		outputMatrixFile1 = outputFileRoot + "\\" + outputMatrixFile1;
		System.out.println("   outputMatrixFile1 " + outputMatrixFile1);
		outputMatrixFile2 = outputFileRoot + "\\" + outputMatrixFile2;
		System.out.println("   outputMatrixFile2 " + outputMatrixFile2);
		// We are interested in the 6th column and the DEPTH column:
		// 6th="Count of As, Cs, Gs, Ts at locus"

		String line = "";
		Scanner sc = new Scanner(new File(inputFile));
		int ct = 0;
		String chrom = "*";
		String id = "";
		int pos = 0;
		double depth = 0;// depth
		double bcA = 0;
		double bcC = 0;
		double bcG = 0;
		double bcT = 0;

		List<String> formatSample;
		List<String> baseCalls;
		String nextFilter;
		boolean currentContigIsInBaseCall_CLUSTER_1_List = false;
		boolean currentContigIsInBaseCall_CLUSTER_2_List = false;
		ArrayList<Double> currentMatrixLine = new ArrayList<Double>();
		String next;

		try {
			// skip the header
			do {
				line = sc.nextLine();
				next = sc.next();

			} while (next.substring(0, 1).equals("#"));
			chrom = next;System.out.println("First chrom: " + chrom);
			
			//treat first vcf line
			if (!currentChromosome.equals(chrom)) {
				currentChromosome = chrom;
				chromNumber++;
				// check if first Contig is in BaseCall List: cluster1
				if (ploestPlotter.continousPloidyContigsNamesCluster1.contains(currentChromosome)) {
					currentContigIsInBaseCall_CLUSTER_1_List = true;
				} else
					currentContigIsInBaseCall_CLUSTER_1_List = false;
				// check if first Contig is in BaseCall List: cluster2
				if (ploestPlotter.continousPloidyContigsNamesCluster2.contains(currentChromosome)) {
					currentContigIsInBaseCall_CLUSTER_2_List = true;
				} else
					currentContigIsInBaseCall_CLUSTER_2_List = false;
			}
			
			// get the values
			while (sc.hasNextLine() ) {
				//if(pos%100000==1)System.out.println("chrom:"+chrom+" pos: " + pos);

		
				if (currentContigIsInBaseCall_CLUSTER_1_List) {// treates variations in contigs from cluster 1
					next = sc.next();
					
					pos = Integer.parseInt(next);// get pos//sc.next(); //skip pos //
					sc.next();// skip id// id=sc.next();
					sc.next();/// skip ref//ref = sc.next();// get reference allele
					sc.next();// skip alt//alt = sc.next();// get alternative allele
					sc.next();// skip qual
					nextFilter = sc.next();// get the filter call. Ambiguous and
											// Deletions are the interesting
											// ones

					formatSample = Arrays.asList(sc.next().split(";"));// splits all fields in the info line

					if (!formatSample.get(0).substring(0, 2).equals("DP")) {// Special vcf line ("STRUCTURAL VARIATION" )
						// DO NOTHING. If need to do something recover code from
						// project SequenceSimulator. Class: VCFManager
					} else {// ...regular vcf line
						depth = (double) Integer.parseInt(formatSample.get(0).substring(3, formatSample.get(0).length()));
						baseCalls = Arrays.asList(formatSample.get(5).substring(3, formatSample.get(5).length()).split(","));
						if (depth > depthThreshold) {
							bcA = Integer.parseInt(baseCalls.get(0)) / depth;
							bcC = Integer.parseInt(baseCalls.get(1)) / depth;
							bcG = Integer.parseInt(baseCalls.get(2)) / depth;
							bcT = Integer.parseInt(baseCalls.get(3)) / depth;
						}
					}

					currentMatrixLine = new ArrayList<Double>();
					currentMatrixLine.add(depth);// 6
					currentMatrixLine.add(bcA);// 7
					currentMatrixLine.add(bcC);// 8
					currentMatrixLine.add(bcG);// 9
					currentMatrixLine.add(bcT);// 10

					if ((nextFilter.substring(0, 3).equals("Amb") || nextFilter.substring(0, 3).equals("Del"))&&(depth > depthThreshold)) {
						vcfMatrix1.add(currentMatrixLine);
						nbOfVarsCluster1++;
					}
					//check next chromosome name
					if (sc.hasNextLine()) { // 'if' to avoid error at end of file
						line = sc.nextLine();
						if (sc.hasNextLine()) {// avoid last line empty (sic!)
							chrom = sc.next();// contig name
							if (!currentChromosome.equals(chrom)) {// if change in chrom
								currentChromosome = chrom;
								chromNumber++;
								//printout current chromosome being processed
								//System.out.println("Allele frequencies analysis. Processing chrom:"+chrom+" in vcf file....");
								// check if new currentContig is in BaseCall List:
								if (ploestPlotter.continousPloidyContigsNamesCluster1.contains(currentChromosome)) {
									currentContigIsInBaseCall_CLUSTER_1_List = true;
								} else{
									currentContigIsInBaseCall_CLUSTER_1_List = false;
									if (ploestPlotter.continousPloidyContigsNamesCluster2.contains(currentChromosome)) {
										currentContigIsInBaseCall_CLUSTER_2_List = true;
									} else{
										currentContigIsInBaseCall_CLUSTER_2_List = false;
									}
								}
							}
						}

					}

					ct++;
				} else if (currentContigIsInBaseCall_CLUSTER_2_List) {// treates variations in contigs from cluster 2
					next = sc.next();
					pos = Integer.parseInt(next);// get pos//sc.next(); //skip pos //

					sc.next();// skip id// id=sc.next();
					sc.next();/// skip ref//ref = sc.next();// get reference allele
					sc.next();// skip alt//alt = sc.next();// get alternative allele
					sc.next();// skip qual
					nextFilter = sc.next();// get the filter call. Ambiguous and Deletions are the interesting ones

					formatSample = Arrays.asList(sc.next().split(";"));// splits all fields in the info line

					if (!formatSample.get(0).substring(0, 2).equals("DP")) {// Special vcf line ("STRUCTURAL VARIATION" )
						// DO NOTHING. If need to do something recover code from
						// project SequenceSimulator. Class: VCFManager
					} else {// ...regular vcf line
						depth = (double) Integer
								.parseInt(formatSample.get(0).substring(3, formatSample.get(0).length()));
						baseCalls = Arrays.asList(formatSample.get(5).substring(3, formatSample.get(5).length()).split(","));
						if (depth > depthThreshold) {
							bcA = Integer.parseInt(baseCalls.get(0)) / depth;
							bcC = Integer.parseInt(baseCalls.get(1)) / depth;
							bcG = Integer.parseInt(baseCalls.get(2)) / depth;
							bcT = Integer.parseInt(baseCalls.get(3)) / depth;
						}
					}

					currentMatrixLine = new ArrayList<Double>();
					currentMatrixLine.add(depth);// 6
					currentMatrixLine.add(bcA);// 7
					currentMatrixLine.add(bcC);// 8
					currentMatrixLine.add(bcG);// 9
					currentMatrixLine.add(bcT);// 10

					if ((nextFilter.substring(0, 3).equals("Amb") || nextFilter.substring(0, 3).equals("Del"))&&(depth > depthThreshold)) {
						vcfMatrix2.add(currentMatrixLine);
						nbOfVarsCluster2++;
					}
					
					//check next chromosome name
					if (sc.hasNextLine()) { // 'if' to avoid error at end of file
						line = sc.nextLine();
						if (sc.hasNextLine()) {// avoid last line empty (sic!)
							chrom = sc.next();// contig name
							if (!currentChromosome.equals(chrom)) {// if change in chrom
								currentChromosome = chrom;
								chromNumber++;	
								//printout current chromosome being processed
								//System.out.println("Allele frequencies analysis. Processing chrom:"+chrom+" in vcf file....");
								// check if new currentContig is in BaseCall List: 
								if (ploestPlotter.continousPloidyContigsNamesCluster1.contains(currentChromosome)) {
									currentContigIsInBaseCall_CLUSTER_1_List = true;
								} else{
									currentContigIsInBaseCall_CLUSTER_1_List = false;
									if (ploestPlotter.continousPloidyContigsNamesCluster2.contains(currentChromosome)) {
										currentContigIsInBaseCall_CLUSTER_2_List = true;
									} else{
										currentContigIsInBaseCall_CLUSTER_2_List = false;
									}
								}
							}
						}
					}
					ct++;
				} else { //we don't care about variations in these cotigs, so ignore but check when a new contig appears

					line = sc.nextLine();
					if (sc.hasNextLine()) { // 'if' to avoid error at end of file
						chrom = sc.next();// contig name
						if (!currentChromosome.equals(chrom)) {// if change in chrom/contig
							currentChromosome = chrom;
							chromNumber++;
							//printout current chromosome being processed							
							//System.out.println("Allele frequencies analysis. Processing chrom:"+chrom+" in vcf file....");
							// check if new currentContig is in BaseCall List: cluster1
							if (ploestPlotter.continousPloidyContigsNamesCluster1.contains(currentChromosome)) {
								currentContigIsInBaseCall_CLUSTER_1_List = true;
							} else{
								currentContigIsInBaseCall_CLUSTER_1_List = false;
								if (ploestPlotter.continousPloidyContigsNamesCluster2.contains(currentChromosome)) {
									currentContigIsInBaseCall_CLUSTER_2_List = true;
								} else{
									currentContigIsInBaseCall_CLUSTER_2_List = false;
								}
							}				
						}
					} else
						ct++;
				}

			}
			printVCFmatrix1(nbOfVarsCluster1);
			printVCFmatrix2(nbOfVarsCluster2);
			if (sc != null)
				sc.close();
		} catch (Exception e) {
			e.printStackTrace();
			System.err.println("error at ct:" + ct + " chrom:" + chrom + " id:" + id + " pos:" + pos + " current:  "
					+ currentMatrixLine);
		}

	}

	

	public void printVCFmatrix1(double nbOfVars) {
		System.out.println("Destination path of outputMatrixFile 1: " + outputMatrixFile1+" size:"+vcfMatrix1.size());
		PrintStream stdout = System.out;
		PrintStream myConsole = null;
		try {
			myConsole = new PrintStream(new File(outputMatrixFile1));
			System.setOut(myConsole);
			ArrayList<Double> currentLine;
			int lineSize = 0;

			if (vcfMatrix1.size() > 0) {
				lineSize = vcfMatrix1.get(0).size();// should always be 5

				for (int i = 0; i < vcfMatrix1.size(); i++) {
					currentLine = vcfMatrix1.get(i);
					//print to text file
					for (int j = 0; j < lineSize - 1; j++) {
						System.out.print(currentLine.get(j) + " ");
					}
					System.out.println(currentLine.get(lineSize - 1));
			
				}
			}
			System.out.println("# Total nb of variation positions =  " + vcfMatrix1.size() + " over "
					+ ploestPlotter.LENGTH_OF_CLUSTER_ONE_CONTIGS + "bp which is "
					+ (100 * (double) vcfMatrix1.size() / ploestPlotter.LENGTH_OF_CLUSTER_ONE_CONTIGS)
					+ " % of the length of all contigs with that ploidy");
			myConsole.close();
		} catch (Exception e) {
			System.err.println("Error trying to write outputMatrixFile ");
			e.printStackTrace();
		} finally {
			if (myConsole != null) {
				myConsole.close();
				System.setOut(stdout);
			}
		}
		System.setOut(stdout);
		
		
	}

	private double[] BaseCallBarchar(ArrayList<ArrayList<Double>> vcfMatrix,int cluster) {

		double[]baseCalls=new double[0];
		if (vcfMatrix.size() > 0) {
			int lineSize = vcfMatrix.get(0).size();// should always be 5
			ArrayList<Double> currentLine;
			double currentPercent;
			ArrayList<Double> baseCallsList=new ArrayList<Double>();
			
			for (int i = 0; i < vcfMatrix.size(); i++) {
				currentLine=vcfMatrix.get(i);
				//print to text file
				for (int j = 1; j < lineSize; j++) {
					currentPercent=currentLine.get(j);
					if(currentPercent>0.1 && currentPercent<0.9){
						baseCallsList.add(currentPercent);
					}
				}
			}
			
			baseCalls=new double[baseCallsList.size()];
			for (int j = 0; j < baseCallsList.size(); j++) {
			
				baseCalls[j]=baseCallsList.get(j);
			}
			
		}
		return baseCalls;
			
	}

	
	
	public void printVCFmatrix2(double nbOfVars2) {
		System.out.println("Destination path of outputMatrixFile2: " + outputMatrixFile2+" size:"+vcfMatrix2.size());
		PrintStream stdout = System.out;
		PrintStream myConsole = null;
		try {
			myConsole = new PrintStream(new File(outputMatrixFile2));
			System.setOut(myConsole);
			ArrayList<Double> currentLine;
			int lineSize = 0;

			if (vcfMatrix2.size() > 0) {
				lineSize = vcfMatrix2.get(0).size();// should always be 5

				for (int i = 0; i < vcfMatrix2.size(); i++) {
					currentLine = vcfMatrix2.get(i);
					for (int j = 0; j < lineSize - 1; j++) {
						System.out.print(currentLine.get(j) + " ");
					}
					System.out.println(currentLine.get(lineSize - 1));
				}
			}
			System.out.println("# Total nb of variation positions =  " + vcfMatrix2.size() + " over "
					+ ploestPlotter.LENGTH_OF_CLUSTER_TWO_CONTIGS + "bp which is "
					+ (100 * (double) vcfMatrix2.size() / ploestPlotter.LENGTH_OF_CLUSTER_TWO_CONTIGS)
					+ " % of the length of all contigs with that ploidy");
			myConsole.close();
		} catch (Exception e) {
			System.err.println("Error trying to write outputMatrixFile2 ");
			e.printStackTrace();
		} finally {
			if (myConsole != null) {
				myConsole.close();
				System.setOut(stdout);
			}
		}
		System.setOut(stdout);
		
		//get the plot values and transform it into histogram values
		//for cluster1
		int[] YbaseCalls1=new int[20];
		double[] baseCalls1 =BaseCallBarchar(vcfMatrix1,1);
		int currentBin=0;
		for (int j = 0; j < baseCalls1.length; j++) {
			currentBin=(int) (baseCalls1[j]/0.05);
			YbaseCalls1[currentBin]++;
		}
		//and for cluster2
		double[] baseCalls2 =BaseCallBarchar(vcfMatrix2,2);
		currentBin=0;
		int[] YbaseCalls2=new int[20];
		for (int j = 0; j < baseCalls2.length; j++) {
			currentBin=(int) (baseCalls2[j]/0.05);
			YbaseCalls2[currentBin]++;
		}
		
		
		
		//find the max value of y axis in YbaseCalls1 and YbaseCalls2 (for comparative reasons we wnt the same y axe range)
		int maxYBaseCall=0;
		for(int i=0;i<YbaseCalls1.length;i++){
			if (YbaseCalls1[i]>maxYBaseCall)maxYBaseCall=YbaseCalls1[i];
		}
		for(int i=0;i<YbaseCalls2.length;i++){
			if (YbaseCalls2[i]>maxYBaseCall)maxYBaseCall=YbaseCalls2[i];
		}
		
		new BarChart(baseCalls1,1,maxYBaseCall);
		new BarChart(baseCalls2,2,maxYBaseCall);
		
	}
}
