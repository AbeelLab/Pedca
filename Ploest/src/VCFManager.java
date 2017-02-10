
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
	String outputFileRoot = Ploest.outputFile + "\\" + Ploest.projectName + "\\BaseCall";
	String endFix = "";
	String outputMatrixFile1 = "Matrix1stCluster.vcf";
	String outputMatrixFile2 = "Matrix2ndCluster.vcf";
	static int depthThreshold=10;//coverage threshold to consider a variation valid
	
	// constructor from vcf file
	public VCFManager(String inputVcf_file) throws FileNotFoundException, InterruptedException {
		File pilonOutFolder = new File(outputFileRoot);
		pilonOutFolder.mkdirs();
		vcfExtractor(inputVcf_file);
	}

	public void vcfExtractor(String inputFile) throws FileNotFoundException, InterruptedException {
		System.out.println("--------------------     vcfExtractor    ----------------------");
		String currentChromosome = "";
		double chromNumber = 0;
		int nbOfVarsCluster1 = 0;
		int nbOfVarsCluster2 = 0;
		
		// solve the paths
		outputMatrixFile1 = outputFileRoot + "\\" + outputMatrixFile1;
		System.out.println("   outputMatrixFile1" + outputMatrixFile1);
		outputMatrixFile2 = outputFileRoot + "\\" + outputMatrixFile2;
		System.out.println("   outputMatrixFile2" + outputMatrixFile2);
		// We are interested in the 6th column and the DEPTH column:
		// 6th="Count of As, Cs, Gs, Ts at locus"

		String line = "";
		Scanner sc = new Scanner(new File(inputFile));
		int ct = 0;
		String chrom = "";
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
			chrom = next;
			
			//treat first vcf line
			if (!currentChromosome.equals(chrom)) {
				currentChromosome = chrom;
				chromNumber++;
				
				// check if first Contig is in BaseCall List: cluster1
				if (NaivePloestPlotter.continousPloidyContigsNamesCluster1.contains(currentChromosome)) {
					currentContigIsInBaseCall_CLUSTER_1_List = true;
				} else
					currentContigIsInBaseCall_CLUSTER_1_List = false;
				// check if first Contig is in BaseCall List: cluster2
				if (NaivePloestPlotter.continousPloidyContigsNamesCluster2.contains(currentChromosome)) {
					currentContigIsInBaseCall_CLUSTER_2_List = true;
				} else
					currentContigIsInBaseCall_CLUSTER_2_List = false;
			}
			
			// get the values
			while (sc.hasNextLine() ) {

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

					if ((nextFilter.substring(0, 3).equals("Amb") || nextFilter.substring(0, 3).equals("Del"))&&(depth > 0.0)) {
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
								// check if new currentContig is in BaseCall List: cluster1
								if (NaivePloestPlotter.continousPloidyContigsNamesCluster1.contains(currentChromosome)) {
									currentContigIsInBaseCall_CLUSTER_1_List = true;
								} else
									currentContigIsInBaseCall_CLUSTER_1_List = false;
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
						if (depth != 0.0) {
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

					if ((nextFilter.substring(0, 3).equals("Amb") || nextFilter.substring(0, 3).equals("Del"))&&(depth > 0.0)) {
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
								// check if new currentContig is in BaseCall List: cluster2
								if (NaivePloestPlotter.continousPloidyContigsNamesCluster2.contains(currentChromosome)) {
									currentContigIsInBaseCall_CLUSTER_2_List = true;
								} else
									currentContigIsInBaseCall_CLUSTER_2_List = false;
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

							// check if new currentContig is in BaseCall List: cluster1
							if (NaivePloestPlotter.continousPloidyContigsNamesCluster1.contains(currentChromosome)) {
								currentContigIsInBaseCall_CLUSTER_1_List = true;
							} else{
								currentContigIsInBaseCall_CLUSTER_1_List = false;
								if (NaivePloestPlotter.continousPloidyContigsNamesCluster2.contains(currentChromosome)) {
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
			System.err.println("error at ct:" + ct + " chrom:" + chrom + " id:" + id + " pos:" + pos + " current:  "
					+ currentMatrixLine);
		}

	}

	

	public void printVCFmatrix1(double nbOfVars) {
		System.out.println("Destination path of outputMatrixFile 1: " + outputMatrixFile1);
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
					+ NaivePloestPlotter.LENGTH_OF_CLUSTER_ONE_CONTIGS + "bp which is "
					+ (100 * (double) vcfMatrix1.size() / NaivePloestPlotter.LENGTH_OF_CLUSTER_ONE_CONTIGS)
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
		System.out.println("End printVCFmatrix1 ");
		BaseCallBarchar(vcfMatrix1,1);
		
	}

	private void BaseCallBarchar(ArrayList<ArrayList<Double>> vcfMatrix,int cluster) {
		System.out.println("BaseCallBarchar called for  cluster "+cluster);
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
			
			double[]baseCalls=new double[baseCallsList.size()];
			for (int j = 0; j < baseCallsList.size(); j++) {
				baseCalls[j]=baseCallsList.get(j);
			}
			System.out.println("BaseCallBarchar call 2");
			BarChart barchart = new BarChart(baseCalls,cluster);
			System.out.println("BaseCallBarchar END");
		}	
	}

	
	
	public void printVCFmatrix2(double nbOfVars2) {
		System.out.println("Destination path of outputMatrixFile2: " + outputMatrixFile2);
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
					+ NaivePloestPlotter.LENGTH_OF_CLUSTER_TWO_CONTIGS + "bp which is "
					+ (100 * (double) vcfMatrix2.size() / NaivePloestPlotter.LENGTH_OF_CLUSTER_TWO_CONTIGS)
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
		BaseCallBarchar(vcfMatrix2,2);
	}
}
