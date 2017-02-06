

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Scanner;


public class VCFManager {
	List<String> humanLines = new ArrayList<String>();
	ArrayList<ArrayList<Double>> vcfMatrix = new ArrayList<ArrayList<Double>>();
	String outputFileRoot=Ploest.outputFile+"\\"+Ploest.projectName+ "\\BaseCall";
	String endFix="";
	String outputHumanFile="" ;
	String outputMatrixFile="Matrix.vcf" ;
	
	
	public VCFManager(String inputVcf_file) throws FileNotFoundException, InterruptedException {//constructor from vcf file
		
		File pilonOutFolder =	 new File(outputFileRoot);
		pilonOutFolder.mkdirs();
		vcfExtractor(inputVcf_file);	
	}

	
	public void vcfExtractor(String inputFile) throws FileNotFoundException, InterruptedException {
		System.out.println("--------------------     vcfExtractor    ----------------------");
		String currentChromosome="";
		double chromNumber=0;
		double nbOfVars=0;
		//solve the paths
			
		outputMatrixFile = outputFileRoot  + "\\"+"MatFileOut.txt";
		System.out.println( "   outputMatrixFile" + outputMatrixFile );		
		// We are interested in the 6th column and the DEPTH column:
		// 6th="Count of As, Cs, Gs, Ts at locus"

		String line = "";
		Scanner sc = new Scanner(new File(inputFile));
		int ct = 0;
	String chrom="";
	String id="";
		int pos = 0;
		double  depth = 0;// depth
		double bcA = 0;
		double bcC = 0;
		double bcG = 0;
		double bcT = 0;

		List<String> formatSample;
		List<String> baseCalls;
		String nextFilter;
		boolean currentContigIsInBaseCallList=false;
		ArrayList<Double> currentMatrixLine = new ArrayList<Double>();
		String next;

		try {
			// skip the header			
			do{
				line = sc.nextLine();
				next=sc.next();
			}while (next.substring(0, 1).equals("#")); 
			chrom=next;
			if(!currentChromosome.equals(chrom)){
				currentChromosome=chrom;
				chromNumber++;
				if(NaivePloestPlotter.continousPloidyContigsNamesWithBasicUnit.contains(currentChromosome)){//check if new currentContig Is In BaseCall List
					currentContigIsInBaseCallList=true;
				}else currentContigIsInBaseCallList=false;
			}
			//get the values
			while (sc.hasNextLine()  /*&& ct < 340*/ ) {
				
				if(currentContigIsInBaseCallList){
					
					next=sc.next();	pos = Integer.parseInt(next);// get pos//sc.next();	//skip pos  //

					sc.next();// skip id// id=sc.next();
					sc.next();///skip ref//ref = sc.next();// get reference allele
					sc.next();//skip alt//alt = sc.next();// get alternative allele
					sc.next();// skip qual
					nextFilter = sc.next();// get the filter call. Ambiguous and Deletions are the interesting ones

					formatSample = Arrays.asList(sc.next().split(";"));// split all fields in the info line

					if (!formatSample.get(0).substring(0, 2).equals("DP")) {// Special vcf line ("STRUCTURAL VARIATION" )
						// DO NOTHING. If need to do something recover code from project SequenceSimulator. Class: VCFManager
					} else {// ...regular vcf line
						depth = (double)Integer.parseInt(formatSample.get(0).substring(3, formatSample.get(0).length()));
						baseCalls = Arrays
								.asList(formatSample.get(5).substring(3, formatSample.get(5).length()).split(","));
						if(depth!=0.0){
							bcA = Integer.parseInt(baseCalls.get(0))/depth;					
							bcC = Integer.parseInt(baseCalls.get(1))/depth; 
							bcG = Integer.parseInt(baseCalls.get(2))/depth; 	
							bcT = Integer.parseInt(baseCalls.get(3))/depth;	
						}
					
					}
					
					currentMatrixLine = new ArrayList<Double>();
					currentMatrixLine.add(depth);//6
					currentMatrixLine.add(bcA);//7
					currentMatrixLine.add(bcC);//8
					currentMatrixLine.add(bcG);//9
					currentMatrixLine.add(bcT);//10
				

					if (nextFilter.substring(0, 3).equals("Amb") || nextFilter.substring(0, 3).equals("Del")) {
						vcfMatrix.add(currentMatrixLine);
						nbOfVars++;
					}
					//previousMatrixLine = currentMatrixLine;
					
					
					if (sc.hasNextLine()) { // 'if' to avoid error at end of file
						line = sc.nextLine();
						chrom=sc.next();// contig name 
						if(!currentChromosome.equals(chrom)){//if change in chrom
							currentChromosome=chrom;
							chromNumber++;
							if(NaivePloestPlotter.continousPloidyContigsNamesWithBasicUnit.contains(currentChromosome)){//check if new currentContig Is In BaseCall List
								currentContigIsInBaseCallList=true;
							}else currentContigIsInBaseCallList=false;
						}
					}
					ct++;
				}else{
				
					line=sc.nextLine();
					if (sc.hasNextLine()) { // 'if' to avoid error at end of file					
						chrom=sc.next();// contig name 
						if(!currentChromosome.equals(chrom)){//if change in chrom
							currentChromosome=chrom;
							chromNumber++;
							if(NaivePloestPlotter.continousPloidyContigsNamesWithBasicUnit.contains(currentChromosome)){//check if new currentContig Is In BaseCall List
								currentContigIsInBaseCallList=true;
							}else currentContigIsInBaseCallList=false;
						}
					}else 
					ct++;
				}
				
			}
		    
			System.out.println("ready for printVCFmatrix:" );

			printVCFmatrix(nbOfVars);

			if (sc != null)
				sc.close();
		} catch (Exception e) {
			System.err.println("error at ct:" + ct +" chrom:"+chrom+ " id:"+id+" pos:" + pos+" current:  "+currentMatrixLine);
		}

	}

	public int bitSequence(String s) {
		char[] cseq = new char[s.length()];
		int chASCII;
		for (int i = 0; i < s.length(); i++) {
			chASCII = (int) Character.toUpperCase(s.charAt(i));
			// System.out.println("processing:"+chASCII+":"+s.charAt(i)+":");
			switch (chASCII) {
			case 65://'A'
				cseq[i] = '1';
				break;
			case 67://'C'
				cseq[i] = '2';
				break;
			case 71://'G'
				cseq[i] = '3';
				break;
			case 84://'T'
				cseq[i] = '4';
				break;
			case 46://'.'
				cseq[i] = '0';
				break;
			default:
				cseq[i] = '9';
				break;
			}
		}

		String strSeq = new String(cseq);

		int number = strToInt(strSeq);
		return number;
	}

	public int codeFilter(String s) {
		List<String> filters = Arrays.asList(s.split(";"));
		char[] cfilt = new char[filters.size()];
		for (int i = 0; i < filters.size(); i++) {
			String filt = filters.get(i);
			switch (filt.toUpperCase()) {
			case "PASS":
				cfilt[i] = '9';
				break;
			case "AMB":
				cfilt[i] = '1';
				break;
			case "LOWCOV":
				cfilt[i] = '5';
				break;
			case "DEL":
				cfilt[i] = '3';
				break;
			default:
				cfilt[i] = '0';
				break;
			}
		}
		return Integer.parseInt(new String(cfilt));
	}

	public static int strToInt(String str) {
		int i = 0;
		int num = 0;
		boolean isNeg = false;

		// Check for negative sign; if it's there, set the isNeg flag
		if (str.charAt(0) == '-') {
			isNeg = true;
			i = 1;
		}

		// Process each character of the string;
		while (i < str.length()) {
			num *= 10;
			num += str.charAt(i++) - '0'; // Minus the ASCII code of '0' to get
			// the value of the charAt(i++).
		}

		if (isNeg)
			num = -num;
		return num;
	}

	public void printVCFmatrix(double nbOfVars) {
		System.out.println("Destination path of outputMatrixFile: "+outputMatrixFile);
		PrintStream stdout = System.out;
		PrintStream myConsole = null;
		try {
			myConsole = new PrintStream(new File(outputMatrixFile));
			System.setOut(myConsole);
			ArrayList<Double> currentLine;
			int lineSize=0;
			
			if(vcfMatrix.size()>0) {
				lineSize=vcfMatrix.get(0).size();//should always be 5
				
				
				for (int i = 0; i < vcfMatrix.size(); i++) {
					currentLine=vcfMatrix.get(i);
					for(int j =0;j<lineSize-2;j++){
						System.out.print(currentLine.get(j)+" ");
					}
					System.out.println(currentLine.get(lineSize-1));
				}
			}
			System.out.println("# Total nb of variation positions =  "+vcfMatrix.size()+" over "+NaivePloestPlotter.LENGTH_OF_BASIC_UNIT_CONTIGS);
			myConsole.close();
		}  catch (Exception e) {
			System.err.println("Error trying to write outputMatrixFile ");
			e.printStackTrace();
		}  finally {
			if (myConsole != null) {	        	
				myConsole.close();
				System.setOut(stdout);  
			}
		}
		System.setOut(stdout);  
	}

	public void printHumanLines() throws FileNotFoundException {

		PrintStream stdout = System.out;
		PrintStream myConsole = null;
		
		try {
			myConsole = new PrintStream(new File(outputHumanFile));
			System.setOut(myConsole);
			for (int i = 0; i < humanLines.size(); i++) {
				System.out.println(humanLines.get(i));
			}
			myConsole.close();
		} catch (Exception e) {
			System.err.println("Error trying to write outputHumanFile");
			e.printStackTrace();
		}  finally {
			if (myConsole != null) {	        	
				myConsole.close();
				System.setOut(stdout);  
			}
		}
		System.setOut(stdout);  
		
	}
}
