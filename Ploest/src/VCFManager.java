

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Scanner;


public class VCFManager {
	List<String> humanLines = new ArrayList<String>();
	ArrayList<ArrayList<Integer>> vcfMatrix = new ArrayList<ArrayList<Integer>>();
	String outputFileRoot=Ploest.outputFile+"\\"+Ploest.projectName;
	String endFix="";
	String outputHumanFile="" ;
	String outputMatrixFile="Matrix.vcf" ;

	
	public VCFManager(String inputFile) throws FileNotFoundException, InterruptedException {//constructor from vcf file
		File pilonFolder =	 new File(outputFileRoot + "\\BaseCall");
		pilonFolder.mkdirs();
		vcfExtractor(inputFile);	
	}

	public VCFManager(String inputFile, boolean isVcfMatrix ) throws FileNotFoundException, InterruptedException {//constructor from matrix file
		if(isVcfMatrix){
			//solve the paths
			int lastIndex = inputFile.lastIndexOf("/");
			if (lastIndex ==-1) {
				lastIndex = inputFile.lastIndexOf("\\");
				if (lastIndex ==-1) {
					lastIndex = 0;
				} 

			}
			outputFileRoot = inputFile.substring(0, lastIndex);
			endFix=inputFile.substring(lastIndex+1, inputFile.length()-4);		
			outputHumanFile = outputFileRoot + "\\"+ endFix + "HumFileOut.txt";
			outputMatrixFile = outputFileRoot  + "\\"+ endFix + "MatFileOut.txt";
			
			//extract the file to vcfMatrix
			String line = "";
			Scanner MatSc = new Scanner(new File(inputFile));
			MatSc.useDelimiter("[ ,\r\n\t]");int ppp=0;
			try{
				int nextInt;
				while (MatSc.hasNextLine()) { ;	
					if (ppp<200000){
					ArrayList<Integer> intArrayLine=new ArrayList<Integer>();
					line = MatSc.nextLine();   // read full line
					Scanner lineScan = new Scanner(line);
					lineScan.useDelimiter("[ ,\r\n\t]");;
					//System.out.println("line="+line);
					for (int i=0;i<=15;i++){
											
						nextInt=lineScan.nextInt();
						//System.out.print("ppp="+ppp+" i:"+i+" nextInt:"+nextInt+" lineScan.next:"+lineScan.next()+"-");	
						//System.out.print(" nextInt:"+nextInt);
						intArrayLine.add(nextInt);
					}
					vcfMatrix.add(intArrayLine);
					ppp++;
					}else {
						break;
					}
									
					//vcfMatrix.add(new ArrayList<Integer>());
				}				
			}catch (Exception e) {
				System.out.println("!*$! Error constructing VCFManager from Matrix file at ppp="+ppp);
			}
			printVCFmatrix();
			
		}
	}
	
	public void vcfExtractor(String inputFile) throws FileNotFoundException, InterruptedException {

		//solve the paths
		int lastIndex = inputFile.lastIndexOf("/");
		if (lastIndex ==-1) {
			lastIndex = inputFile.lastIndexOf("\\");
			if (lastIndex ==-1) {
				lastIndex = 0;
			} 

		}
		outputFileRoot = inputFile.substring(0, lastIndex);
		endFix=inputFile.substring(lastIndex+1, inputFile.length()-4);		
		outputHumanFile = outputFileRoot + "\\"+ endFix + "HumFileOut.txt";
		outputMatrixFile = outputFileRoot  + "\\"+ endFix + "MatFileOut.txt";
		//System.out.println("outputHumanFile" +outputHumanFile + "   outputMatrixFile" + outputMatrixFile );		
		// We are interested in the 6th and 7th column:
		// 6th="Count of As, Cs, Gs, Ts at locus"
		// 7th="="Percentage of As, Cs, Gs, Ts weighted by Q & MQ at locus"
		String line = "";
		Scanner sc = new Scanner(new File(inputFile));
		int ct = 0;
		int chrom=0;
		int pos = 0;
		int depth = 0;// depth
		int bcA = 0;
		int bcC = 0;
		int bcG = 0;
		int bcT = 0;
		int qpA = 0;
		int qpC = 0;
		int qpG = 0;
		int qpT = 0;
		int sv;// structural variance (boolean)
		int svLength;// length of the sv
		String ref = "";// reference allele
		String alt;// alternative allele
		//System.out.println("POS" + "\t" + "REF" + "\t" + "ALT" + "\t" + "FILTER" + "\t" + "INFO");
		List<String> formatSample;
		List<String> baseCalls;
		List<String> qPercentages;
		String nextFilter;
		String currentHumanLine = "";
		String previousHumanLine = "";
		ArrayList<Integer> currentMatrixLine;
		ArrayList<Integer> previousMatrixLine = new ArrayList<Integer>(16);
		String next;
		try {
			// skip the header			
			do{
				line = sc.nextLine();
				next=sc.next();
			}while (next.substring(0, 1).equals("#")); 

			chrom=Integer.parseInt(next.substring(5,11));

			//get the values
			while (sc.hasNextLine() /* && ct < 34000 */) {
				next=sc.next();				
				pos = Integer.parseInt(next);// get pos
				sc.next();// skip id
				ref = sc.next();// get reference allele
				alt = sc.next();// get alternative allele
				sc.next();// skip qual
				nextFilter = sc.next();// get the filter call. Ambiguous and
				// Deletions are the interesting ones

				formatSample = Arrays.asList(sc.next().split(";"));// split all
				// fields in
				// the info
				// line

				if (!formatSample.get(0).substring(0, 2).equals("DP")) {// Special vcf line ("STRUCTURAL VARIATION" )

					depth = 0;
					bcA = 0;
					bcC = 0;
					bcG = 0;
					bcT = 0;
					qpA = 0;
					qpC = 0;
					qpG = 0;
					qpT = 0;
					sv = 1;
					svLength = Integer.parseInt(formatSample.get(1).substring(6, formatSample.get(1).length())) ;
					//System.out.println(" formatSample.get(2):"+formatSample+ " LENGTH:" + svLength);
					currentHumanLine = "chrom:"+chrom+"pos" + pos + " " + "ref:" + ref + " " + "alt:" + alt + " " + nextFilter + " "
							+ formatSample.get(0) + " LENGTH:" + svLength;
				} else {// ...regular vcf line
					depth = Integer.parseInt(formatSample.get(0).substring(3, formatSample.get(0).length()));
					baseCalls = Arrays
							.asList(formatSample.get(5).substring(3, formatSample.get(5).length()).split(","));

					bcA = Integer.parseInt(baseCalls.get(0));
					bcC = Integer.parseInt(baseCalls.get(1));
					bcG = Integer.parseInt(baseCalls.get(2));
					bcT = Integer.parseInt(baseCalls.get(3));

					qPercentages = Arrays
							.asList(formatSample.get(6).substring(3, formatSample.get(6).length()).split(","));

					qpA = Integer.parseInt(qPercentages.get(0));
					qpC = Integer.parseInt(qPercentages.get(1));
					qpG = Integer.parseInt(qPercentages.get(2));
					qpT = Integer.parseInt(qPercentages.get(3));
					sv = 0;
					svLength = 0;
					currentHumanLine = "chrom:"+chrom+" pos" + pos + " ref:" + ref + " alt:" + alt + " " + nextFilter + " "
							+ formatSample.get(0).substring(3, formatSample.get(0).length()) + " "
							+ formatSample.get(5).substring(3, formatSample.get(5).length()) + " "
							+ formatSample.get(6).substring(3, formatSample.get(6).length());
				}
				currentMatrixLine = new ArrayList<Integer>();
				currentMatrixLine.add(chrom);//1
				currentMatrixLine.add(pos);//2
				currentMatrixLine.add(bitSequence(ref));//3
				currentMatrixLine.add(bitSequence(alt));//4
				currentMatrixLine.add(codeFilter(nextFilter));//5
				currentMatrixLine.add(depth);//6
				currentMatrixLine.add(bcA);//7
				currentMatrixLine.add(bcC);//8
				currentMatrixLine.add(bcG);//9
				currentMatrixLine.add(bcT);//10
				currentMatrixLine.add(qpA);//11
				currentMatrixLine.add(qpC);//12
				currentMatrixLine.add(qpG);//13
				currentMatrixLine.add(qpT);//14

				currentMatrixLine.add(sv);//15
				currentMatrixLine.add(svLength);//16
				/*
				// if deletion or ambiguous only
				if (nextFilter.substring(0, 3).equals("Amb") || nextFilter.substring(0, 3).equals("Del")) {

					codeFilter(nextFilter);
					if (nextFilter.substring(0, 3).equals("Del")
							&& !previousMatrixLine.equals(vcfMatrix.get(vcfMatrix.size() - 1))) {
						// store the previous line avoiding repeats
						vcfMatrix.add(previousMatrixLine);
						humanLines.add(previousHumanLine);
					}
					humanLines.add(currentHumanLine);
					vcfMatrix.add(currentMatrixLine);
				}*/
				humanLines.add(currentHumanLine);
				vcfMatrix.add(currentMatrixLine);

				previousMatrixLine = currentMatrixLine;
				previousHumanLine = currentHumanLine;



				line = sc.nextLine();
				if (sc.hasNextLine()) { // 'if' to avoid error at end of file
					//String eraseme=sc.next().substring(5,11);
					//System.out.println("sc.next():"+eraseme+":");
					chrom=Integer.parseInt(sc.next().substring(5,11));// contig name 
					//System.out.println("chrom:"+chrom);
				}
				ct++;
			}

			printHumanLines();
			printVCFmatrix();


			if (sc != null)
				sc.close();
		} catch (Exception e) {
			System.out.println("error at ct:" + ct + " pos:" + pos+" current:  "+currentHumanLine);
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

	public void printVCFmatrix() {
		System.out.println("Destination path of outputMatrixFile: "+outputMatrixFile);
		PrintStream stdout = System.out;
		PrintStream myConsole = null;
		try {
			myConsole = new PrintStream(new File(outputMatrixFile));
			System.setOut(myConsole);
			for (int i = 0; i < vcfMatrix.size(); i++) {
				System.out.println(vcfMatrix.get(i));
			}
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
