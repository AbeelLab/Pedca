package abeellab.pedca;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.io.UnsupportedEncodingException;

import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMFileWriter;
import net.sf.samtools.SAMFileWriterFactory;
import net.sf.samtools.SAMLineParser;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;
import net.sf.samtools.SAMFileReader.ValidationStringency;
import java.util.*;

import org.jfree.data.xy.XYDataset;

import jMEF.PVector;

public class SamParser {
	static String debuggingTarget="cerevisiaeS288cchromosomeI";
	
	
	int nbSeq=0;// nb of sequences in the FileHeader
	Map<String, ContigData> contigsList;// Map of ContigDatas(value) and their  name (key)
	static int[] readCounts;
	
	static PVector[] fitPoints;	//remains for ploestpotter but it is not longer used
	static int[] intFitPoints;//all the points of all windows positions (coverages) in all contigs. Used to fit the read counts distribution chart

	static boolean RUN_SECOND_ROUND=false;
	static String stringSecondRound="";
	static int windowLength ;
	List<String> contArrList;
	static int readsDistributionMaxCoverage;//max coverage found in all contigs (to be used in x axis reads counts)
    static float readDistributionMaxY=0;//max normalized value in the y axis (to be used in y axis reads counts)
	PedcaPlotter plotter;//ploidy estimation  plott and pdf gaussian fit data
	NaivePedcaPlotter myploter;
	static BarChart barchart;
	
	static int totalDataPoints=0;//total amount of datapoints (all contigs considered)


	public SamParser(String inputFile, String outputfile)
			throws FileNotFoundException, UnsupportedEncodingException {
		
		
		//reset all static variables
		readCounts=null;
		intFitPoints=null;
		RUN_SECOND_ROUND=false;
		stringSecondRound="";
		windowLength =0;
		readsDistributionMaxCoverage=0;
		readDistributionMaxY=0;
		barchart=null;
		

		this.windowLength=Pedca.windowLength;
		
		File samFileInput=new File(inputFile);
		SAMFileReader inputSam = new SAMFileReader(samFileInput);
	
		nbSeq = inputSam.getFileHeader().getSequenceDictionary().size();// nb of sequences in the FileHeader
		contigsList = new HashMap<String, ContigData>(nbSeq);// Map of  ContigDatas(value) and their name (key)
		int[] contigLengths = new int[nbSeq];//length of all contigs
		int contigLength=0;
		String contigName="";
		// fill the contigsList
		for (int i = 0; i < nbSeq; i++) {
			contigLength=inputSam.getFileHeader().getSequenceDictionary().getSequence(i).getSequenceLength();
			contigName=inputSam.getFileHeader().getSequenceDictionary().getSequence(i).getSequenceName();
			//Printout contig names and lengths
			//System.out.println(contigName+","+contigLength/1000+" Kbp");
			contigsList.put(contigName,new ContigData(contigName,contigLength));
			contigLengths[i]=contigLength;
		}
		computeNmeasure(contigLengths,50);

		inputSam.setValidationStringency(ValidationStringency.SILENT);
		SAMRecordIterator iter = inputSam.iterator();

		int i = 1;
		System.out.println("Analyzing "+contigsList.size()+" contigs");
		String refName = "";// for debugging a bad line in the .sam file
		int alStart;

		while (iter.hasNext()) {//iterates the sam file

			SAMRecord rec = iter.next();
			refName = rec.getReferenceName();
			alStart = rec.getAlignmentStart();

			try {
				contigsList.get(refName).setPos(alStart);//storing the starting positions in the corresponing contig
			} catch (Exception e) {
			}
			i++;

		}
		iter.close();
		inputSam.close();

		windowSlideContigList();
		barchart = new BarChart(readCounts);
		readDistributionMaxY=barchart.maxY;
		
		myploter = new NaivePedcaPlotter(contigsList,readsDistributionMaxCoverage, barchart.normReadCounts);//plotter = new PloestPlotter(contigsList,maxWindows);
		
		if (RUN_SECOND_ROUND){//runs a second round with a new window length to solve smallests contigs
			System.out.println("-*-*-*-*-*-*-RUN_SECOND_ROUND*-*-*-*-*-*-*-*-*-*-*");
		
			Map<String, ContigData> newContigsList=new HashMap <String, ContigData> ();
			for (int c=0;c<myploter.unsolvedPloidyContigs.size();c++){
				newContigsList.put(myploter.unsolvedPloidyContigs.get(c).contigName, myploter.unsolvedPloidyContigs.get(c));
			}
			contigsList=newContigsList;
			windowLength=Pedca.windowLength;
			windowSlideContigList();
			//barchart = new BarChart(readCounts);//We don't need to plot this again
			readDistributionMaxY=barchart.maxY;
		    myploter.naivePloestPlotter2ndRound(contigsList,readsDistributionMaxCoverage, barchart.normReadCounts);
		}
		
		myploter.rt.writer.close();
		if(myploter.rt.writer2ndRun!=null)myploter.rt.writer2ndRun.close();
	}


	private void computeNmeasure(int[] contigLengths, int x){
		Arrays.sort(contigLengths);

		int sum=0;
		for (int i=contigLengths.length-1;i>=0;i--){
			sum+=contigLengths[i];
		}
		System.out.println("Total length of all contigs: "+sum+" bp");
		int target=sum*x/100;
		int cumulativeSum=0;
		for (int i=contigLengths.length-1;i>0;i--){
			cumulativeSum+=contigLengths[i];
			if (cumulativeSum>target){
				System.out.println("N"+x+" of input contigs:"+contigLengths[i]+" bp");
				break;
			}
		}
		
	}


	public SamParser(File fin, String outputfile)throws IOException {

		this.windowLength=Pedca.windowLength;
		List<String> listOfInputFiles=new ArrayList<String>();
		List<SAMFileReader> listOfSAMreaders=new ArrayList<SAMFileReader>();
		FileInputStream fis = new FileInputStream(fin);

		//Construct BufferedReader from InputStreamReader, and fill list of inputFiles
		BufferedReader br = new BufferedReader(new InputStreamReader(fis));	 
		String line = null;
		try {
			while ((line = br.readLine()) != null) {
				listOfInputFiles.add(line);
			}
		} catch (IOException e1) {
			e1.printStackTrace();
		}	 
		br.close();
		//determine nbOfSeq in all inputs
		SAMFileReader inputSam=null;
		for (int i=0;i<listOfInputFiles.size();i++){

			inputSam = new SAMFileReader(new File(listOfInputFiles.get(i)));
			listOfSAMreaders.add(inputSam);
			nbSeq = inputSam.getFileHeader().getSequenceDictionary().size();// nb of
			// sequences
			// in
			// the
			// FileHeader


		}
		contigsList = new HashMap<String, ContigData>(nbSeq*2);// Map of ContigDatas  name (key) and the atual ContigData (value)

		// fill the contigsList

		for (int sr=0;sr<listOfSAMreaders.size();sr++){//for all SAMreaders

			inputSam = listOfSAMreaders.get(sr);
			//fill contisList
			for (int i = 0; i < inputSam.getFileHeader().getSequenceDictionary().size(); i++) {
				contigsList.put(inputSam.getFileHeader().getSequenceDictionary().getSequence(i).getSequenceName(),
						new ContigData(inputSam.getFileHeader().getSequenceDictionary().getSequence(i).getSequenceName(),
								inputSam.getFileHeader().getSequenceDictionary().getSequence(i).getSequenceLength()));
			}

			inputSam.setValidationStringency(ValidationStringency.SILENT);
			SAMRecordIterator iter = inputSam.iterator();


			System.out.println("Analyzing "+contigsList.size()+" contigs for "+sr+ " libraries");
			String refName = "";// for debugging a bad line in the .sam file
			int alStart;

			while (iter.hasNext()) {//iterates the sam file

				SAMRecord rec = iter.next();
				refName = rec.getReferenceName();
				alStart = rec.getAlignmentStart();
				try {
					contigsList.get(refName).setPos(alStart);//storing the starting positions in the corresponing contig
				} catch (Exception e) {
				}


			}
			iter.close();
			inputSam.close();
		}
		windowSlideContigList();
		barchart = new BarChart(readCounts);
		plotter = new PedcaPlotter(contigsList,readsDistributionMaxCoverage);
	}


	public  int findMinimumContigLength(ArrayList<Integer> contLengt){ 
		  int minL=contLengt.get(0);
		  for (int i = 0;i< contLengt.size();i++){
			  if (contLengt.get(i)<minL)minL=contLengt.get(i);
		  }
		   
		  return minL; 
	} 
	
	
	 public double n98(ArrayList<Integer> lenghts, double n){ 
		  Collections.sort(lenghts); 
		  Collections.reverse(lenghts); 
		  int total = 0; 
		  int partial =0; 
		  int index=0; 
		  for(int i : lenghts){ 
		   total+=i; 
		  } 
		  //Defined accordingly to Assemlathon 2 definition. 
		  //For the Assemplaton 1 definition "<="---> "<"; 
		  while(index < lenghts.size() && partial+lenghts.get(index) <= (total*n)){ 
		   partial += lenghts.get(index); 
		   index++; 
		  } 
		   
		  return lenghts.get(index); 
	} 
	
	
	
	public void windowSlideContigList() throws FileNotFoundException, UnsupportedEncodingException {
		
		findMaxWindows();

		int nbZeroVals=0;
		//fill readCounts and count totalDataPoints substracting the zero values
		for (int i = 0; i < contArrList.size(); i++) {//for each contig
			ContigData currentContig = contigsList.get(contArrList.get(i));
			for (int w = 0; w < currentContig.windPos.size(); w++) {//store all readCounts of every window position 
				if(currentContig.windPos.get(w)!=0){
					readCounts[currentContig.windPos.get(w)] += 1;
				}else {
					nbZeroVals++;
					
				}
			}
			totalDataPoints+= (currentContig.windPos.size()-nbZeroVals);
			nbZeroVals=0;
		}
		
		readCounts[0] =readCounts[1] ;//trick to avoid false peak at 0
		
		//fill fit points
		int ind=0;
		intFitPoints=new int[totalDataPoints+nbZeroVals];
		for (int i = 0; i < contArrList.size(); i++) {//for each contig
			ContigData currentContig = contigsList.get(contArrList.get(i));
			for (int w = 0; w < currentContig.windPos.size(); w++) {//store all readCounts of every window position 
				if(currentContig.windPos.get(w)!=0){
					intFitPoints[ind++]=currentContig.windPos.get(w);//coverage per window position
				}			
			}
		}

		insureReadCountsRange();//excludes outliers values (very high and non significant) from readCount

		if(RUN_SECOND_ROUND)  {
			stringSecondRound="_2nd_Round_";		
		}

		//PRINT OUT readCounts
		/*
		PrintWriter writer = new PrintWriter(Ploest.outputFile + "//" + Ploest.projectName+ "//readCounts"+stringSecondRound+".txt", "UTF-8");
		for (int r = 0; r < readCounts.length; r++) {
			writer.println(r + " " + readCounts[r] + "; ");
		}
		writer.println();
		writer.close();
		 */

	}

	private void insureReadCountsRange() {
		int sum=0;
		for (int i = 0; i < readCounts.length; i++) {
			sum+=readCounts[i];
		}
		int space=sum;
		int midPoint;
		for (midPoint=0;midPoint<readCounts.length;midPoint++){
			space-=readCounts[midPoint];
			if (space<(sum*0.03))break;//select range that takes 95% of all datapoints in readcounts
		}

		//redo readCounts with proper range if necesary
		int SAFE_RANGE=(int) Math.ceil(midPoint*1.2);
		if((readCounts.length-midPoint)>SAFE_RANGE){

			System.out.println("keeping original ReadCounts range? no, reaffecting ReadCounts range. Old:"+readCounts.length+ " new:"+SAFE_RANGE);
			totalDataPoints=0;
			readCounts = new int[(int)SAFE_RANGE];
			int nbNullVals=0;
			//fill readCounts without zero values nor extralarge values
			for (int i = 0; i < contArrList.size(); i++) {//for each contig
				ContigData currentContig = contigsList.get(contArrList.get(i));
				for (int w = 0; w < currentContig.windPos.size(); w++) {//store all readCounts of every window position 
					if (currentContig.windPos.get(w)<SAFE_RANGE && currentContig.windPos.get(w)!=0){//but only if they are inside the new range
						readCounts[currentContig.windPos.get(w)] += 1;//readCounts[0] is goign to be kept zero since it doesn't add info to our problem
					}else if (currentContig.windPos.get(w)>=SAFE_RANGE){
						currentContig.windPos.set(w, null);
						nbNullVals++;
					}
				}

				totalDataPoints+= (currentContig.windPos.size()-nbNullVals);
				nbNullVals=0;
			}


			//fill fit points (dataset to be fitted)
			int ind=0;
			intFitPoints=new int[SamParser.totalDataPoints];
			for (int i = 0; i < contArrList.size(); i++) {//for each contig
				ContigData currentContig = contigsList.get(contArrList.get(i));
				for (int w = 0; w < currentContig.windPos.size(); w++) {//store all readCounts of every window position 
					if(currentContig.windPos.get(w)!=null){//extra big values have been set to zero in previoous loop
						intFitPoints[ind++]=currentContig.windPos.get(w);//coverage per window position 
					}			
				}

			}
			



		}else System.out.println("Keeping original ReadCounts range? yes, keeping:"+readCounts.length+ " instead of changing to mid point :"+SAFE_RANGE);
		//System.out.println("sum:"+sum+ " midPoint:"+midPoint+ " readCounts size:"+readCounts.length);
	}

	//to correctly range the y axis of readCounts
	public void findMaxWindows() throws FileNotFoundException, UnsupportedEncodingException {

		int currentMax=0;
		contArrList = new ArrayList<String>(contigsList.keySet());
		//int sumOfAllWindCoverages=0;

		for (int i = 0; i < contArrList.size(); i++) {
			try {
				ContigData currentContig = contigsList.get(contArrList.get(i));
				currentMax=currentContig.windPos(windowLength);				
				if (currentMax>readsDistributionMaxCoverage)readsDistributionMaxCoverage= currentMax;
			} catch (FileNotFoundException | UnsupportedEncodingException e) {
				e.printStackTrace();
			}
			//totalDataPoints+=maxWindows+1;
		}
		readCounts = new int[(int)Math.ceil(readsDistributionMaxCoverage*1.05)];
		System.out.println("findMaxWindows maxWindows:"+readsDistributionMaxCoverage);
	}


	public void printContigList() {

		for (int i = 0; i < contArrList.size(); i++) {
			System.out.println(contArrList.get(i));
			for (int k = 0; k < 200; k++) {// here I print oby the first 200
				// positions
				System.out.print(contigsList.get(contArrList.get(i)).startPos[k] + "; ");
			}
			System.out.println();
		}
	}



	
	
}
