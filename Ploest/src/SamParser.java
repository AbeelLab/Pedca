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
	int nbSeq=0;// nb of sequences in the FileHeader
	Map<String, ContigData> contigsList;// Map of ContigDatas(value) and their
	// name (key)
	static int[] readCounts;

	static PVector[] fitPoints;//all the points of all windows positions (coverages) in all contigs. Used to fit the read counts distribution chart

	static int windowLength ;
	List<String> contArrList;
	static int maxWindows;//max coverage found in all contigs (to be used in x axis reads counts)
	static float maxY=0;//max normalized value in the y axis (to be used in y axis reads counts)
	PloestPlotter plotter;//ploidy estimation  plott and pdf gaussian fit data
	NaivePloestPlotter myploter;
	static BarChart barchart;
	
	static int totalDataPoints=0;//total amount of datapoints (all contigs considered)


	public SamParser(String inputFile, String outputfile)
			throws FileNotFoundException, UnsupportedEncodingException {
		
		this.windowLength=Ploest.windowLength;
		SAMFileReader inputSam = new SAMFileReader(new File(inputFile));
		//readCounts = new int[nc];
		nbSeq = inputSam.getFileHeader().getSequenceDictionary().size();// nb of
		// sequences
		// in
		// the
		// FileHeader
		contigsList = new HashMap<String, ContigData>(nbSeq);// Map of  ContigDatas(value) and their name (key)

		// fill the contigsList
		for (int i = 0; i < nbSeq; i++) {
			contigsList.put(inputSam.getFileHeader().getSequenceDictionary().getSequence(i).getSequenceName(),
					new ContigData(inputSam.getFileHeader().getSequenceDictionary().getSequence(i).getSequenceName(),
							inputSam.getFileHeader().getSequenceDictionary().getSequence(i).getSequenceLength()));
		}

		inputSam.setValidationStringency(ValidationStringency.SILENT);
		SAMRecordIterator iter = inputSam.iterator();
	//PrintWriter writer = new PrintWriter(Ploest.outputFile + "//" + Ploest.projectName+ "//"+Ploest.projectName+"SamParsed.txt", "UTF-8");
	//String line = "";
		int i = 1;
		System.out.println("Analyzing "+contigsList.size()+" contigs");
		String refName = "";// for debugging a bad line in the .sam file
		int alStart;

		while (iter.hasNext()) {//iterates the sam file

			SAMRecord rec = iter.next();
			refName = rec.getReferenceName();
			alStart = rec.getAlignmentStart();
	//line = (i + " " + rec.getReadName() + " " + refName + " " + alStart + " " + rec.getAlignmentEnd() + " "
	//+ rec.getReadLength() + " " + rec.getMappingQuality() + " ;");
	//writer.println(line);

			try {
				contigsList.get(refName).setPos(alStart);//storing the starting positions in the corresponing contig
			} catch (Exception e) {
			}
			i++;

		}
		iter.close();
		inputSam.close();
	//writer.close();
		windowSlideContigList();
		barchart = new BarChart(readCounts);
		maxY=barchart.maxY;
		myploter = new NaivePloestPlotter(contigsList,maxWindows, barchart.normReadCounts);//plotter = new PloestPlotter(contigsList,maxWindows);
		
	}


	public SamParser(File fin, String outputfile)throws IOException {

		this.windowLength=Ploest.windowLength;
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
		plotter = new PloestPlotter(contigsList,maxWindows);
	}




	public void cleanContigList(){//To clean outliers -- TURNED OUT TO BE A BAD iDEA (could still be useful to show windows with irregular coverage -transforms it to zero in the plot-)
		//CLEAN READ COUNTS with non significant presence (outliers)
		//first clean in readCounts
		int ctrDELETEME=0;
		double minThreshold=0.0035;
		int totalSumOfReads=0;

		for (int r = 0; r < readCounts.length; r++) {//get total sum of reads
			totalSumOfReads+=readCounts[r];
		}
		int minYCorrectedValue=1;
		System.out.print("Relative readcounts (!=0) above threshold: ");
		for (int r = 0; r < readCounts.length; r++) {//check if relative readcounts is above threshold
			if(readCounts[r]!=0 && ((double)readCounts[r])/totalSumOfReads < minThreshold){//if a read count !=0 is still below threshold
				System.out.print(" "+r);
				readCounts[r]=minYCorrectedValue;
			}
			readCounts[0]=0;//to avoid a false peak at 0 which tells us nothing about the real read counts distribution
		}
		System.out.println();
		//next clean in contigsList
		ctrDELETEME=0;
		int nbOFChangesDELETEME=0;//nb of changed contigs
		boolean hasChanged=false; 
		for (int i = 0; i < contArrList.size(); i++) {//for each contig
			ContigData currentContig = contigsList.get(contArrList.get(i));//get the contig from contigsList
			for (int w = 0; w < currentContig.windPos.size(); w++) {//check all readCounts of every window position 
				if(readCounts[currentContig.windPos.get(w)]==minYCorrectedValue){
					currentContig.windPos.set(w, 0);//[w]=0;//(int) (currentContig.windPos[w]/10);//instead of putting 0 we reduce keeping the proportions, so that the curve mins are respected
					ctrDELETEME++;
					hasChanged=true;
				}
			}
			if(hasChanged){
				contigsList.replace(currentContig.contigName, currentContig);	
				nbOFChangesDELETEME++;
			}
		}
		System.out.println("Window position with ReadCounts above threshold:"+ctrDELETEME);
		System.out.println("Contigs Changed inContigList:"+nbOFChangesDELETEME);
		//WARNING: DELETE THIS. IT IS A TEST: REDO READCOUNTS

		totalSumOfReads=0;
		contArrList = new ArrayList<String>(contigsList.keySet());
		for (int i = 0; i < contArrList.size(); i++) {//for each contig
			ContigData currentContig = contigsList.get(contArrList.get(i));
			for (int w = 0; w < currentContig.windPos.size(); w++) {//store all readCounts of every window position 
				readCounts[currentContig.windPos.get(w)] ++;
			}
		}
		readCounts[0]=0;//to avoid a false peak at 0 which tells us nothing about the real read counts distribution

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
				}else nbZeroVals++;				
			}
			totalDataPoints+= (currentContig.windPos.size()-nbZeroVals);
			nbZeroVals=0;
		}
		System.out.println(" totalDataPoints="+totalDataPoints);
		//fill fit points
		int ind=0;
		fitPoints=new PVector[SamParser.totalDataPoints];
		for (int i = 0; i < contArrList.size(); i++) {//for each contig
			ContigData currentContig = contigsList.get(contArrList.get(i));
			for (int w = 0; w < currentContig.windPos.size(); w++) {//store all readCounts of every window position 
				if(currentContig.windPos.get(w)!=0){
					//readCounts[currentContig.windPos.get(w)] += 1;
					PVector curVec=new PVector(1);
					curVec.array[0]=currentContig.windPos.get(w);//coverage per window position
					//if(w>(currentContig.windPos.size()-20))System.out.print(" p:"+currentContig.windPos.get(w));
					fitPoints[ind++]=curVec; 

				}			
			}	
			//System.out.println();
		}

		insureReadCountsRange();//excludes outliers values (very high and non significant) form readCount
		//cleanContigList();//Clean outliers BAD IDEA
		//PRINT OUT readCounts
		PrintWriter writer = new PrintWriter(Ploest.outputFile + "//" + Ploest.projectName+ "//readCounts.txt", "UTF-8");
		for (int r = 0; r < readCounts.length; r++) {
			writer.println(r + " " + readCounts[r] + "; ");
		}
		writer.println();
		writer.close();


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
			if (space<(sum*0.01))break;//select range that takes 99% of all datapoints in readcounts
		}

		//redo readCounts with proper range if necesary
		int SAFE_RANGE=2;
		if((readCounts.length-midPoint)>(SAFE_RANGE*midPoint)){

			System.out.println("reaffecting ReadCounts range. Old:"+readCounts.length+ " new:"+(int)Math.ceil(midPoint*SAFE_RANGE));
			totalDataPoints=0;
			readCounts = new int[(int)Math.ceil(midPoint*SAFE_RANGE)];
			int nbNullVals=0;
			//fill readCounts without zero values nor extralarge values
			for (int i = 0; i < contArrList.size(); i++) {//for each contig
				ContigData currentContig = contigsList.get(contArrList.get(i));
				for (int w = 0; w < currentContig.windPos.size(); w++) {//store all readCounts of every window position 
					if (currentContig.windPos.get(w)<midPoint*SAFE_RANGE && currentContig.windPos.get(w)!=0){//but only if they are inside the new range
						readCounts[currentContig.windPos.get(w)] += 1;//readCounts[0] is goign to be kept zero since it doesn't add info to our problem
					}else if (currentContig.windPos.get(w)>=midPoint*SAFE_RANGE){
						currentContig.windPos.set(w, null);
						nbNullVals++;
					}
				}
				totalDataPoints+= (currentContig.windPos.size()-nbNullVals);
				nbNullVals=0;
			}


			//fill fit points (dataset to be fitted)
			int ind=0;
			fitPoints=new PVector[SamParser.totalDataPoints];
			for (int i = 0; i < contArrList.size(); i++) {//for each contig
				ContigData currentContig = contigsList.get(contArrList.get(i));
				for (int w = 0; w < currentContig.windPos.size(); w++) {//store all readCounts of every window position 
					if(currentContig.windPos.get(w)!=null){//extra big values have been set to zero in previoous loop
						PVector curVec=new PVector(1);
						curVec.array[0]=currentContig.windPos.get(w);//coverage per window position
						fitPoints[ind++]=curVec; 
					}			
				}		
			}



		}else System.out.println("keeping original ReadCounts range. keeping:"+readCounts.length+ " instead of midpoint calculated:"+(int)Math.ceil(midPoint*1.1));
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
				if (currentMax>maxWindows)maxWindows= currentMax;
			} catch (FileNotFoundException | UnsupportedEncodingException e) {
				e.printStackTrace();
			}
			//totalDataPoints+=maxWindows+1;
		}
		readCounts = new int[(int)Math.ceil(maxWindows*1.05)];
		System.out.println("findMaxWindows maxWindows:"+maxWindows);
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
