import java.io.File;
import java.io.FileNotFoundException;
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
	int nbSeq;// nb of sequences in the FileHeader
	Map<String, ContigData> contigsList;// Map of ContigDatas(value) and their
										// name (key)
	static int[] readCounts;
	//List<Integer> points=new ArrayList<Integer>();
	
	int windowLength ;
	List<String> contArrList;
	int maxWindows;
	PloestPlotter plotter;//ploidy estimation  plott and pdf gaussian fit data
	static BarChart barchart;

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
		contigsList = new HashMap<String, ContigData>(nbSeq);// Map of
																// ContigDatas(value)
																// and their
																// name (key)

		// fill the contigsList
		for (int i = 0; i < nbSeq; i++) {
			contigsList.put(inputSam.getFileHeader().getSequenceDictionary().getSequence(i).getSequenceName(),
					new ContigData(inputSam.getFileHeader().getSequenceDictionary().getSequence(i).getSequenceName(),
							inputSam.getFileHeader().getSequenceDictionary().getSequence(i).getSequenceLength()));
		}

		inputSam.setValidationStringency(ValidationStringency.SILENT);
		SAMRecordIterator iter = inputSam.iterator();
		PrintWriter writer = new PrintWriter(Ploest.outputFile + "//" + Ploest.projectName+ "//"+Ploest.projectName+"SamParsed.txt", "UTF-8");
		String line = "";
		int i = 1;
		System.out.println("Analyzing "+contigsList.size()+" contigs");
		String refName = "";// for debugging a bad line in the .sam file
		int alStart;

		while (iter.hasNext()) {//iterates the sam file

			SAMRecord rec = iter.next();
			refName = rec.getReferenceName();
			alStart = rec.getAlignmentStart();
			line = (i + " " + rec.getReadName() + " " + refName + " " + alStart + " " + rec.getAlignmentEnd() + " "
					+ rec.getReadLength() + " " + rec.getMappingQuality() + " ;");
			writer.println(line);

			try {
				contigsList.get(refName).setPos(alStart);//storing the starting positions in the corresponing contig
			} catch (Exception e) {
			}
			i++;

		}
		iter.close();
		inputSam.close();
		writer.close();
		windowSlideContigList();
		barchart = new BarChart(readCounts);
		PloestPlotter plotter = new PloestPlotter(contigsList,maxWindows);
	}

	

	



	private void barchartWithFit(PloestPlotter plotter ) {
		for (int r=0;r<plotter.gmPDF.length;r++){
			barchart.BarChartWithFit(plotter.gmPDF[r],r);
		}		
	}

	
	public void cleanContigList(){//To clean outliers -- TURNED OUT TO BE A BAD iDEA
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
			for (int w = 0; w < currentContig.windPos.length; w++) {//check all readCounts of every window position 
				if(readCounts[currentContig.windPos[w]]==minYCorrectedValue){
					currentContig.windPos[w]=0;//(int) (currentContig.windPos[w]/10);//instead of putting 0 we reduce keeping the proportions, so that the curve mins are respected
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
			for (int w = 0; w < currentContig.windPos.length; w++) {//store all readCounts of every window position 
				readCounts[currentContig.windPos[w]] ++;
			}
		}
		readCounts[0]=0;//to avoid a false peak at 0 which tells us nothing about the real read counts distribution
	
	}


	public void windowSlideContigList() throws FileNotFoundException, UnsupportedEncodingException {
		findMaxWindows();
		
		
		for (int i = 0; i < contArrList.size(); i++) {//for each contig
			ContigData currentContig = contigsList.get(contArrList.get(i));
			for (int w = 0; w < currentContig.windPos.length; w++) {//store all readCounts of every window position 
				readCounts[currentContig.windPos[w]] += 1;
			}
		}
		
		//cleanContigList();//Clean outliers BAD IDEA
		
		PrintWriter writer = new PrintWriter(Ploest.outputFile + "//" + Ploest.projectName+ "//readCounts.txt", "UTF-8");
		for (int r = 0; r < readCounts.length; r++) {
			writer.print(r + " " + readCounts[r] + "; ");
		}
		writer.println();
		writer.close();

		
	}

	//to correctly range the y axis of readCounts
	public void findMaxWindows() throws FileNotFoundException, UnsupportedEncodingException {
	
		int currentMax=0;
		contArrList = new ArrayList<String>(contigsList.keySet());
		 
		for (int i = 0; i < contArrList.size(); i++) {
			try {
				ContigData currentContig = contigsList.get(contArrList.get(i));
				currentMax=currentContig.windPos(windowLength);				
				if (currentMax>maxWindows)maxWindows= currentMax;
			} catch (FileNotFoundException | UnsupportedEncodingException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
			//totalDataPoints+=maxWindows+1;
		}
		readCounts = new int[(int)Math.ceil(maxWindows*1.05)];
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
