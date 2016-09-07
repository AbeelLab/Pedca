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
		System.out.println(contigsList.size()+" Contigs :" + contigsList.keySet());
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
		PloestPlotter plotter = new PloestPlotter(contigsList,maxWindows);
		
		barchartWithFit(plotter);
		//FitGaussian fitGauss=new FitGaussian(plotter.data);
		// printContigList();

	}

	

	private void barchartWithFit(PloestPlotter plotter ) {
		BarChart barchart = new BarChart(readCounts);
		//System.out.println("plotter.gmPDF size:"+plotter.gmPDF.length);
		for (int r=0;r<plotter.gmPDF.length;r++){
			//System.out.println("r:"+r);
			barchart.BarChartWithFit(plotter.gmPDF[r],r);
		}
		
		
		
	}



	public void windowSlideContigList() throws FileNotFoundException, UnsupportedEncodingException {
		findMaxWindows();
		int totalSumOfReads=0;
		for (int i = 0; i < contArrList.size(); i++) {//for each contig
			ContigData currentContig = contigsList.get(contArrList.get(i));
			for (int w = 0; w < currentContig.windPos.length; w++) {//store all readCounts of every window position 
				readCounts[currentContig.windPos[w]] += 1;
			}
		}
		/*
		//CLEAN READ COUNTS with non significant presence (outliers)
		//first clean in readCounts
		for (int r = 0; r < readCounts.length; r++) {
			totalSumOfReads+=readCounts[r];
		}
	
		for (int r = 0; r < readCounts.length; r++) {
			if(((double)readCounts[r])/totalSumOfReads<0.0035){
				readCounts[r]=0;
			}
		}
		//next clean in contigsList
		for (int i = 0; i < contArrList.size(); i++) {//for each contig
			ContigData currentContig = contigsList.get(contArrList.get(i));//get the contig from contigsList
			for (int w = 0; w < currentContig.windPos.length; w++) {//check all readCounts of every window position 
				if(readCounts[currentContig.windPos[w]]==0){
					currentContig.windPos[w]=0;					
				}
			}
			contigsList.replace(currentContig.contigName, currentContig);
				
		}
	
		*/
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
		System.out.println("readCounts SIZE:"+readCounts.length);
		
		 
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
