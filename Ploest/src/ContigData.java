import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.io.UnsupportedEncodingException;

public class ContigData {

	String contigName;
	int maxLength;
	int[] startPos;
	int[] windPos;
	int windLength;
	int maxWindows;
	static double COV_RATE=10;// ratio by which the average coverage of each sliding window  is mutiplied
						  //the bigger, the more detailed will be the readsDistribution bar chart
	

	public ContigData(String name, int length) {
		contigName = name;
		maxLength = length;
		startPos = new int[maxLength];
	}

	public void setPos(int p) {
		startPos[p] = startPos[p] + 1;
	}

	public String getContigName() {
		return contigName;
	}

	public int windPos(int wl) throws FileNotFoundException, UnsupportedEncodingException {//computes the windPos vector, storing 
			//at each window position, the number of reads found in that window. Returns the maximum value found at that window
											
		int max=0;
		windLength = wl;
		windPos = new int[(startPos.length / wl) * 2];
		int stIndex = 0;// index in startPos array
		int wdIndex = 0;
		int wsum = 0;// window sum
		PrintWriter writer = new PrintWriter(Ploest.outputFile + "//" + Ploest.projectName+"//windPositionsTest.txt", "UTF-8");
		String line="";
		while (stIndex < (startPos.length-windLength) && (wdIndex<windPos.length-1)) {
			for (int i = 0; i < windLength; i++) {
				wsum += startPos[stIndex++];		
			}
			
			windPos[wdIndex++] =(int) (wsum /(windLength/COV_RATE));// relative average of coverage over
													// the range of the window;
			if (windPos[(wdIndex-1)]>max)max=windPos[(wdIndex-1)];
			
			line="wdIndex:"+(wdIndex-1)+" wsum:"+wsum+" avg:"+(windPos[wdIndex-1]);
			stIndex = stIndex - (wl / 2);// window slides over all positions for
											// a length of wl , but a new window
											// is computed after each wl/2;
			wsum=0;
			writer.println(line);
		}
		writer.close();
		maxWindows=max;
		
		return max;
	}
	
	
	
}
