import java.io.File;
import java.io.FileNotFoundException;
import java.util.Scanner;

import org.jfree.data.xy.XYDataset;

import dataFitters.GaussianDattaFiter;
import jMEF.PVector;

public class FitGaussian {
	PVector[]         points  ;
	static int MAXNBMODELS=10;

	public FitGaussian(String dataPointsFilePath ){
		
		PVector[] points   = readDataFromFile(dataPointsFilePath);

		//DattaFiter df=new DattaFiter();
		
		for(int nbModels=1;nbModels<MAXNBMODELS;nbModels++){
		
	
		}
	}
	
	
public PVector[] readDataFromFile(String filePath){
		
	    File file = new File(filePath); 
	    PVector[] result=null;
	    try {
	        Scanner scanner = new Scanner(file); 

	        scanner.nextLine(); // to ignore the first line which has the header

	       
	        result=new PVector[7269];
	        int ind=0;
	        while (scanner.hasNextLine()) {
	            String line = scanner.nextLine();
	            if (ind>7250)System.out.println("line:"+line+" ind:"+ind);
	            // Do something with these values
	            PVector curVec=new PVector(1);
	            curVec.array[0]=Double.parseDouble(line);
	            result[ind++]=curVec;
	         
	        }
	        scanner.close();
	    } catch (FileNotFoundException e) {
	        e.printStackTrace();
	    }
	    return result;
	}
	
}
