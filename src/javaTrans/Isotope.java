package javaTrans;

import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.stream.Stream;

class Isotope {
	int name;
	double numDensity;
	
	ArrayList<Double> chi;
	ArrayList<Double> nuFission;
	ArrayList<Double> totalXS;
	ArrayList<ArrayList<ArrayList<Double>>> sKernal;
	ArrayList<ArrayList<Double>> fFactor;
	ArrayList<Double> newFFactors;
	ArrayList<Double> sigmaB;
	ArrayList<Double> scatterXS;
	ArrayList<Double> absorbXS;
	
	double change;
	double absorbC;
	double nuFissionC;
	
	void readIsoInformation(){
		Stream<String> lines = null;
		try {
			lines = Files.lines(Paths.get(String.valueOf(name)+".xs"), StandardCharsets.UTF_8);
		} catch (IOException e) {
			e.printStackTrace();
		}
		lines.forEach(System.out::println);
	}
	
	void 
}
