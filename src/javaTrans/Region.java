package javaTrans;

import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.stream.Stream;

import org.jblas.DoubleMatrix;

class Region {
	double xUpper;
	double xLower;
	
	double yUpper;
	double yLower;
	
	double zUpper;
	double zLower;
	
	ArrayList<Mesh> meshPoints;
	ArrayList<Isotope> isos;
	
	String regionType;
	
	double maxMesh;
	double meshSize;
	double meshNumber;
	
	ArrayList<ArrayList<ArrayList<Double>>> sKernal;
	ArrayList<Double> totalXS;
	ArrayList<Double> absorbXS;
	ArrayList<Double> scatterXS;
	ArrayList<Double> chi;
	DoubleMatrix fissionM;
	
	double totalC;
	double absorbC;
	double nuFissionC;
	
	ArrayList<Double> finalFlux;
	double finalFF;
	double fluxFF;
	
	void readIsos(){
		Stream<String> lines = null;
		try {
			lines = Files.lines(Paths.get("Isos.txt"), StandardCharsets.UTF_8);
		} catch (IOException e) {
			e.printStackTrace();
		}
		lines.forEach(System.out::println);
	}
	
}
