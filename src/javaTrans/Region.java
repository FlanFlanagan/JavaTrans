package javaTrans;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;

import org.jblas.DoubleMatrix;

public class Region {
	double xUpper;
	double xLower;
	
	double yUpper;
	double yLower;
	
	double zUpper;
	double zLower;
	
	ArrayList<Mesh> meshPoints = new ArrayList<Mesh>();
	ArrayList<Isotope> isos = new ArrayList<Isotope>();
	Map<Integer, Double> isotopes = new HashMap<Integer, Double>();	
	String regionType;
	
	double maxMesh;
	double meshSize;
	double meshNumber;
	
	ArrayList<ArrayList<ArrayList<Double>>> sKernal;
	ArrayList<Double> totalXS = new ArrayList<Double>();
	ArrayList<Double> absorbXS = new ArrayList<Double>();
	ArrayList<Double> scatterXS = new ArrayList<Double>();
	ArrayList<Double> chi = new ArrayList<Double>();
	DoubleMatrix fissionM;
	
	double totalC;
	double absorbC;
	double nuFissionC;
	
	ArrayList<Double> finalFlux = new ArrayList<Double>();
	double finalFF;
	double fluxFF;
	
	
}
