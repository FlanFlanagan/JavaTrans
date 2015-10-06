package javaTrans;

import java.util.ArrayList;

class Mesh {
	double xLower;
	double xUpper;
	double xPosition;
	
	double yLower;
	double yUpper;	
	double yPosition;
	
	double zLower;
	double zUpper;	
	double zPosition;
	
	ArrayList<ArrayList<ArrayList<Double>>> sourceTerm;
	ArrayList<ArrayList<ArrayList<ArrayList<Double>>>> flux;
	ArrayList<ArrayList<ArrayList<ArrayList<Double>>>> fluxTotal;
	ArrayList<ArrayList<ArrayList<ArrayList<Double>>>> fluxPTotal;
	ArrayList<Double> energyFlux;
	double FinalFlux;	
}
