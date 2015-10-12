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
	
	
	Mesh(Constants conts){
		/**TODO double check the ArrayLists here. It might be easier to just store i-1/2, i, 1+1/2 */
		for(int i = 0; i < conts.dimensions; i++){
			ArrayList<ArrayList<ArrayList<Double>>> edgesCurrent = new ArrayList<ArrayList<ArrayList<Double>>>();
			ArrayList<ArrayList<ArrayList<Double>>> edgesOld = new ArrayList<ArrayList<ArrayList<Double>>>();
			for(int j = 0; j < conts.dimensions; j++){
				ArrayList<ArrayList<Double>> sidesCurrent = new ArrayList<ArrayList<Double>>();
				ArrayList<ArrayList<Double>> sidesOld = new ArrayList<ArrayList<Double>>();
				for(int k = 0; k < conts.ordinates; k++){
					ArrayList<Double> binsCurrent = new ArrayList<Double>();
					ArrayList<Double> binsOld = new ArrayList<Double>();
					for(int bin = 0; bin < conts.groups; bin++){
						binsCurrent.add(0.);
						binsOld.add(0.);
					}
					sidesCurrent.add(binsCurrent);
					sidesOld.add(binsOld);
				}
				edgesCurrent.add(sidesCurrent);
				edgesOld.add(sidesOld);
			}
			this.fluxTotal.add(edgesCurrent);
			this.fluxPTotal.add(edgesOld);
		}
	}
	
	void buildfluxes(){
		
	}
}
