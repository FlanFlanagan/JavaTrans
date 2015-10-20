package javaTrans;

import java.util.ArrayList;

class Mesh {
	double xLower;
	double xUpper;
	double xPosition;
	double xInt;
	
	double yLower;
	double yUpper;	
	double yPosition;
	double yInt;
	
	double zLower;
	double zUpper;	
	double zPosition;
	double zInt;
	
	ArrayList<ArrayList<ArrayList<Double>>> sourceTerm = new ArrayList<ArrayList<ArrayList<Double>>>();
	ArrayList<ArrayList<ArrayList<ArrayList<Double>>>> flux = new ArrayList<ArrayList<ArrayList<ArrayList<Double>>>>();
	ArrayList<ArrayList<ArrayList<ArrayList<Double>>>> fluxTotal = new ArrayList<ArrayList<ArrayList<ArrayList<Double>>>>();
	ArrayList<ArrayList<ArrayList<ArrayList<Double>>>> fluxPTotal = new ArrayList<ArrayList<ArrayList<ArrayList<Double>>>>();
	ArrayList<Double> energyFlux;
	double FinalFlux;
	
	
	Mesh(Constants conts, double size, int xPos){
		/**TODO double check the ArrayLists here. It might be easier to just store i-1/2, i, 1+1/2 */
		buildFluxes(conts);
		buildSource(conts);
		
		this.xLower = xPos * size;
		this.xPosition = this.xLower + size / 2.;
		this.xUpper = this.xLower + size;
	}
	
	void buildFluxes(Constants conts){
		for(int i = 0; i < conts.dimensions; i++){
			ArrayList<ArrayList<ArrayList<Double>>> edgesCurrent = new ArrayList<ArrayList<ArrayList<Double>>>();
			ArrayList<ArrayList<ArrayList<Double>>> edgesOld = new ArrayList<ArrayList<ArrayList<Double>>>();
			ArrayList<ArrayList<ArrayList<Double>>> edgesFlux = new ArrayList<ArrayList<ArrayList<Double>>>();
			for(int j = 0; j < conts.dimensions; j++){
				ArrayList<ArrayList<Double>> sidesCurrent = new ArrayList<ArrayList<Double>>();
				ArrayList<ArrayList<Double>> sidesOld = new ArrayList<ArrayList<Double>>();
				ArrayList<ArrayList<Double>> sidesFlux = new ArrayList<ArrayList<Double>>();
				for(int k = 0; k < conts.ordinates; k++){
					ArrayList<Double> binsCurrent = new ArrayList<Double>();
					ArrayList<Double> binsOld = new ArrayList<Double>();
					ArrayList<Double> binsFlux = new ArrayList<Double>();
					for(int bin = 0; bin < conts.groups; bin++){
						binsCurrent.add(0.);
						binsOld.add(10000.);
						binsFlux.add(0.);
					}
					sidesCurrent.add(binsCurrent);
					sidesOld.add(binsOld);
					sidesFlux.add(binsFlux);
				}
				edgesCurrent.add(sidesCurrent);
				edgesOld.add(sidesOld);
				edgesFlux.add(sidesFlux);	
			}
			this.fluxTotal.add(edgesCurrent);
			this.fluxPTotal.add(edgesOld);
			this.flux.add(edgesFlux);
		}
	}
	
	void buildSource(Constants conts){
		for(int j = 0; j < conts.dimensions; j++){
			ArrayList<ArrayList<Double>> sourceCurrent = new ArrayList<ArrayList<Double>>();
			for(int k = 0; k < conts.ordinates; k++){
				ArrayList<Double> binsCurrent = new ArrayList<Double>();
				for(int bin = 0; bin < conts.groups; bin++){
					binsCurrent.add(conts.source / ((conts.ordinates * conts.wew.get(k))/2));
				}
				sourceCurrent.add(binsCurrent);
			}
			this.sourceTerm.add(sourceCurrent);	
		}
	}
	
	void zeroSource(){
		for(int i = 0; i < this.sourceTerm.size(); i++){
			for(int j = 0; j < this.sourceTerm.get(i).size(); j++){
				for(int k = 0; k < this.sourceTerm.get(i).get(j).size(); k++){
					this.sourceTerm.get(i).get(j).set(k, 0.);
				}
			}
		}
	}
	
	void calcScatter(Constants conts){
		for(int m = 0; m < conts.mew.size(); m++){
			for(int e = 0; e < conts.eBins; e++){
				for(int l = 0; l < conts.legendre; l++){
					for(int e2 = 0; e2 < conts.eBins; e2++){
						for(int m2 = 0; m2 < conts.eBins; m2++){
							double leg = 0.5 * (2 * l + 1) ;
							this.sourceTerm.get(0).get(m).set(e, this.sourceTerm.get(0).get(m).get(e) + )
							
							//REGIONS[si].meshpoints[imesh].source[0][0][mew1][eng] = REGIONS[si].meshpoints[imesh].source[0][0][mew1][eng] + (0.5)*((2*l+1)*tr1::legendre(l, mew[mew1]))*(REGIONS[si].SKERNAL[eng][i][l])*(Wmew[mew2]*tr1::legendre(l, mew[mew2]) * REGIONS[si].meshpoints[imesh].flux[1][0][0][mew2][i]);
						}
					}
				}
			}
		}
	}
}
