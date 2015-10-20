package javaTrans;

import java.util.ArrayList;

import org.opensourcephysics.numerics.specialfunctions.Legendre;
import org.opensourcephysics.resources.numerics.specialfunctions;

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
	
	double sizeX;
	Region region;
	
	ArrayList<ArrayList<ArrayList<Double>>> sourceTerm = new ArrayList<ArrayList<ArrayList<Double>>>();
	ArrayList<ArrayList<ArrayList<ArrayList<Double>>>> flux = new ArrayList<ArrayList<ArrayList<ArrayList<Double>>>>();
	ArrayList<ArrayList<ArrayList<ArrayList<Double>>>> fluxTotal = new ArrayList<ArrayList<ArrayList<ArrayList<Double>>>>();
	ArrayList<ArrayList<ArrayList<ArrayList<Double>>>> fluxPTotal = new ArrayList<ArrayList<ArrayList<ArrayList<Double>>>>();
	ArrayList<Double> energyFlux;
	double FinalFlux;
	
	
	Mesh(Constants conts, double size, int xPos, Region region){
		/**TODO double check the ArrayLists here. It might be easier to just store i-1/2, i, 1+1/2 */
		buildFluxes(conts);
		buildSource(conts);
		this.sizeX = size;
		this.xLower = xPos * size;
		this.xPosition = this.xLower + size / 2.;
		this.xUpper = this.xLower + size;
		this.region = region;
	}
	
	void buildFluxes(Constants conts){
		for(int i = 0; i < conts.edges; i++){
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
	void setFlux(int edge, int mew, int e, double flux){
		this.flux.get(edge).get(0).get(mew).set(e, flux);
	}
	
	
	void calcFluxCenterXR(int mew, int e, Constants conts){
		double left = this.sizeX * this.sourceTerm.get(0).get(mew).get(e);
		double mewNum = 2 * conts.mew.get(mew) * this.flux.get(0).get(0).get(mew).get(e);
		double denom = 2 * (conts.mew.get(mew) + this.sizeX * this.region.totalXS.get(e));
		this.flux.get(1).get(0).get(mew).set(e, (left + mewNum) / denom);
	}
	
	void calcFluxCenterXL(int mew, int e, Constants conts){
		double left = this.sizeX * this.sourceTerm.get(0).get(mew).get(e);
		double mewNum = 2 * conts.mew.get(mew) * this.flux.get(2).get(0).get(mew).get(e);
		double denom = -2 * (conts.mew.get(mew) + this.sizeX * this.region.totalXS.get(e));
		this.flux.get(1).get(0).get(mew).set(e, (left - mewNum) / denom);
	}
	
	void calcFluxRightX(int mew, int e, Constants conts){
		/*double left = 2*this.meshPoints.get(m).flux.get(1).get(0).get(mew).get(e);
		double right = this.meshPoints.get(m).flux.get(0).get(0).get(mew).get(e);
		this.meshPoints.get(m).flux.get(2).get(0).get(mew).set(e, ((2 * left) - right));*/
		this.flux.get(2).get(0).get(mew).set(e, ((2*this.flux.get(1).get(0).get(mew).get(e)) - this.flux.get(0).get(0).get(mew).get(e)));
	}
	
	void calcFluxLeftX(int mew, int e, Constants conts){
		/*double left = 2*this.meshPoints.get(m).flux.get(1).get(0).get(mew).get(e);
		double right = this.meshPoints.get(m).flux.get(2).get(0).get(mew).get(e);
		this.meshPoints.get(m).flux.get(0).get(0).get(mew).set(e, ((2 * left) - right));*/
		this.flux.get(0).get(0).get(mew).set(e, ((2*this.flux.get(1).get(0).get(mew).get(e)) - this.flux.get(2).get(0).get(mew).get(e)));		
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
						for(int m2 = 0; m2 < conts.mew.size(); m2++){
							/*double leg1 = 0.5 * (2 * l + 1) * Legendre.evaluate(l, conts.mew.get(m));
							double leg2 = this.region.sKernal.get(e).get(e2).get(l) * conts.wew.get(m2) * Legendre.evaluate(l, conts.mew.get(m2)); 
							double flux = this.flux.get(1).get(0).get(m2).get(e2);
							this.sourceTerm.get(0).get(m).set(e, this.sourceTerm.get(0).get(m).get(e) + leg1 * leg2 * flux);*/
							this.sourceTerm.get(0).get(m).set(e, this.sourceTerm.get(0).get(m).get(e) + (0.5) * (2 * l + 1) * Legendre.evaluate(l, conts.mew.get(m)) * this.region.sKernal.get(e).get(e2).get(l) * conts.wew.get(m2) * Legendre.evaluate(l, conts.mew.get(m2)) * this.flux.get(1).get(0).get(m2).get(e2));
						}
					}
				}
			}
		}
	}
	
	void sumTotal(int edge, int mew, int e){
		this.fluxTotal.get(edge).get(0).get(mew).set(e, this.fluxTotal.get(edge).get(0).get(mew).get(e) + this.flux.get(edge).get(0).get(mew).get(e));
	}
	
	boolean convengenceCheck(Constants conts){
		for(int m = 0; m < conts.mew.size(); m++){
			for(int e = 0; e < conts.eBins; e++){
				if(this.flux.get(1).get(0).get(m).get(e) > conts.convergence * this.fluxTotal.get(1).get(0).get(m).get(e)){
					return false;
				}
			}
		}
		return true;
	}
	
	void sumTotalFlux(){
		for(int e = 0; e < region.conts.eBins; e++){
			double totalEFlux = 0;
			for(int m = 0; m < region.conts.mew.size(); m++){
				totalEFlux += this.fluxTotal.get(1).get(0).get(m).get(e);
			}
			this.energyFlux.add(totalEFlux);
		}
	}
}
