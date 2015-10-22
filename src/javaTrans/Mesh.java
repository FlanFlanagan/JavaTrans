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
	
	double sizeX;
	Region region;
	
	ArrayList<ArrayList<ArrayList<Double>>> sourceTerm = new ArrayList<ArrayList<ArrayList<Double>>>();
	ArrayList<ArrayList<ArrayList<ArrayList<Double>>>> flux = new ArrayList<ArrayList<ArrayList<ArrayList<Double>>>>();
	ArrayList<ArrayList<ArrayList<ArrayList<Double>>>> fluxTotal = new ArrayList<ArrayList<ArrayList<ArrayList<Double>>>>();
	ArrayList<ArrayList<ArrayList<ArrayList<Double>>>> fluxPTotal = new ArrayList<ArrayList<ArrayList<ArrayList<Double>>>>();
	ArrayList<Double> energyFlux = new ArrayList<Double>();
	double[][][] sourceTArray;
	double[][][][] fluxArray;
	double[][][][] fluxTArray;
	double[][][][] fluxPTArray;
	double FinalFlux;
	
	Mesh(Constants conts, double size, int pos, double xLower2, Region region){
		/**TODO double check the ArrayLists here. It might be easier to just store i-1/2, i, 1+1/2 */
		buildFluxes(conts);
		buildSource(conts);
		source2Array();
		this.fluxArray = flux2Array(this.flux);
		this.fluxTArray = flux2Array(this.fluxTotal);
		this.fluxPTArray = flux2Array(this.fluxPTotal);
		this.sizeX = size;
		this.xLower = xLower2 + pos * size;
		this.xPosition = this.xLower + size / 2.;
		this.xUpper = this.xLower + size;
		this.region = region;
	}
	
	double[][][][] flux2Array(ArrayList<ArrayList<ArrayList<ArrayList<Double>>>> arrayList){
		double[][][][] tempArray = new double[arrayList.size()][][][];
		for(int i = 0, ii = arrayList.size(); i < ii; i++){
			double[][][] tempArray1 = new double[arrayList.get(i).size()][][];
			for(int j = 0, ji = arrayList.get(i).size(); j < ji; j++){
				double[][] tempArray2 = new double[arrayList.get(i).get(j).size()][];
				for(int k = 0, ki = arrayList.get(i).get(j).size(); k < ki; k++){
					double[] tempArray3 = new double[arrayList.get(i).get(j).get(k).size()];
					for(int m = 0, mi = arrayList.get(i).get(j).get(k).size(); m < mi; m++){
						tempArray3[m] = arrayList.get(i).get(j).get(k).get(m);
					}
					tempArray2[k] = tempArray3;
				}
				tempArray1[j] = tempArray2;
			}
			tempArray[i] = tempArray1;
		}
		return tempArray;
	}
	
	void source2Array(){
		double[][][] tempArray = new double[this.sourceTerm.size()][][];
		for(int i = 0, si = this.sourceTerm.size(); i < si; i++){
			double[][] tempArray1 = new double[this.sourceTerm.get(i).size()][];
			for(int j = 0, ji = this.sourceTerm.get(i).size(); j < ji; j++){
				double[] tempArray2 = new double[this.sourceTerm.get(i).get(j).size()];
				for(int k = 0, ki = this.sourceTerm.get(i).get(j).size(); k < ki; k++){
					tempArray2[k] = this.sourceTerm.get(i).get(j).get(k);
				}
				tempArray1[j] = tempArray2;
			}
			tempArray[i] = tempArray1;
		}
		this.sourceTArray = tempArray;
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
					for(int bin = 0; bin < conts.eBins; bin++){
						binsCurrent.add(0.);
						binsOld.add(0.);
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
				for(int bin = 0; bin < conts.eBins; bin++){
					binsCurrent.add(conts.source / ((conts.ordinates * conts.wew[k])/2));
				}
				sourceCurrent.add(binsCurrent);
			}
			this.sourceTerm.add(sourceCurrent);	
		}
	}
	
	void setFlux(int edge, int mew, int e, double flux){
		this.fluxArray[edge][0][mew][e] = flux;
	}
	
	void zeroFlux(){
		for(int i = 0, ii = this.fluxArray.length; i < ii; i++){
			for(int j = 0, ji = this.fluxArray[i].length; j < ji; j++){
				for(int k = 0, ki = this.fluxArray[i][j].length; k < ki; k++){
					for(int m = 0, mi = this.fluxArray[i][j][k].length; m < mi; m++){
						this.fluxArray[i][j][k][m] = 0.;
					}
				}
			}
		}
	}
	
	void calcFluxCenterXR(int mew, int e, Constants conts){
		double left = this.sizeX * this.sourceTArray[0][mew][e];
		double mewNum = 2 * conts.mew[mew] * this.fluxArray[0][0][mew][e];
		double denom = 2 * conts.mew[mew] + this.sizeX * this.region.totalXS.get(e);
		this.fluxArray[1][0][mew][e] = (left + mewNum) / denom;
	}
	
	void calcFluxCenterXL(int mew, int e, Constants conts){
		double left = this.sizeX * this.sourceTArray[0][mew][e];
		double mewNum = 2 * conts.mew[mew] * this.fluxArray[2][0][mew][e];
		double denom = -2 * conts.mew[mew] + this.sizeX * this.region.totalXS.get(e);
		this.fluxArray[1][0][mew][e] = (left - mewNum) / denom;
	}
	
	void calcFluxRightX(int mew, int e, Constants conts){
		this.fluxArray[2][0][mew][e] = 2*this.fluxArray[1][0][mew][e] - this.fluxArray[0][0][mew][e];
	}
	
	void calcFluxLeftX(int mew, int e, Constants conts){
		this.fluxArray[0][0][mew][e] = 2*this.fluxArray[1][0][mew][e] - this.fluxArray[2][0][mew][e];		
	}
	
	void zeroSource(){
		for(int i = 0; i < this.sourceTArray.length; i++){
			for(int j = 0; j < this.sourceTArray[i].length; j++){
				for(int k = 0; k < this.sourceTArray[i][j].length; k++){
					this.sourceTArray[i][j][k] = 0.;
				}
			}
		}
	}
	
	void setSource(int mew, int e, double flux){
		this.sourceTArray[0][mew][e] = flux;
	}
	
	void calcScatter(Constants conts){
		for(int l = 0; l < conts.legendre; l++){
			for(int e = 0; e < conts.eBins; e++){
				for(int m = 0, mew = conts.mew.length; m < mew; m++){
					for(int e2 = 0; e2 < conts.eBins; e2++){
						for(int m2 = 0; m2 < mew; m2++){
							this.sourceTArray[0][m][e] += (0.5) * ((2 * l + 1) * conts.lgdr[l][m]) * this.region.sKernalArray[e][e2][l] * (conts.wew[m2] * conts.lgdr[l][m2] * this.fluxArray[1][0][m2][e2]);
						}
					}
				}
			}
		}
	}
	
	void printSource(){
		for(int i = 0; i < this.sourceTArray[0].length; i++){
			for(int j = 0; j < this.sourceTArray[0][i].length; j++){
				System.out.print(String.format("%-20s" , this.sourceTArray[0][i][j] + " "));
			}
			System.out.print("\n");
		}
	}
	
	void sumTotal(int edge, int mew, int e){
		this.fluxTArray[edge][0][mew][e] += this.fluxArray[edge][0][mew][e];
	}
	
	boolean convengenceCheck(Constants conts){
		for(int m = 0, mew = conts.mew.length; m < mew; m++){
			for(int e = 0; e < conts.eBins; e++){
				if(this.fluxArray[1][0][m][e] > conts.convergence * this.fluxTArray[1][0][m][e]){
					return false;
				}
			}
		}
		return true;
	}
	
	void sumTotalEFlux(){
		for(int e = 0; e < region.conts.eBins; e++){
			double totalEFlux = 0;
			for(int m = 0, mew = region.conts.mew.length; m < mew; m++){
				totalEFlux += this.fluxTArray[1][0][m][e];
			}
			this.energyFlux.add(totalEFlux);
		}
	}
	
	void sumTotalFlux(){
		double total = 0;
		for(int e = 0, etot = this.energyFlux.size(); e < etot; e++){
			total += this.energyFlux.get(e);
		}
		System.out.println(this.xPosition + " : " + total);
	}
	
	void printEFlux(int e){
		System.out.println(this.xPosition + " : " + this.energyFlux.get(e));
	}
}
