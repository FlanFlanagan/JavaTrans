package javaTrans;

import java.util.Arrays;

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
	
	double[] energyFlux;
	double[][][] sourceTArray;
	double[][][][] fluxArray;
	double[][][][] fluxTArray;
	double[][][][] fluxPTArray;
	double finalFlux;
	
	Mesh(Constants conts, double size, int pos, double xLower2, Region region){
		source2Array(conts);
		this.fluxArray = flux2Array(conts);
		this.fluxTArray = flux2Array(conts);
		this.fluxPTArray = flux2Array(conts);
		this.energyFlux = new double[conts.eBins];
		this.sizeX = size;
		this.xLower = xLower2 + pos * size;
		this.xPosition = this.xLower + size / 2.;
		this.xUpper = this.xLower + size;
		this.region = region;
	}
	
	double[][][][] flux2Array(Constants conts){
		double[][][][] tempArray = new double[conts.edges][conts.dimensions][conts.ordinates][conts.eBins];
		return tempArray;
	}
	
	void source2Array(Constants conts){
		this.sourceTArray = new double[conts.dimensions][conts.ordinates][conts.eBins];
	}
	
	void setFlux(int edge, int mew, int e, double flux){
		this.fluxArray[edge][0][mew][e] = flux;
	}
	
	void zeroFlux(){
		for(int i = 0, ii = this.fluxArray.length; i < ii; i++){
			for(int j = 0, ji = this.fluxArray[i].length; j < ji; j++){
				for(int k = 0, ki = this.fluxArray[i][j].length; k < ki; k++){
					Arrays.fill(this.fluxArray[i][j][k], 0.);
				}
			}
		}
	}
	
	void zeroTotalFlux(){
		for(int i = 0, ii = this.fluxTArray.length; i < ii; i++){
			for(int j = 0, ji = this.fluxTArray[i].length; j < ji; j++){
				for(int k = 0, ki = this.fluxTArray[i][j].length; k < ki; k++){
					for(int m = 0, mi = this.fluxTArray[i][j][k].length; m < mi; m++){
						this.fluxTArray[i][j][k][m] = 0.;
					}
				}
			}
		}
	}
	
	void calcFluxCenterXR(int mew, int e, Constants conts){
		double left = this.sizeX * this.sourceTArray[0][mew][e];
		double mewNum = 2 * conts.mew[mew] * this.fluxArray[0][0][mew][e];
		double denom = 2 * conts.mew[mew] + this.sizeX * this.region.totalXS[e];
		this.fluxArray[1][0][mew][e] = (left + mewNum) / denom;
	}
	
	void calcFluxCenterXL(int mew, int e, Constants conts){
		double left = this.sizeX * this.sourceTArray[0][mew][e];
		double mewNum = 2 * conts.mew[mew] * this.fluxArray[2][0][mew][e];
		double denom = -2 * conts.mew[mew] + this.sizeX * this.region.totalXS[e];
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
				Arrays.fill(this.sourceTArray[i][j], 0.);;
			}
		}
	}
	
	void setSource(int mew, int e, double flux){
		this.sourceTArray[0][mew][e] = flux;
	}
	
	void calcScatter(Constants conts){
		double phiL;
		double lValue;
		double skernal;
		for(int l = 0; l < conts.legendre; l++){
			lValue = (2 * l + 1);
			for(int e = 0; e < conts.eBins; e++){
				phiL = calcPhiL(conts, l, e);
				for(int e2 = 0; e2 < conts.eBins; e2++){
					skernal = this.region.sKernalArray[e2][e][l];
					if(skernal == 0){
						continue;
					}
					for(int m2 = 0, mew = conts.mew.length; m2 < mew; m2++){
						this.sourceTArray[0][m2][e2] += lValue * conts.lgdr[l][m2] * skernal * phiL;
					}
				}
			}
		}
	}
	
	double calcPhiL(Constants conts, int l, int e){
		double flux = 0;
		for(int m =0, mew = conts.mew.length; m < mew; m++){
			flux += 0.5 * conts.wew[m] * conts.lgdr[l][m] * this.fluxArray[1][0][m][e];
		}
		return flux;
	}
	
	void calcFissionSource(Constants conts){
		double eFlux;
		double tempChi;
		double nufiss;
		for(int m = 0, mew = conts.mew.length; m < mew; m++){
			for(int e = 0; e < conts.eBins; e++){
				for(int e2 = 0; e2 < conts.eBins; e2++){
					eFlux = this.energyFlux[e2];
					for(Isotope iso:this.region.projectIsos){
						if(iso.fissile){
							tempChi = iso.chi.get(e)/this.region.criticality;
							nufiss = iso.nuFission.get(e2);
							if(this.region.isotopes.containsKey(iso.name)){
								this.sourceTArray[0][m][e] +=  tempChi * nufiss * this.region.isotopes.get(iso.name) * eFlux;
							}
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
	
	void sumTotalEFlux(Constants conts){
		for(int e = 0; e < region.conts.eBins; e++){
			double totalEFlux = 0;
			for(int m = 0, mew = region.conts.mew.length; m < mew; m++){
				totalEFlux += this.fluxTArray[1][0][m][e]*conts.wew[m];
			}
			this.energyFlux[e] = totalEFlux;
		}
	}
	
	void sumTotalFlux(){
		this.finalFlux = Arrays.stream(this.energyFlux).sum();
	}
	
	
	void oldTotalUpdate(){
		this.fluxPTArray = this.fluxTArray.clone();
	}
	
	void printEFlux(int e){
		System.out.println(this.xPosition + " : " + this.energyFlux[e]);
	}
}
