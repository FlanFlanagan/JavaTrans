package javaTrans;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;

public class Region {
	double xUpper;
	double xLower;
	
	double yUpper;
	double yLower;
	
	double zUpper;
	double zLower;
	
	ArrayList<Mesh> meshPoints = new ArrayList<Mesh>();
	Map<Integer, Double> isotopes = new HashMap<Integer, Double>();	
	String regionType;
	double fissileNumDen = 0;
	double criticality = 1.1;
	double maxMesh;
	double meshSize;
	int meshNumber;
	double maxXS;
	Constants conts;
	ArrayList<Isotope> projectIsos;
	
	double[][][] sKernalArray; 
	double[] totalXS;
	double[] absorbXS;
	double[] scatterXS;
	double[] chi;
	double[] nufission;
	
	double totalC;
	double absorbC;
	double nuFissionC;
	
	double[] finalFlux;
	double finalFF;
	double fluxFF;
	
	Region(double xLow, double xUp, String type, ArrayList<Isotope> projectIsos, Constants conts, Map<Integer, Double> isotopes){
		this.xLower = xLow;
		this.xUpper = xUp;
		this.regionType = type;
		this.isotopes = isotopes;
		this.conts = conts;
		this.projectIsos = projectIsos;
		this.zeroXS();
		this.sumFissile(projectIsos);
		this.buildXS(projectIsos);
		this.buildMeshValues();
		this.buildMesh();
	}
	
	void sumFissile(ArrayList<Isotope> projectIsos){
		for(Isotope iso:projectIsos){
			if(isotopes.containsKey(iso.name)){
				if(iso.fissile){
					fissileNumDen += isotopes.get(iso.name);
				}
			}
		}
	}
		
	void buildTotal(ArrayList<Isotope> projectIsos){
		for(Isotope iso: projectIsos){
			if(this.isotopes.containsKey(iso.name)){
				for(int i = 0; i < iso.eBins; i++){
					this.totalXS[i] += this.isotopes.get(iso.name) * iso.totalXS.get(i);
				}	
			}
		}
	}
	
	void buildChi(ArrayList<Isotope> projectIsos){
		this.chi = new double[this.conts.eBins];
		for(Isotope iso: projectIsos){
			if(isotopes.containsKey(iso.name) && iso.fissile == true){
				for(int i = 0; i < iso.eBins; i++){
					this.chi[i] +=this.isotopes.get(iso.name)/this.fissileNumDen * iso.chi.get(i);
				}	
			}
		}
	}

	void buildSkernal(ArrayList<Isotope> projectIsos){
		for(Isotope iso: projectIsos){
			if(isotopes.containsKey(iso.name)){
				for(int i = 0; i < iso.eBins; i++){
					for(int j = 0; j < iso.eBins; j++){
						for(int lgdr = 0; lgdr < iso.legdrNum; lgdr ++){
							this.sKernalArray[i][j][lgdr] += this.isotopes.get(iso.name) * iso.sKernal.get(i).get(j).get(lgdr);
						}
					}
				}	
			}
		}
	}
	
	void printSkernalArray(){
		for(int i = 0; i < this.sKernalArray.length; i++){
			for(int j = 0; j < this.sKernalArray[i].length; j++){
				System.out.print(String.format("%-20s" , this.sKernalArray[i][j][1]));
			}
			System.out.print('\n');
		}
	}
	
	void buildScatter(ArrayList<Isotope> projectIsos){
		for(Isotope iso: projectIsos){
			if(isotopes.containsKey(iso.name)){
				for(int i = 0; i < iso.eBins; i++){
					this.scatterXS[i] += (this.isotopes.get(iso.name) * iso.scatterXS[i]);
				}	
			}
		}
	}
	
	void buildAbsorb(ArrayList<Isotope> projectIsos){
		for(Isotope iso: projectIsos){
			if(isotopes.containsKey(iso.name)){
				for(int i = 0; i < iso.eBins; i++){
					this.absorbXS[i] += this.isotopes.get(iso.name) * iso.absorbXS[i];
				}	
			}
		}
	}
	
	void buildNuFission(){
		double[] tempArray = new double[this.totalXS.length];
		for(Isotope iso: projectIsos){
			if(isotopes.containsKey(iso.name) && iso.fissile){
				for(int i = 0, max = iso.eBins; i < max; i++){
					tempArray[i] += this.isotopes.get(iso.name) * iso.nuFission.get(i);
				}
			}
		}
		this.nufission = tempArray;
	}
	
	void buildXS(ArrayList<Isotope> projectIsos){
		this.buildTotal(projectIsos);
		this.buildChi(projectIsos);
		this.buildScatter(projectIsos);
		this.buildSkernal(projectIsos);
		this.buildAbsorb(projectIsos);
		this.buildNuFission();
	}
	
	void calcMaxXS(){
		this.maxXS = this.totalXS[0];
		for(int i = 0; i < this.totalXS.length; i++){
			if(this.totalXS[i] > maxXS){
				this.maxXS = this.totalXS[i];
			}
		}
	}
	
	void calcMaxMesh(){
		this.maxMesh = (2 * (1/this.maxXS) * this.conts.minMew);
	}
	
	void calcMeshSize(){
		this.meshSize = this.xUpper - this.xLower;
		while(this.meshSize > this.maxMesh){
			this.meshSize /= 2;
		}
	}
	
	void calcMeshNumber(){
		double xDist = this.xUpper - this.xLower;
		this.meshNumber = (int) (xDist/this.meshSize);
	}
	
	void buildMeshValues(){
		this.calcMaxXS();
		this.calcMaxMesh();
		this.calcMeshSize();
		//this.meshSize = 0.4;
		this.calcMeshNumber();
	}
	
	void buildMesh(){
		while(this.meshPoints.size() < this.meshNumber){
			this.meshPoints.add(new Mesh(this.conts, this.meshSize, this.meshPoints.size(), this.xLower, this));
		}
	}
		
	void zeroXS(){
		this.totalXS = new double[conts.eBins];
		this.absorbXS = new double[conts.eBins];
		this.scatterXS = new double[conts.eBins];
		this.chi = new double[conts.eBins];
		this.finalFlux = new double[conts.eBins];
		/*Arrays.fill(this.totalXS, 0);
		Arrays.fill(this.absorbXS, 0);
		Arrays.fill(this.scatterXS, 0);
		Arrays.fill(this.chi, 0);
		this.finalFlux.add(0.);*/
		this.sKernalArray = new double[this.conts.eBins][this.conts.eBins][this.conts.legendre];
	}

	void printinfo() {
		System.out.println("PRINTING REGION INFORMATION");
		System.out.println("X start: " + this.xLower + "X end: " + this.xUpper + "X");
		System.out.println("Region Type: " + this.regionType);
		System.out.println("Fissile Number Density: " + this.fissileNumDen);
		System.out.println(this.meshNumber + " mesh points of " + this.meshSize + "cm");
	}
	
	void sweepRight(){
		for(int m = 0; m < this.meshNumber; m++){
			for(int mew = 0; mew < this.conts.mew.length/2; mew++){
				for(int e = 0; e < this.conts.eBins; e++){
					this.meshPoints.get(m).calcFluxCenterXR(mew, e, this.conts);
					this.meshPoints.get(m).calcFluxRightX(mew, e, this.conts);
					if(m < meshPoints.size()-1){
						this.meshPoints.get(m+1).setFlux(0, mew, e, this.meshPoints.get(m).fluxArray[2][0][mew][e]);
					}	
					this.meshPoints.get(m).sumTotal(1, mew, e);
				}
			}
		}
	}
	
	void sweepLeft(){
		for(int m = this.meshNumber-1; m >= 0; m--){
			for(int mew = this.conts.mew.length/2; mew < this.conts.mew.length; mew++){
				for(int e = 0; e < this.conts.eBins; e++){
					this.meshPoints.get(m).calcFluxCenterXL(mew, e, this.conts);
					this.meshPoints.get(m).calcFluxLeftX(mew, e, this.conts);
					if(m > 0){
						this.meshPoints.get(m-1).setFlux(2, mew, e, this.meshPoints.get(m).fluxArray[0][0][mew][e]);
					}	
					this.meshPoints.get(m).sumTotal(1, mew, e);
				}
			}
		}
	}
	
	void sourceCalc(){
		for(Mesh mesh: this.meshPoints){
			mesh.zeroSource();
			mesh.calcScatter(this.conts);
		}
		
	}
	
	void zeroFlux(){
		for(Mesh mesh: this.meshPoints){
			mesh.zeroFlux();
		}
	}
	
	void zeroTotalFlux(){
		for(Mesh mesh: this.meshPoints){
			mesh.zeroTotalFlux();
		}
	}
	
	boolean convergenceCheck(){
		for(int mesh = 0; mesh < this.meshNumber; mesh++){
			if(!this.meshPoints.get(mesh).convengenceCheck(conts)){
				return false;
			}
		}
		return true;
	}
	
	void setBeamSource(double flux){
		double absorption = this.absorbXS[0];
		for(Mesh mesh: this.meshPoints){
			mesh.setSource(0, 0, flux*Math.exp(-absorption*mesh.xPosition));
			System.out.println(mesh.xPosition + ": " + mesh.sourceTArray[0][0][0]);
		}
	}
}
