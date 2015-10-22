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
	double maxMesh;
	double meshSize;
	int meshNumber;
	double maxXS;
	Constants conts;
	
	ArrayList<ArrayList<ArrayList<Double>>> sKernal = new ArrayList<ArrayList<ArrayList<Double>>>();
	double[][][] sKernalArray; 
	
	ArrayList<Double> totalXS = new ArrayList<Double>();
	ArrayList<Double> absorbXS = new ArrayList<Double>();
	ArrayList<Double> scatterXS = new ArrayList<Double>();
	ArrayList<Double> chi = new ArrayList<Double>();
	
	double totalC;
	double absorbC;
	double nuFissionC;
	
	ArrayList<Double> finalFlux = new ArrayList<Double>();
	double finalFF;
	double fluxFF;
	
	Region(double xLow, double xUp, String type, ArrayList<Isotope> projectIsos, Constants conts, Map<Integer, Double> isotopes){
		this.xLower = xLow;
		this.xUpper = xUp;
		this.regionType = type;
		this.isotopes = isotopes;
		this.conts = conts;
		this.zeroXS();
		this.sumFissile(projectIsos);
		this.buildXS(projectIsos);
		this.skernal2Array();
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
					this.totalXS.set(i, (this.totalXS.get(i) + (this.isotopes.get(iso.name) * iso.totalXS.get(i))));
				}	
			}
		}
	}
	
	void buildChi(ArrayList<Isotope> projectIsos){
		for(Isotope iso: projectIsos){
			if(isotopes.containsKey(iso.name) && iso.fissile == true){
				for(int i = 0; i < iso.eBins; i++){
					this.chi.set(i, (this.chi.get(i) + (this.isotopes.get(iso.name)/this.fissileNumDen * iso.chi.get(i))));
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
							this.sKernal.get(i).get(j).set(lgdr, (this.isotopes.get(iso.name) * iso.sKernal.get(i).get(j).get(lgdr)));
						}
					}
				}	
			}
		}
	}
	
	void skernal2Array(){
		double[][][] tempArray = new double[this.sKernal.size()][][];
		for(int i = 0, si = this.sKernal.size(); i < si; i++){
			double[][] tempArray1 = new double[this.sKernal.get(i).size()][];
			for(int j = 0, ji = this.sKernal.get(i).size(); j < ji; j++){
				double[] tempArray2 = new double[this.sKernal.get(i).get(j).size()];
				for(int k = 0, ki = this.sKernal.get(i).get(j).size(); k < ki; k++){
					tempArray2[k] = this.sKernal.get(i).get(j).get(k);
				}
				tempArray1[j] = tempArray2;
			}
			tempArray[i] = tempArray1;
		}
		this.sKernalArray = tempArray;
		//printSkernalArray();
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
					this.scatterXS.set(i, (this.scatterXS.get(i) + (this.isotopes.get(iso.name) * iso.scatterXS.get(i))));
				}	
			}
		}
	}
	
	void buildAbsorb(ArrayList<Isotope> projectIsos){
		for(Isotope iso: projectIsos){
			if(isotopes.containsKey(iso.name)){
				for(int i = 0; i < iso.eBins; i++){
					this.absorbXS.set(i, (this.absorbXS.get(i) + (this.isotopes.get(iso.name) * iso.absorbXS.get(i))));
				}	
			}
		}
	}
	
	void buildXS(ArrayList<Isotope> projectIsos){
		this.buildTotal(projectIsos);
		this.buildChi(projectIsos);
		this.buildScatter(projectIsos);
		this.buildSkernal(projectIsos);
		this.buildAbsorb(projectIsos);
	}
	
	void calcMaxXS(){
		this.maxXS = this.totalXS.get(0);
		for(int i = 0; i < this.totalXS.size(); i++){
			if(this.totalXS.get(i) > maxXS){
				this.maxXS = this.totalXS.get(i);
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
		this.calcMeshNumber();
	}
	
	void buildMesh(){
		while(this.meshPoints.size() < this.meshNumber){
			this.meshPoints.add(new Mesh(this.conts, this.meshSize, this.meshPoints.size(), this.xLower, this));
		}
	}
		
	void zeroXS(){
		for(int i = 0; i < conts.eBins; i++){
			this.totalXS.add(0.);
			this.absorbXS.add(0.);
			this.scatterXS.add(0.);
			this.chi.add(0.);
			this.finalFlux.add(0.);
		}
		for(int i = 0; i < this.conts.eBins; i++){
			ArrayList<ArrayList<Double>> tempList1 = new ArrayList<ArrayList<Double>>();
			for(int j = 0; j < this.conts.eBins; j++){
				ArrayList<Double> tempList2 = new ArrayList<Double>();
				for(int k = 0; k < this.conts.legendre; k++){
					tempList2.add(0.);
				}
				tempList1.add(tempList2);
			}
			this.sKernal.add(tempList1);
		}
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
	
	boolean convergenceCheck(){
		for(int mesh = 0; mesh < this.meshNumber; mesh++){
			if(!this.meshPoints.get(mesh).convengenceCheck(conts)){
				return false;
			}
		}
		return true;
	}
	
	void setBeamSource(double flux){
		double absorption = this.absorbXS.get(0);
		for(Mesh mesh: this.meshPoints){
			mesh.setSource(0, 0, flux*Math.exp(-absorption*mesh.xPosition));
			//System.out.println(mesh.xPosition + ": " + mesh.sourceTArray[0][0][0]);
		}
	}
}
