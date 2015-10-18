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
	Map<Integer, Double> isotopes = new HashMap<Integer, Double>();	
	String regionType;
	double fissileNumDen = 0;
	double maxMesh;
	double meshSize;
	int meshNumber;
	double maxXS;
	
	ArrayList<ArrayList<ArrayList<Double>>> sKernal = new ArrayList<ArrayList<ArrayList<Double>>>();
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
		this.zeroXS(conts);
		this.sumFissile(projectIsos);
		this.buildXS(projectIsos);
		this.buildMeshValues(conts);
		this.buildMesh(conts);
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
							this.sKernal.get(i).get(j).set(lgdr, this.sKernal.get(i).get(j).get(lgdr) + (this.isotopes.get(iso.name) * iso.sKernal.get(i).get(j).get(lgdr)));
						}
					}
				}	
			}
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
	
	void calcMaxMesh(Constants conts){
		this.maxMesh = (2 * (1/this.maxXS*1E24) * conts.minMew);
	}
	
	void calcMeshSize(){
		this.meshSize = this.xUpper - this.xLower;
		while(this.meshSize > this.maxMesh){
			meshSize /= 2;
		}
	}
	
	void calcMeshNumber(){
		double xDist = this.xUpper - this.xLower;
		this.meshNumber = (int) (xDist/this.meshSize);
	}
	
	void buildMeshValues(Constants conts){
		this.calcMaxXS();
		this.calcMaxMesh(conts);
		this.calcMeshSize();
		this.calcMeshNumber();
	}
	
	void buildMesh(Constants conts){
		while(this.meshPoints.size() < this.meshNumber){
			this.meshPoints.add(new Mesh(conts, this.meshSize, this.meshPoints.size()));
		}
	}
		
	void zeroXS(Constants conts){
		for(int i = 0; i < conts.eBins; i++){
			this.totalXS.add(0.);
			this.absorbXS.add(0.);
			this.scatterXS.add(0.);
			this.chi.add(0.);
			this.finalFlux.add(0.);
		}
		for(int i = 0; i < conts.eBins; i++){
			ArrayList<ArrayList<Double>> tempList1 = new ArrayList<ArrayList<Double>>();
			for(int j = 0; j < conts.eBins; j++){
				ArrayList<Double> tempList2 = new ArrayList<Double>();
				for(int k = 0; k < conts.legendre; k++){
					tempList2.add(0.);
				}
				tempList1.add(tempList2);
			}
			this.sKernal.add(tempList1);
		}
	}

	public void printinfo() {
		System.out.println("PRINTING REGION INFORMATION");
		System.out.println("X start: " + this.xLower + "X end: " + this.xUpper);
		System.out.println("Region Type: " + this.regionType);
		System.out.println("Fissile Number Density: " + this.fissileNumDen);
		System.out.println(this.meshNumber + " mesh points of " + this.meshSize + "cm");
	}
}
