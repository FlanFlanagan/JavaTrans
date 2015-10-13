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
	int fissileNumber = 0;
	double maxMesh;
	double meshSize;
	double meshNumber;
	double maxXS;
	
	ArrayList<ArrayList<ArrayList<Double>>> sKernal;
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
	
	Region(double xLow, double xUp, String type, ArrayList<Isotope> projectIsos, Constants conts){
		this.xLower = xLow;
		this.xUpper = xUp;
		this.regionType = type;
		
		this.sumFissile(projectIsos);
		this.buildXS(projectIsos);
		
		this.buildMesh(conts);
	}
	
	void sumFissile(ArrayList<Isotope> projectIsos){
		for(Isotope iso:projectIsos){
			if(isotopes.containsKey(iso.name)){
				if(iso.fissile){
					fissileNumber += 1;
				}
			}
		}
	}
	
	void buildTotal(ArrayList<Isotope> projectIsos){
		for(Isotope iso: projectIsos){
			if(isotopes.containsKey(iso.name)){
				for(int i = 0; i < iso.eBins; i++){
					this.totalXS.set(i, (this.totalXS.get(i) + (isotopes.get(iso.name) * iso.totalXS.get(i))));
				}	
			}
		}
	}
	
	void buildChi(ArrayList<Isotope> projectIsos){
		for(Isotope iso: projectIsos){
			if(isotopes.containsKey(iso.name)){
				for(int i = 0; i < iso.eBins; i++){
					this.chi.set(i, (this.chi.get(i) + (this.isotopes.get(iso.name)/this.fissileNumber * iso.chi.get(i))));
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
		this.maxMesh = (2 * (1/this.maxXS) * conts.minMew);
	}
	
	void calcMeshSize(){
		this.meshSize = this.xUpper - this.xLower;
		while(this.meshSize > this.maxMesh){
			meshSize /= 2;
		}
	}
	
	void buildMesh(Constants conts){
		while(this.meshPoints.size() < this.meshNumber){
			this.meshPoints.add(new Mesh(conts, this.meshSize, this.meshPoints.size()));
		}
	}
	
	void buildMeshValues(Constants conts){
		this.calcMaxXS();
		this.calcMaxMesh(conts);
		this.calcMeshSize();
	}
}
