package javaTrans;

import java.io.IOException;
import java.util.ArrayList;

public class Domain {
	String sourceType = "beam";
	Constants conts;
	ArrayList<Isotope> projectIsos = new ArrayList<Isotope>();
	ArrayList<Region> regions = new ArrayList<Region>();
	double criticality;
	double oldCrit =1.;
    double currentFission = 0;
    double oldFission;
    double leakRight = 0;
    double leakLeft = 0;

	Domain(int ordinates){
		this.conts = new Constants(ordinates);
		ProjectBuilder.buildProblem(projectIsos, regions, conts);
	}
	
	void runProblemReflect(){
		setSource("volumeIsoFlux", 15);
		nufissionCalc();
		updateOldFission();
		for(Region reg: this.regions){
			for(Mesh mesh: reg.meshPoints){
				mesh.sumTotalEFlux(this.conts);
				mesh.sumTotalFlux();
				mesh.calcFissionSource(this.conts);
			}
		}
		int count2 = 0;
		while(convergenceTest(count2, 20)){
			int count1 = 0;
			while(convergenceTest(count1, 150)){
				int count = 0;			
				while(convergenceTest(count, 10)){
					sweepRight();
					for(int m = 0, mew = conts.mew.length; m < mew/2; m++){
						for(int e = 0; e < conts.eBins; e++){
							regions.get(0).meshPoints.get(regions.get(0).meshNumber-1).fluxArray[2][0][mew-m-1][e] = regions.get(0).meshPoints.get(regions.get(0).meshNumber-1).fluxArray[2][0][m][e];
						}
					}
					sweepLeft();
					for(int m = conts.mew.length-1, mew = conts.mew.length; m >= mew/2; m--){
						for(int e = 0; e < conts.eBins; e++){
							regions.get(0).meshPoints.get(0).fluxArray[0][0][mew-m-1][e] = regions.get(0).meshPoints.get(0).fluxArray[0][0][m][e];
						}
					}
					count++;
				}
				for(int i = 0; i < regions.size(); i++){
					regions.get(i).sourceCalc();
				}
				count1++;
			}
			for(Region reg: this.regions){
				for(Mesh mesh: reg.meshPoints){
					mesh.sumTotalEFlux(conts);
					//mesh.sumTotalFlux();
					mesh.calcFissionSource(conts);
				}
			}
			nufissionCalc();
			criticalityCalc();
			updateOldFission();
			updateOldCrit();
			zeroTotalFlux();
			count2++;
		}
		for(int i = 0; i < conts.eBins; i++){
			System.out.println(this.regions.get(0).meshPoints.get(0).energyFlux[i]);
		}
		/*try {
			PlotTools.printMeshPlots(regions, "JustG");
			PlotTools.printPolarCenter(regions.get(0).meshPoints.get(regions.get(0).meshNumber-1), 0, 1, "right");
			PlotTools.printPolarCenter(regions.get(0).meshPoints.get(0), 0, 1, "left");
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}*/
	}
	/**
	 * 
	 */
	void runProblem(){
		int count = 0;
		setSource(this.sourceType, 100000.);
		while(convergenceTest(this.regions, count)){
			sweepLeft();
			sweepRight();
			updateScatterSource();
			regions.get(0).meshPoints.get(0).setFlux(0, 0, 0, 0);
			count++;
		}
		for(Region reg: regions){
			for(Mesh mesh: reg.meshPoints){
				mesh.sumTotalEFlux(reg.conts);
			}
		}
		try {
			PlotTools.printMeshPlots(regions, "JustG");
			PlotTools.printPolarCenter(regions.get(0).meshPoints.get(regions.get(0).meshNumber-1), 9, 1, "right");
			PlotTools.printPolarCenter(regions.get(0).meshPoints.get(0), 9, 1, "left");
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
	
	static void regLeft(Region reg1, Region reg2){
		for(int mew = 0; mew < reg1.conts.mew.length/2; mew++){
			for(int e = 0; e < reg1.conts.eBins; e++){
				reg2.meshPoints.get(0).fluxArray[0][0][mew][e] = reg1.meshPoints.get(reg1.meshNumber-1).fluxArray[2][0][mew][e];
			}
		}	
	}
	
	static void regRight(Region reg1, Region reg2){
		for(int mew = reg1.conts.mew.length/2; mew < reg1.conts.mew.length; mew++){
			for(int e = 0; e < reg1.conts.eBins; e++){
				reg2.meshPoints.get(0).fluxArray[2][0][mew][e] = reg1.meshPoints.get(reg1.meshNumber-1).fluxArray[0][0][mew][e];
			}
		}	
	}
	
	static boolean convergenceTest(ArrayList<Region> regions, int number){
		if(number == 0){
			return true;
		}
		for(Region reg: regions){
			if(!reg.convergenceCheck()){
				return true;
			}
		}
		return false;
	}
	
	static boolean convergenceTest(int number, int iterations){
		if(number < iterations || number == 0){
			return true;
		}
		return false;
	}

	void nufissionCalc(){
		for(Region reg: this.regions){
			for(int eng = 0, bins = this.conts.eBins; eng < bins; eng++){
				for(Mesh mesh: reg.meshPoints){
					for(int m = 0, mew = this.conts.mew.length; m < mew; m++){
						this.currentFission += mesh.fluxTArray[1][0][m][eng] * reg.nufission[eng];
					}
				}
			}
		}
	}

	void updateOldFission(){
		this.oldFission = this.currentFission;
		this.currentFission = 0;
	}
	
	void updateOldCrit(){
		this.oldCrit = this.criticality;
		this.criticality = 0;
	}
	
	void criticalityCalc(){
		System.out.println(this.currentFission + " " + this.oldFission);
		this.criticality = this.oldCrit * (this.currentFission/this.oldFission);
		for(Region reg: this.regions){
			reg.criticality = this.criticality;
		}
		System.out.println(this.criticality);
	}
	
	boolean criticalityConvergence(){
		if(Math.abs(this.criticality-this.oldCrit) > this.conts.convergence*this.criticality){
			return false;
		}
		return true;
	}
	
	void leakageLeft(){
		for(int m = this.conts.mew.length/2, mew = this.conts.mew.length; m < mew; m++){
			for(int e = 0, eng = this.conts.eBins; e < eng; e++){
				this.leakLeft = (0.5) * conts.wew[m] * conts.mew[m] *this.regions.get(0).meshPoints.get(0).fluxTArray[0][0][m][e];
			}
		}
	}
	
	void leakageRight(){
		for(int m = 0, mew = this.conts.mew.length/2; m < mew; m++){
			for(int e = 0, eng = this.conts.eBins; e < eng; e++){
				this.leakRight = (0.5) * conts.wew[m] * conts.mew[m] * this.regions.get(this.regions.size()-1).meshPoints.get(this.regions.get(0).meshNumber-1).fluxTArray[2][0][m][e];
			}
		}
	}
	
	void setBeamSource(double flux){
		for(Region reg: this.regions){
			reg.setBeamSource(flux);
		}
	}
	
	void zeroTotalFlux(){
		for(Region reg: this.regions){
			reg.zeroTotalFlux();
		}
	}
	
	void sweepLeft(){
		for(int i = 0; i < this.regions.size(); i++){
			this.regions.get(i).sweepLeft();
			if(i < this.regions.size()-1){
				regLeft(this.regions.get(i), this.regions.get(i+1));
			}
		}
	}
	
	void sweepRight(){
		for(int i = regions.size()-1; i >= 0; i--){
			regions.get(i).sweepRight();
			if(i > 0){
				regRight(regions.get(i), regions.get(i-1));
			}
		}
	}
	
	void updateScatterSource(){
		for(int i = 0; i < regions.size(); i++){
			regions.get(i).sourceCalc();
		}
	}
	
	void setSource(String sourceType, double flux){
		switch(sourceType){
		case "beam":
			regions.get(0).meshPoints.get(0).setFlux(0, 0, 0, flux /*/ this.conts.wew[0]*/);
			break;
		case "beamAtten":
			regions.get(0).setBeamSource(flux);
			break;
		case "volumeIsoSource":
			for(Region reg: this.regions){
				reg.setVolIsoSource(100);
			}
			break;
		case "volumeIsoFlux":
			for(Region reg: this.regions){
				reg.setVolIsoFlux(100);
			}
			break;
		default:
			System.out.println("Not a source type you dingus.");
			break;
		}
	}
	
}
