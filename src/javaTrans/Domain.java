package javaTrans;

import java.io.IOException;
import java.util.ArrayList;

public class Domain {
	Constants conts;
	ArrayList<Isotope> projectIsos = new ArrayList<Isotope>();
	ArrayList<Region> regions = new ArrayList<Region>();
	double criticality;
	double oldCrit =1.;
    double currentFission = 0;
    double oldFission = 10;
    double leakRight = 0;
    double leakLeft = 0;

	Domain(int ordinates){
		this.conts = new Constants(ordinates);
		ProjectBuilder.buildProblem(projectIsos, regions, conts);
	}
	
	void runProblemReflect(){
		for(Region reg: regions){
			for(Mesh mesh: reg.meshPoints){
				mesh.setSource(0, 0, 100);
			}
		}
		int count2 = 0;
		while(count2 < 50){
			int count1 = 0;
			while(count1 < 50){
				int count = 0;			
				//while(convergenceTest(regions) != true || count == 0){
				while(count < 50){
					for(int i = 0; i < regions.size(); i++){
						regions.get(i).sweepRight();
						if(i < regions.size()-1){
							regLeft(regions.get(i), regions.get(i+1));
						}
					}
					for(int m = 0, mew = conts.mew.length; m < mew/2; m++){
						for(int e = 0; e < conts.eBins; e++){
							regions.get(0).meshPoints.get(regions.get(0).meshNumber-1).fluxArray[2][0][mew-m-1][e] = regions.get(0).meshPoints.get(regions.get(0).meshNumber-1).fluxArray[2][0][m][e];
						}
					}
					for(int i = regions.size()-1; i >= 0; i--){
						regions.get(i).sweepLeft();
						if(i > 0){
							regRight(regions.get(i), regions.get(i-1));
						}
					}
					for(int m = 0, mew = conts.mew.length; m > mew/2; m--){
						for(int e = 0; e < conts.eBins; e++){
							regions.get(0).meshPoints.get(0).fluxArray[0][0][m][e] = regions.get(0).meshPoints.get(0).fluxArray[0][0][mew-m-1][e];
						}
					}
					count++;
				}
				for(int i = 0; i < regions.size(); i++){
					regions.get(i).sourceCalc();
				}
				count1++;
			}
			for(Region reg: regions){
				for(Mesh mesh: reg.meshPoints){
					mesh.sumTotalEFlux(conts);
					mesh.sumTotalFlux();
					mesh.calcFissionSource(conts);
				}
			}
			nufissionCalc();
			criticalityCalc();
			updateOldFission();
			
			count2++;
		}
		try {
			PlotTools.printMeshPlots(regions, "JustG");
			//PlotTools.printPolar(regions.get(1).meshPoints.get(regions.get(1).meshNumber-1), 9, 1);
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
	void runProblem(){
		int count = 0;
		//regions.get(0).setBeamSource(1000000);
		regions.get(0).meshPoints.get(0).setFlux(0, 0, 0, 1000000);
		//while(convergenceTest(regions) != true || count == 0){
		while(count < 180){
			for(int i = 0; i < regions.size(); i++){
				regions.get(i).sweepLeft();
				if(i < regions.size()-1){
					regLeft(regions.get(i), regions.get(i+1));
				}
			}
			for(int i = regions.size()-1; i >= 0; i--){
				regions.get(i).sweepRight();
				if(i > 0){
					regRight(regions.get(i), regions.get(i-1));
				}
			}
			for(int i = 0; i < regions.size(); i++){
				regions.get(i).sourceCalc();
			}
			regions.get(0).meshPoints.get(0).setFlux(0, 0, 0, 0);
			count++;
		}
		for(Region reg: regions){
			for(Mesh mesh: reg.meshPoints){
				mesh.sumTotalEFlux(reg.conts);
				mesh.sumTotalFlux();
			}
		}
		try {
			PlotTools.printMeshPlots(regions, "JustG");
			//PlotTools.printPolar(regions.get(1).meshPoints.get(regions.get(1).meshNumber-1), 9, 1);
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
	
	static boolean convergenceTest(ArrayList<Region> regions){
		for(Region reg: regions){
			if(!reg.convergenceCheck()){
				return false;
			}
		}
		return true;
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
	}
	
	void criticalityCalc(){
		this.criticality = this.oldCrit * (this.currentFission/this.oldFission);
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
				this.leakLeft = (0.5) * conts.wew[m] * conts.mew[m] * regions.get(0).meshPoints.get(0).fluxTArray[0][0][m][e];
			}
		}
	}
	
	void leakageRight(){
		for(int m = 0, mew = this.conts.mew.length/2; m < mew; m++){
			for(int e = 0, eng = this.conts.eBins; e < eng; e++){
				this.leakRight = (0.5) * conts.wew[m] * conts.mew[m] * regions.get(regions.size()-1).meshPoints.get(regions.get(0).meshNumber-1).fluxTArray[2][0][m][e];
			}
		}
	}
	
	void setBeamSource(double flux){
		for(Region reg: this.regions){
			reg.setBeamSource(flux);
		}
	}
	
	void zeroTotalFlux(){
		for(Region reg: regions){
			reg.zeroTotalFlux();
		}
	}
}
