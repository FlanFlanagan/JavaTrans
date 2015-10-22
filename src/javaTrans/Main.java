package javaTrans;

import java.util.ArrayList;


public class Main {
	public static void main(String[] args){
		double time1 = System.nanoTime();
		Constants conts = new Constants();
		ArrayList<Isotope> projectIsos = new ArrayList<Isotope>();
		ArrayList<Region> regions = new ArrayList<Region>();
		ProjectBuilder.buildProblem(projectIsos, regions, conts);
		int count = 0;
		regions.get(0).setBeamSource(1000000);
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
			count++;
			System.out.println(count);
		}
		double time2 = System.nanoTime();
		for(Region reg: regions){
			for(int i = 0; i < reg.meshPoints.size(); i++){
				reg.meshPoints.get(i).sumTotalEFlux();
				reg.meshPoints.get(i).sumTotalFlux();
			}
		}
		System.out.println("Done. Finished run in: " + String.valueOf((time2-time1)/1E9) + " seconds and interations = " + count);
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
}
