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
		while(regions.get(0).convergenceCheck() == false || count == 0){
			regions.get(0).sweepLeft();
			regions.get(0).sweepRight();
			regions.get(0).sourceCalc();
			if(count > 5){
				System.out.println(count);
				//System.out.println("Iterations over 100, stopping simulation");
				//break;
			}
			System.out.println(count);
			count++;
		}
		double time2 = System.nanoTime();
		System.out.println("Done. Finished run in: " + String.valueOf((time2-time1)/1E9) + " seconds");
	}
}
