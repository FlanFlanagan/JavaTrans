package javaTrans;

import java.util.ArrayList;

public class Main {
	public static void main(String[] args){
		Constants conts = new Constants();
		ArrayList<Isotope> projectIsos = new ArrayList<Isotope>();
		ArrayList<Region> regions = new ArrayList<Region>();
		ProjectBuilder.buildProblem(projectIsos, regions, conts);
		System.out.println("Done.");
	}
}
