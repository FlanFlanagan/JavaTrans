package javaTrans;

import java.util.ArrayList;

public class Main {
	public static void main(String[] args){
		ArrayList<Isotope> projectIsos = new ArrayList<Isotope>();
		IsoReader.readIsos(projectIsos);
		projectIsos.get(0).printInfo();
	}
}
