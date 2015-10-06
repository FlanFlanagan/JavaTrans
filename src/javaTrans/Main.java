package javaTrans;

import java.util.ArrayList;

public class Main {
	public static void main(String[] args){
		ArrayList<Region> regions = new ArrayList<Region>();
		regions.add(new Region());
		for(Region reg: regions){
			reg.readIsos();
		}
	}
}
