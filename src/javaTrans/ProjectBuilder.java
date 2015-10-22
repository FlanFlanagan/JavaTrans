package javaTrans;

import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;
import java.util.Vector;
import java.util.stream.Stream;

public class ProjectBuilder {
	
	public static void buildProblem(ArrayList<Isotope> isoArray, ArrayList<Region> regions, Constants conts){
		Stream<String> lines = null;
		try {
			lines = Files.lines(Paths.get("Isos/input.temp"), StandardCharsets.UTF_8);
		} catch (IOException e) {
			e.printStackTrace();
		}
		Vector<String> strings = new Vector<String>();
		lines.forEach((line) -> strings.add(line));
		readInput(isoArray, regions, conts, strings);
		
	}
	
	public static void readInput(ArrayList<Isotope> isoArray, ArrayList<Region> regions, Constants conts, Vector<String> strings){
		for(int i = 0; i < strings.size(); i++){
			if(strings.get(i).trim().split(" ")[0].equalsIgnoreCase("isos")){
				i = readIsos(isoArray, conts, strings, i+1);
			}
			if(strings.get(i).trim().split(" ")[0].equalsIgnoreCase("region")){
				Map<Integer, Double> isoNumDens = new HashMap<Integer, Double>();
				double xLower = Double.parseDouble(strings.get(i).trim().split(" ")[1]);
				double xUpper = Double.parseDouble(strings.get(i).trim().split(" ")[2]);
				String type = strings.get(i).trim().split(" ")[3];
				i = readRegion(isoNumDens, strings, i+1);
				Region reg = new Region(xLower, xUpper, type, isoArray, conts, isoNumDens);
				reg.printinfo();
				regions.add(reg);
			}
		}
	}
	
	public static int readIsos(ArrayList<Isotope> isos, Constants conts, Vector<String> strings, int i){
		while(!strings.get(i).equalsIgnoreCase("isos end")){
			Isotope iso = new Isotope(strings.get(i));
			isos.add(iso);
			//iso.printInfo();
			conts.eBins = iso.eBins;
			i++;
		}
		return i;
	}

	public static int readRegion(Map<Integer, Double> isoNumDens, Vector<String> strings, int i){
		while(!strings.get(i).equalsIgnoreCase("region end")){
			int isoNumber = Integer.parseInt(strings.get(i).split(" ")[0]);
			double numDensity = Double.parseDouble(strings.get(i).split(" ")[1]);
			isoNumDens.put(isoNumber, numDensity*1E-24);
			i++;
		}
		return i;
	}
}
