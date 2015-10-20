package javaTrans;

import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Vector;
import java.util.stream.Stream;

class Isotope {
	
	int name;
	int eBins;
	int legdrNum;
	boolean fissile = false;
	boolean background = false;
	
	ArrayList<Double> chi = new ArrayList<Double>();
	ArrayList<Double> nuFission = new ArrayList<Double>();
	ArrayList<Double> totalXS = new ArrayList<Double>();
	ArrayList<ArrayList<ArrayList<Double>>> sKernal = new ArrayList<ArrayList<ArrayList<Double>>>();
	ArrayList<ArrayList<Double>> fFactor = new ArrayList<ArrayList<Double>>();
	ArrayList<Double> newFFactors = new ArrayList<Double>();
	ArrayList<Double> sigmaB = new ArrayList<Double>();
	ArrayList<Double> scatterXS = new ArrayList<Double>();
	ArrayList<Double> absorbXS = new ArrayList<Double>();
	
	double change;
	double absorbC;
	double nuFissionC;
	
	
	Isotope(String string){
		this.name = Integer.parseInt(string.split(" ")[0]);
		System.out.println("Building Isotope information for: " + this.name);
		readIsoInformation();
		if(nuFission.size() > 0){
			this.fissile = true;
		}
		if(Integer.parseInt(string.split(" ")[1]) == 1){
			this.background = true;
		}
		buildScatter();
		buildAbsorb();
		System.out.println("Finished building Isotope information for: " + this.name);
	}
	
	void readIsoInformation(){
		Stream<String> lines = null;
		try {
			lines = Files.lines(Paths.get("Isos/" + String.valueOf(this.name) + ".xs"), StandardCharsets.UTF_8);
		} catch (IOException e) {
			e.printStackTrace();
		}
		Vector<String> strings = new Vector<String>();
		lines.forEach((line) -> strings.add(line));
		buildIsoInformation(strings);
		
	}
	
	void buildIsoInformation(Vector<String> strings){
		//removes any empty lines
		for(int i = 0; i < strings.size(); i++){
			strings.set(i, strings.get(i).trim().replaceAll("\\s+", " "));
			if(strings.get(i).trim().equalsIgnoreCase("")){
				strings.remove(i);
				i--;
			}
		}
		for(int i = 0; i < strings.size(); i++){
			if(strings.get(i).equalsIgnoreCase("total")){
				this.eBins = readTotal(strings, i+1);
			}
			switch (strings.get(i).toLowerCase()){
			case "ffactor":
				readFFactor(strings, i+1);
				break;
			case "chi":
				readChi(strings, i+1);
				break;
			case "nufission":
				readNuFission(strings, i+1);
				break;
			case "skernel":
				readSkernal(strings, i+1);
				break;
			default:
				break;
			}
		}
	}
	
	int readTotal(Vector<String> strings, int i){
		for(;i < strings.size(); i++){
			String[] tempStrings = strings.get(i).split(" ");
			if(tempStrings.length > 1){
				this.totalXS.add(Double.parseDouble(tempStrings[1]));
			} else {
				Collections.reverse(this.totalXS);
				System.out.println("Total length: " + this.totalXS.size());
				return this.totalXS.size();
			}
		}
		Collections.reverse(this.totalXS);
		return this.totalXS.size();
	}
	
	void readFFactor(Vector<String> strings, int i){
		int limit = i + eBins;
		for(;i < limit; i++){
			ArrayList<Double> tempArray = new ArrayList<Double>();
			String[] tempStrings = strings.get(i).split(" ");
			for(int j = 1; j <tempStrings.length; j++){
				tempArray.add(Double.parseDouble(tempStrings[j]));
			} 
			this.fFactor.add(tempArray);
		}
		Collections.reverse(this.fFactor);
	}
	
	void readChi(Vector<String> strings, int i){
		if(strings.get(i).trim().equalsIgnoreCase("nufission")){
			return;
		}
		int limit = i + eBins;
		for(;i < limit; i++){
			String[] tempStrings = strings.get(i).split(" ");
			this.chi.add(Double.parseDouble(tempStrings[1]));
		}
		Collections.reverse(this.chi);
	}
	
	void readNuFission(Vector<String> strings, int i){
		if(strings.get(i).trim().equalsIgnoreCase("skernel")){
			return;
		}
		int limit = i + eBins;
		for(;i < limit; i++){
			String[] tempStrings = strings.get(i).split(" ");
			this.nuFission.add(Double.parseDouble(tempStrings[1]));
		}
		Collections.reverse(this.nuFission);
	}
	
	void readSkernal(Vector<String> strings, int i){
		int limit = i + eBins;
		//size up group to group lists
		for(int j = 0; j < eBins; j++){
			ArrayList<ArrayList<Double>> tempArray = new ArrayList<ArrayList<Double>>();
			for(int k = 0; k < eBins; k++){
				tempArray.add(new ArrayList<Double>());
			}
			this.sKernal.add(tempArray);
		}
		//add known legendre moments
		for(;i < limit; i++){
			String[] tempStrings = strings.get(i).split(" ");
			int toGroup = eBins - Integer.parseInt(tempStrings[0]);
			int fromGroup = eBins - Integer.parseInt(tempStrings[1]);
			ArrayList<Double> legendre = new ArrayList<Double>();
			for(int j = 2; j < tempStrings.length; j++){
				legendre.add(Double.parseDouble(tempStrings[j]));
			}
			this.legdrNum = legendre.size();
			this.sKernal.get(toGroup).set(fromGroup, legendre);
		}
		// fill in empty legendre moments
		for(int k = 0; k < this.sKernal.size(); k++){
			for(int j = 0; j < this.sKernal.get(k).size(); j++){
				while(this.sKernal.get(k).get(j).size() < this.legdrNum){
					this.sKernal.get(k).get(j).add(0.);
				}
			}
		}
	}
	
	void buildScatter(){
		for(int i = 0; i < this.eBins; i++){
			this.scatterXS.add(0.);
		}
		for(int i = 0; i < this.eBins; i++){
			for(int j = 0; j < this.eBins; j++){
				this.scatterXS.set(i, this.scatterXS.get(i) + this.sKernal.get(j).get(i).get(0));
			}
		}
	}
	
	void buildAbsorb(){
		for(int i = 0; i < this.eBins; i++){
			this.absorbXS.add(0.);
		}
		for(int i = 0; i < this.eBins; i++){
			for(int j = 0; j < this.eBins; j++){
				this.absorbXS.set(i, this.totalXS.get(i) - this.scatterXS.get(i));
			}
		}
	}
	
	void printInfo(){
		System.out.println("\nISOTOPE INFORMATION FOR " + this.name);
		System.out.println("TOTAL: " + this.totalXS);
		System.out.println("CHI: " + this.chi);
		System.out.println("NUFISSION: " + this.nuFission);
		System.out.println("FFACTORS: ");
		for(int i = 0; i < this.fFactor.size(); i++){
			System.out.println(this.fFactor.get(i));
		}
		System.out.println("SKERNAL: ");
		for(int i = 0; i < this.sKernal.size(); i++){
			for(int j = 0; j < this.sKernal.get(i).size(); j++){
				System.out.print(String.format("%-12s" , this.sKernal.get(i).get(j).get(0)));
			}
			System.out.print('\n');
		}
		System.out.println(this.name + " fissionable status: " + this.fissile);
		System.out.println(this.name + " background status: " + this.background + "\n");
	}
}
