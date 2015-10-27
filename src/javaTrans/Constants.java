package javaTrans;

import java.util.Arrays;
import org.opensourcephysics.numerics.specialfunctions.Legendre;

public class Constants {
	int ordinates;
	int dimensions = 1;
	double convergence = 0.001;
	int legendre = 9;
	int edges = 3;
	double source = 0.;
	int eBins;
	
	
	double[] mew;
	double[] wew;
	double minMew;
	
	double[][] lgdr;
	
	Constants(int ordinates){
		buildOrdinates(ordinates);
		mewTest();
		buildLgdr();
	}
	
	
	void buildOrdinates(int l){
		double[] tempMews = Legendre.getPolynomial(l).rootsReal();
		Arrays.sort(tempMews);
		double[] mews = new double[tempMews.length];
		for(int i = 0; i < mews.length; i++){
			mews[i] = tempMews[mews.length-1-i];
		}
		double[] wews = new double[mews.length];
		for(int i = 0; i < mews.length; i++){
			double legDer = Legendre.getPolynomial(l).derivative().evaluate(mews[i]);
			double legMinus = Legendre.getPolynomial(l-1).evaluate(mews[i]);
			wews[i] = 2 / (l * legMinus * legDer);
		}
		this.ordinates = mews.length;
		this.mew = mews;
		this.wew = wews;
		/*for(int i = 0; i < this.mew.length; i++){
			System.out.println(mew[i] + " "+ wew[i]);
		}*/
		
	}
	void mewTest(){
		this.minMew = mew[this.ordinates/2-1];
	}
	
	void buildLgdr(){
		double[][] tempArray = new double[this.legendre][this.ordinates];
		for(int l = 0; l < legendre; l++){
			for(int m = 0, mew = this.mew.length; m < mew; m++){
				tempArray[l][m] = Legendre.evaluate(l, this.mew[m]);
			}
		}
		this.lgdr = tempArray;
	}
	
	/*void oridinateReader(){
		double[] tempMew = new double[this.ordinates];
		double[] tempWew = new double[this.ordinates];
		Stream<String> lines = null;
		try {
			lines = Files.lines(Paths.get("Ordinates/" + this.ordinates + ".ord"), StandardCharsets.UTF_8);
		} catch (IOException e) {
			e.printStackTrace();
		}
		Vector<String> strings = new Vector<String>();
		lines.forEach((line) -> strings.add(line));
		for(int i = 0; i < strings.size(); i++){
			strings.set(i, strings.get(i).trim().replaceAll("\\s+", " "));
			String[] tempString = strings.get(i).split(" ");
			tempMew[i] = Double.parseDouble(tempString[0]);
			tempMew[this.ordinates-1-i] = Double.parseDouble(tempString[0])*-1;
			tempWew[i] = Double.parseDouble(tempString[1]);
			tempWew[this.ordinates-1-i] = Double.parseDouble(tempString[1]);
		}
		this.mew = tempMew;
		this.wew = tempWew;
	}*/
}
