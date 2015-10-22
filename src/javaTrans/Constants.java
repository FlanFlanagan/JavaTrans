package javaTrans;

import org.opensourcephysics.numerics.specialfunctions.Legendre;

public class Constants {
	int ordinates = 8;
	int dimensions = 1;
	double convergence = 0.0001;
	int legendre = 9;
	int edges = 3;
	double source = 0.;
	int eBins;
	
	
	Double[] mew;
	Double[] wew;
	double minMew;
	
	Double[][] lgdr;
	
	Double[] mew2 = {0.5774, -0.5774};
	Double[] wew2 = {1., 1.};
	Double[] mew8 = {0.9603, 0.7967, 0.5255, 0.1834, -.01834, -0.5255, -0.7967, -0.9603};
	Double[] wew8 = {0.1021 , 0.2224, 0.3137, 0.3627, 0.3627, 0.3137, 0.2224, 0.1021};
	
	Constants(){
		mewTest();
		buildLgdr();
	}
	
	void mewTest(){
		if(ordinates == 2){
			mew = mew2;
			wew = wew2;
			minMew = 0.5571;
		} else {
			mew = mew8;
			wew = wew8;
			minMew = 0.1834;
		}
	}
	
	void buildLgdr(){
		Double[][] tempArray = new Double[this.legendre][this.mew.length];
		for(int l = 0; l < legendre; l++){
			for(int m = 0, mew = this.mew.length; m < mew; m++){
				tempArray[l][m] = Legendre.evaluate(l, this.mew[m]);
			}
		}
		this.lgdr = tempArray;
	}
}
