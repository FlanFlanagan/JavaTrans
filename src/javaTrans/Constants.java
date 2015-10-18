package javaTrans;

import java.util.ArrayList;

public class Constants {
	int ordinates = 2;
	int dimensions = 1;
	int groups = 10; /** TODO have this derive from Isotopes */
	double convergence = 0.005;
	int legendre = 9;
	double source = 100000.;
	int eBins;
	
	
	ArrayList<Double> mew = new ArrayList<Double>();
	ArrayList<Double> wew = new ArrayList<Double>();
	double minMew;
	
	ArrayList<Double> mew2 = new ArrayList<Double>(){
		{
			this.add(0.5774);
			this.add(-0.5774);
		}
	};
	ArrayList<Double> wew2 = new ArrayList<Double>(){
		{
			this.add(1.);
			this.add(1.);
		}
	};
	
	ArrayList<Double> mew8 = new ArrayList<Double>(){
		{
			this.add(0.9603);
			this.add(0.7967);
			this.add(0.5255);
			this.add(0.1834);
			this.add(-0.1834);
			this.add(-0.5255);
			this.add(-0.7967);
			this.add(-0.9603);
		}
	};
	ArrayList<Double> wew8 = new ArrayList<Double>(){
		{
			this.add(0.1021);
			this.add(0.2224);
			this.add(0.3137);
			this.add(0.3627);
			this.add(0.3627);
			this.add(0.3137);
			this.add(0.2224);
			this.add(0.1021);
		}
	};
	
	Constants(){
		mewTest();
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
}
