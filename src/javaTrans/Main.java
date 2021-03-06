package javaTrans;

import java.util.Arrays;


public class Main {
	public static void main(String[] args){
		double time1 = System.nanoTime();
		Domain domain = new Domain(8);
		double time2 = System.nanoTime();
		System.out.println("Domain build time = " + (time2-time1)/1E9);
		//run(domain);
		speedTest(domain);
	}
	
	static void speedTest(Domain domain){
		double[] time = new double[20];
		int count = 0;
		while(count < time.length){
			double time1 = System.nanoTime();
			domain.runProblem();
			double time2 = System.nanoTime();
			time[count] = (time2 - time1)/1.E9;
			count++;
		}
		System.out.println("Done. Finished " + count + " runs. Average run time was " + Arrays.stream(time).average().getAsDouble());
		System.out.println("Test results");
		for(int i = 0; i < count; i++){
			System.out.println(i + ":" + time[i]);
		}
	}
	
	static void run(Domain domain){
		double time1 = System.nanoTime();
		//domain.runProblemReflect();
		domain.runProblem();
		double time2 = System.nanoTime();
		System.out.println("Done. Finished run in: " + String.valueOf((time2-time1)/1E9) + " seconds");
	}
}
