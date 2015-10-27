package javaTrans;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;

import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartUtilities;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.plot.PolarPlot;
import org.jfree.chart.renderer.DefaultPolarItemRenderer;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;

public class PlotTools {
	static void printMeshPlots(ArrayList<Region> regions, String name) throws IOException{
		XYSeriesCollection total = new XYSeriesCollection();
		for(int e = 0, eng = regions.get(0).conts.eBins; e < eng; e++){
			XYSeries temp = new XYSeries(String.valueOf(e));
			for(Region reg: regions){
				for(Mesh mesh: reg.meshPoints){
					temp.add(mesh.xPosition, mesh.energyFlux[e]);
				}
			}
			total.addSeries(temp);
		}
		XYSeries temp = new XYSeries(String.valueOf("Total"));
		for(Region reg: regions){
			for(Mesh mesh: reg.meshPoints){
				temp.add(mesh.xPosition, mesh.finalFlux);
			}
		}
		total.addSeries(temp);
		
		JFreeChart xylineChart = ChartFactory.createXYLineChart(
				"Fluxes", 
				"X",
				"Flux", 
				total,
				PlotOrientation.VERTICAL, 
				true, true, false);
		File XYChart = new File(name+".jpeg"); 
		ChartUtilities.saveChartAsJPEG( XYChart, xylineChart, 2560, 1920);
	}
	
	static double maxFluxCalc(Mesh mesh, int edge){
		double maxFlux = 0;
		for(int m = 0, mew = mesh.region.conts.mew.length; m < mew; m++){
			for(int e = 0, eng = mesh.region.conts.eBins; e < eng; e++){
				if(mesh.fluxTArray[edge][0][m][e] > maxFlux){maxFlux = mesh.fluxTArray[edge][0][m][e];}
			}
		}
		return maxFlux;
	}
	
	static void printPolarCenter(Mesh mesh, int energy, int edge, String name){
		double maxFlux = maxFluxCalc(mesh, edge);
		XYSeries temp = new XYSeries(energy);
		for(int m = 0, mew = mesh.region.conts.mew.length; m < mew; m++){
			temp.add(Math.acos(mesh.region.conts.mew[m])*(360/(Math.PI*2)), mesh.fluxTArray[edge][0][m][energy]/maxFlux);
		}
		XYSeriesCollection data = new XYSeriesCollection();
		data.addSeries(temp);
		JFreeChart polar = ChartFactory.createPolarChart(
				"Ordinate Fluxes",
				data, 
				true, 
				true,
				false);
		PolarPlot plot = (PolarPlot) polar.getPlot();
		DefaultPolarItemRenderer render = (DefaultPolarItemRenderer) plot.getRenderer();
		render.setSeriesFilled(2, true);
		File XYChart = new File(name+".jpeg"); 
		try {
			ChartUtilities.saveChartAsJPEG(XYChart, polar, 2560, 1920);
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
	static void printPolarRight(Mesh mesh, int energy, int edge, String name){
		double maxFlux = maxFluxCalc(mesh, 2);
		XYSeries temp = new XYSeries(energy);
		for(int m = 0, mew = mesh.region.conts.mew.length; m < mew; m++){
			temp.add(Math.acos(mesh.region.conts.mew[m])*(360/(Math.PI*2)), mesh.fluxTArray[edge][0][m][energy]/maxFlux);
		}
		XYSeriesCollection data = new XYSeriesCollection();
		data.addSeries(temp);
		JFreeChart polar = ChartFactory.createPolarChart(
				"Ordinate Fluxes",
				data, 
				true, 
				true,
				false);
		PolarPlot plot = (PolarPlot) polar.getPlot();
		DefaultPolarItemRenderer render = (DefaultPolarItemRenderer) plot.getRenderer();
		render.setSeriesFilled(2, true);
		File XYChart = new File(name+".jpeg"); 
		try {
			ChartUtilities.saveChartAsJPEG(XYChart, polar, 2560, 1920);
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
	static void printPolarLeft(Mesh mesh, int energy, int edge, String name){
		double maxFlux = maxFluxCalc(mesh, 1);
		XYSeries temp = new XYSeries(energy);
		for(int m = 0, mew = mesh.region.conts.mew.length; m < mew; m++){
			temp.add(Math.acos(mesh.region.conts.mew[m])*(360/(Math.PI*2)), mesh.fluxTArray[edge][0][m][energy]/maxFlux);
		}
		XYSeriesCollection data = new XYSeriesCollection();
		data.addSeries(temp);
		JFreeChart polar = ChartFactory.createPolarChart(
				"Ordinate Fluxes",
				data, 
				true, 
				true,
				false);
		PolarPlot plot = (PolarPlot) polar.getPlot();
		DefaultPolarItemRenderer render = (DefaultPolarItemRenderer) plot.getRenderer();
		render.setSeriesFilled(2, true);
		File XYChart = new File(name+".jpeg"); 
		try {
			ChartUtilities.saveChartAsJPEG(XYChart, polar, 2560, 1920);
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
}
