package javaTrans;

import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Vector;
import java.util.stream.Stream;

public class IsoReader {
	
	public static void readIsos(ArrayList<Isotope> isoArray){
		Stream<String> lines = null;
		try {
			lines = Files.lines(Paths.get("Isos/Isos.txt"), StandardCharsets.UTF_8);
		} catch (IOException e) {
			e.printStackTrace();
		}
		Vector<String> strings = new Vector<String>();
		lines.forEach((line) -> strings.add(line));
	}
	
	public static void buildProblem(ArrayList<Isotope> isoArray, Vector<String> strings){
		
	}
}
