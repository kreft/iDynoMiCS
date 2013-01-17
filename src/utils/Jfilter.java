/**
 * Project iDynoMiCS (copyright -> see Idynomics.java)
 */


package utils;
import java.io.*;

public class Jfilter implements FilenameFilter {

	String _extension;

	public Jfilter(String extension) {
		_extension = "."+extension;
	}

	public boolean accept(File directory, String filename) {

		if (filename.endsWith(_extension)) return true;
		return false;

	}
}
