/**
 * Project iDynoMiCS (copyright -> see Idynomics.java) 
 */

/**
 * 
 * @since March, 2004 
 * @version 1.0
 * @author  * @author João Xavier (xavierj@mskcc.org), Memorial Sloan-Kettering Cancer Center (NY, USA)
 * 
 */

package utils;

import java.io.*;
import java.util.zip.*;

public class ZipArchive {

	static final int        BUFFER = 2048;
	private ZipOutputStream out;
	private String          _zipFileName;

	/**
	 * 
	 */
	public ZipArchive(String zipFileName) throws IOException {
		// create file
		_zipFileName = zipFileName;
		out = new ZipOutputStream(new BufferedOutputStream(new FileOutputStream(_zipFileName)));
	}

	/**
	 * Add a file to a zip archive
	 * @param fn
	 * @param zipArchive
	 */
	public final void addToZipArchiveAndDelete(File f) throws IOException {
		int count;
		byte data[] = new byte[BUFFER];

		BufferedInputStream origin = new BufferedInputStream(new FileInputStream(f), BUFFER);
		//out = new ZipOutputStream(new BufferedOutputStream(new FileOutputStream(_zipFileName)));
		out.putNextEntry(new ZipEntry(f.getName()));

		while ((count = origin.read(data, 0, BUFFER))!=-1) {
			out.write(data, 0, count);
		}
		out.closeEntry();
		origin.close();
		f.delete();
	}

	public final void addToZipArchiveAndKeepLast(File f, File fNew) throws IOException {
		int count;
		byte data[] = new byte[BUFFER];

		BufferedInputStream origin = new BufferedInputStream(new FileInputStream(f), BUFFER);
		//out = new ZipOutputStream(new BufferedOutputStream(new FileOutputStream(_zipFileName)));
		out.putNextEntry(new ZipEntry(f.getName()));

		while ((count = origin.read(data, 0, BUFFER))!=-1) {
			out.write(data, 0, count);
		}
		origin.close();
		//out.close();
		f.renameTo(fNew);
	}

	/**
	 * Closes the zip archive
	 * @throws IOException
	 */
	public void close() throws IOException {
		out.close();
	}
}
