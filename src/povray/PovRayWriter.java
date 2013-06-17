/**
 * \package povray
 * \brief Package of classes used to create output files that can be processed using POV-Ray software.
 * 
 * Package of classes used to create output files that can be processed using POV-Ray software. This package is part of iDynoMiCS v1.2, 
 * governed by the CeCILL license under French law and abides by the rules of distribution of free software.  You can use, modify and/ 
 * or redistribute iDynoMiCS under the terms of the CeCILL license as circulated by CEA, CNRS and INRIA at the following URL  
 * "http://www.cecill.info".
 */
package povray;

import de.schlichtherle.io.File;

// import java.io.File;
import java.io.IOException;
import java.io.Serializable;

import utils.LogFile;

import simulator.Simulator;

/**
 * \brief Class that writes POV-Ray files that can be used to visualise simulation results
 * 
 * Class that writes POV-Ray files that can be used to visualise simulation results
 * 
 * @author João Xavier (xavierj@mskcc.org), Memorial Sloan-Kettering Cancer Center (NY, USA)
 * @author Laurent Lardon (lardonl@supagro.inra.fr), INRA, France
 */
public class PovRayWriter implements Serializable 
{

	/**
	 * Serial version used for the serialisation of the class
	 */
	private static final long serialVersionUID = 1L;

	/**
	 * A 3D scene for the PovRay rendering engine
	 */
	private Povray3DScene     _povRay;
	
	/**
	 * Directory in which the POV-Ray output files should be stored
	 */
	private String            dir;

	/**
	 * Archive (.zip) file in which all POV-Ray output files will be stored
	 */
	private File              _vArchive;

	/**
	 * \brief Initialises a POV-Ray writer object to produce simulation statistics that can be presented using POV-Ray
	 * 
	 * Initialises a POV-Ray writer object to produce simulation statistics that can be presented using POV-Ray. Creates the required 
	 * output file paths and the archive in which these files will be stored
	 * 
	 * @param aSim	The simulation object used to simulate the conditions specified in the protocol file
	 * @param outPath	A string stating the directory in which POV-Ray output files should be stored
	 */
	public void initPovRay(Simulator aSim, String outPath) 
	{
		try 
		{
			dir = outPath+File.separator;
			// Create the PovRay description object
			_povRay = new Povray3DScene(aSim, aSim.world.domainList.get(0).getName());
			_povRay.writePovrayIncFiles(dir+"lastIter"+File.separator);			

			// Create the archive file in which these
			_vArchive = new File(dir+"povray"+".zip");
			File incFile=new File(dir+"lastIter"+File.separator+"sceneheader.inc");
			incFile.copyTo(new File(dir+"povray.zip"+File.separator+"sceneheader.inc"));

			incFile=new File(dir+"lastIter"+File.separator+"scenefooter.inc");
			incFile.copyTo(new File(dir+"povray.zip"+File.separator+"scenefooter.inc"));
			File.update(_vArchive);

		} 
		catch (Exception e) 
		{
			LogFile.writeError(e.getLocalizedMessage(), "PovRayWriter.initPovRay()");
		}
	}

	/**
	 * \brief Writes model state information to the POV-Ray file for representation in POV-Ray
	 * 
	 * Writes model state information to the POV-Ray file for representation in POV-Ray. Utilises the current run iteration within the 
	 * file name to state when this snapshot was taken (added by BVM 27.1.09)
	 * 
	 * @param fileIndex	The current run iteration. Used in the creation of the file name
	 */
	public void write(int fileIndex) 
	{
		try 
		{
			// Create the povray file
			File f = new File(_povRay.writeModelState(dir+"lastIter"+File.separator+"it(last).pov"));

			// Copy the povray file inside the archive
			f.copyTo(new File(dir+"povray.zip"+File.separator+"it("+fileIndex+").pov"));
			File.update(_vArchive);

		} 
		catch (IOException e) 
		{
			System.out.println("Error trying to write povRayFile");
		}
	}

}