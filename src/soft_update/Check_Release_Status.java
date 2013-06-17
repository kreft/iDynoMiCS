/**
 * \package soft_update
 * \brief Package of classes required to detect whether a new release of iDynomics is available, and to perform the update
 * 
 * This package contains classes that are used to detect whether a new release of iDynomics is available, and to perform a software 
 * update that downloads this new version and makes the required changes that ensures the user is running the latest version. This 
 * package is part of iDynoMiCS v1.2, governed by the CeCILL license under French law and abides by the rules of distribution of free 
 * software.  You can use, modify and/ or redistribute iDynoMiCS under the terms of the CeCILL license as circulated by CEA, CNRS and 
 * INRIA at the following URL  "http://www.cecill.info".
 */
package soft_update;

import java.io.InputStream;
import java.net.URL;
import java.net.UnknownHostException;

/**
 * \brief Checks if a new release of iDynomics is available
 * 
 * This class contains methods that check whether there is a new version of iDynomics available. This is done by checking the content 
 * of a file, version.html, that resides in the iDynomics repository. If there is an update available, a further file is downloaded, 
 * history.html, that contains a detailed description of the changes made in the update.
 * 
 * @author Thomas Otero - Created the original code and made this available online. This has been modified from his code.
 * @author Kieran Alden - amended for our use
 *
 */
public class Check_Release_Status 
{
	/**
	 * Link to the version.html page that contains the number of the latest release of iDynomics
	 */
	private final static String versionURL = "http://www.biosciences-labs.bham.ac.uk/kreftlab/iDynomics/updater_check/version.html";
	
	/**
	 * Link to the file containing the list of differences between the new and old versions
	 */
	private final static String versionChanges = "http://www.biosciences-labs.bham.ac.uk/kreftlab/iDynomics/updater_check/history.html";
	
	/**
	 * \brief Gets the content of url.html - the latest released version number
	 * 
	 * Returns the latest release version number by retrieving the content of the url.html file specified by the versionURL parameter.
	 * This is used to determine if an update is required
	 * 
	 * @return	A text string containing the content of url.html
	 * @throws Exception	Exception generated if this web address cannot be found
	 */
	public static String getLatestVersion() throws Exception
	{
		// Get the data from the web page
		String data = getData(versionURL);
		
		// Now, if the user is not connected to the internet, this will have returned "999"
		// Check whether this is the case before processing
		
		if(data.equals("999"))
		{
			// Return the same 999, thus the main simulation knows there is no internet and no new upgrade
			return "999";
		}
		else
		{
			// Return just the version number, removing the version tags
		    return data.substring(data.indexOf("[version]")+9,data.indexOf("[/version]"));
		}
		
		
	}
	
	/**
	 * \brief Returns the text detailing what is new about this version, for display in the user info screen
	 * 
	 * This method returns the content of the file history.html, specified by the parameter versionChanges. This is then displayed to 
	 * the user when they are asked whether they wish to perform an update or not
	 * 
	 * @return	A text string containing the content of history.html
	 * @throws Exception	Exception generated if this web address cannot be found
	 */
	public static String getWhatsNew() throws Exception
	{
		// Get the file from the server
		String data = getData(versionChanges);
		
		// Get just the changes, removing the tags
		return data.substring(data.indexOf("[history]")+9,data.indexOf("[/history]"));
	}
	
	/**
	 * \brief Retrieves the data either in the history or version html files. Used by getWhatsNew and getLatestVersion
	 * 
	 * Reads in the html files specified by the address parameter, and returns the text contained in those files as a string. Used by 
	 * the version checking method and method that gets the info on the differences between the software versions
	 * 
	 * @param address	Web address of the file to be processed
	 * @return	A string containing the text in that file
	 * @throws Exception	An exception generated if that file cannot be found or read
	 */
	private static String getData(String address) throws Exception
	{
		
		// Have to be careful here as the simulation may not be running on a machine that has internet
		// So lets catch that exception
		
		try
		{
			URL url = new URL(address);
			InputStream html = null;
				 
			html = url.openStream();
				         
			int c = 0;
			StringBuffer buffer = new StringBuffer("");
			 
			while(c != -1) 
			{
				c = html.read();
				             
				buffer.append((char)c);
			}
			
			return buffer.toString();
		}
		catch(UnknownHostException e)
		{
			// We're going to return a version number of 999 - this will tell the other routine that there is no internet connection
			// and iDynoMiCS simulation should start
			return "999";
		}
		
	}
	
	/**
	 * \brief Test stub for testing the methods that get the version number and new release details
	 * 
	 * Test stub for testing the methods that get the version number and new release details
	 * 
	 * @param args	Arguments passed to the test stub - none required in this case
	 */
	public static void main(String[] args) 
	{
		try 
		{
			// Test whether the online file is read correctly
			System.out.println(Check_Release_Status.getLatestVersion());
		}
		catch (Exception ex) 
		{
			ex.printStackTrace();
		}
	}

}
