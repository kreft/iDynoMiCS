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

import java.awt.BorderLayout;
import java.awt.FlowLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.net.MalformedURLException;
import java.net.URL;
import java.net.URLConnection;
import java.util.Enumeration;
import java.util.zip.ZipEntry;
import java.util.zip.ZipFile;
import javax.swing.JButton;
import javax.swing.JFrame;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JTextArea;

/**
 * \brief Performs an update to the iDynomics tool if a new version has been detected
 * 
 * This class should be thought of as a standalone class. It may be packaged within iDynomics, and called by iDynomics, but is run 
 * independently as its own jar file. The purpose of this is to download the latest version from the web repository, extract the files, 
 * update the currently downloaded version to the new version, and then delete files used in the update. As this is changing class files, 
 * one of which could be the main simulator, the simulation should be shut down when this is done. Therefore, iDynomics calls this jar and 
 * exits - the update can be performed, and the user should then restart iDynomics.
 * 
 * This code has been amended from a freely available example update procedure written by Thomas Otero, documented here:
 * http://www.dreamincode.net/forums/topic/190944-creating-an-updater-in-java/ 
 * 
 * @author Thomas Otero - Created the original code and made this available online. This has been modified from his code.
 * @author Kieran Alden - Altered this code to fit iDynomics specifically and added the required documentation to suit doxygen output
 */
public class Update_IDynomics extends JFrame 
{
	/**
	 * Path to the web address where url.html can be found - this file contains the web address of the update package
	 */
	public String pathToURLhtml = "http://www.biosciences-labs.bham.ac.uk/kreftlab/iDynomics/updater_check/url.html";

	/**
	 * Thread that performs the update - used in the download method
	 */
	private Thread worker;
	
	/**
	 * Path to the temporary update folder that will be created within the simulation folder during the process
	 */
	private final String root = "update/";

	/**
	 * Text area that will be populated with relevant text about the procedure, to tell the user whats going on
	 */
	private JTextArea outText;
	
	/**
	 * Button to cancel the update, and to exit the process at the end
	 */
	private JButton cancel;
	
	/**
	 * Scroll panel that will be used by the status screen
	 */
	private JScrollPane sp;
	
	/**
	 * Panel used by status frame - specifically to contain the text area
	 */
	private JPanel pan1;
	
	/**
	 * Panel used by status frame - specifically to contain the buttons
	 */
	private JPanel pan2;

	/**
	 * \brief Constructor to launch the procedure of an iDynomics update
	 * 
	 * This constructor begins the update procedure by calling the methods to create a status window that specifies update progress, 
	 * and the methods to perform the download of the update and copying of the relevant files
	 */
	public Update_IDynomics() 
	{
		// Create the status window 
		initComponents();
		
		// Tell the user what is happening
		outText.setText("Contacting Download Server...");
		
		// Perform the update
		download();
	}

	/**
	 * \brief Creates a status window to show the progress of the update procedure
	 * 
	 * This method creates a Java Swing window which is used to give the user updates on the progress of the iDynomics update. This 
	 * panel can be updated via the other methods
	 */
	private void initComponents() 
	{
		// Set the action that should occur if this window is shut
		setDefaultCloseOperation(javax.swing.WindowConstants.EXIT_ON_CLOSE);

		// Create some panels - one for the text that is updated and another for the cancel/exit button
		pan1 = new JPanel();
		pan1.setLayout(new BorderLayout());

		pan2 = new JPanel();
		pan2.setLayout(new FlowLayout());

		// Create a text area for the update progress lines
		outText = new JTextArea();
		sp = new JScrollPane();
		sp.setViewportView(outText);

		// Create a cancel button. This is altered to 'Exit' at the end of the procedure
		cancel = new JButton("Cancel Update");
		cancel.addActionListener(new ActionListener() 
		{
			public void actionPerformed(ActionEvent e) {
				System.exit(0);
			}
		});

		// Add the button to the panel
		pan2.add(cancel);
		pan1.add(sp, BorderLayout.CENTER);
		pan1.add(pan2, BorderLayout.SOUTH);

		// Add the panel to the frame
		add(pan1);
		pack();
		this.setSize(500, 400);
	}

	/**
	 * \brief Performs the update of the iDynomics simulation tool, calling relevant methods to do so
	 * 
	 * This method is called to perform the update of iDynomics, calling the required methods to download the update, unzip it, 
	 * copy the files over, and delete the temporary update files. The user is updated on the progress in a status window as the update 
	 * proceeds
	 */
	private void download() 
	{
		worker = new Thread(new Runnable() {
			public void run() {
				try 
				{
					// Download the update file
					downloadFile(getDownloadLinkFromHost());
					
					// Unzip it
					unzip();
					
					// Copy the files over, replacing former versions or adding new files
					copyFiles(new File(root), new File("").getAbsolutePath());
					
					// Delete the update files
					cleanup();
					
					// Tell the user the update is finished, and change the text on the cancel button
					outText.setText(outText.getText() + "\nUpdate Finished!");
					cancel.setText("Exit iDynoMiCS Updater");
				} 
				catch (Exception ex) 
				{
					ex.printStackTrace();
					JOptionPane.showMessageDialog(null,
							"An error occured while preforming update!");
				}
			}
		});
		worker.start();
	}

	/**
	 * \brief Removes the update files downloaded from the server.
	 * 
	 * This method removes the update.zip file and the extracted update folder from the simulation folder once the update has been 
	 * successful.
	 */
	private void cleanup() 
	{
		// Tell the user what is happening
		outText.setText(outText.getText() + "\nPreforming clean up...");
		
		// Remove the zip file
		File f = new File("update.zip");
		f.delete();
		
		// Remove the extracted files by calling the remove method
		remove(new File(root));
		new File(root).delete();
	}

	/**
	 * \brief Part of the Cleanup procedure - removes files that were extracted from the update.zip package
	 * 
	 * This method is part of the cleanup procedure, and removes all files within the 'update' folder created when the 
	 * update zip file was extracted
	 * 
	 * @param updateDirectory	Path to the directory to be removed
	 */
	private void remove(File updateDirectory) {
		File[] files = updateDirectory.listFiles();
		for (File ff : files) {
			if (ff.isDirectory()) {
				remove(ff);
				ff.delete();
			} else {
				ff.delete();
			}
		}
	}

	/**
	 * \brief Copy updated files to iDynomics simulation repository
	 * 
	 * This method examines all the files in the update and replaces the current files in the simulation directory with these files.
	 * Each directory is examined in turn, and this method called recursively to update all the files in that directory.
	 * 
	 * @param dirFiles	The name of the update directory being examined
	 * @param dir	The absolute path to the simulation directory
	 * @throws IOException	Exception thrown should there be a problem copying these files
	 */
	private void copyFiles(File dirFiles, String dir) throws IOException 
	{
		/**
		 * An array of the file names in this update directory
		 */
		File[] files = dirFiles.listFiles();
		
		// Iterate through each of these files
		for (File ff : files) 
		{
			// Check if this is an internal directory
			if (ff.isDirectory()) 
			{
				// Make that directory if it currently does not exist in the simulator
				new File(dir + "/" + ff.getName()).mkdir();
				// Call this method again with all the files in this directory
				copyFiles(ff, dir + "/" + ff.getName());
			}
			else 
			{
				// Copy the updated files in place of the old ones.
				copy(ff.getAbsolutePath(), dir + "/" + ff.getName());
			}

		}
		
	}

	/**
	 * \brief Replaces old version of a file with the updated version from the repository
	 * 
	 * This method replaces an old version of a simulation file with an up to date version
	 * 
	 * @param srFile	Path to the new updated file
	 * @param dtFile	Path to the file to be replaced
	 * @throws FileNotFoundException	Exception should this file not be found	
	 * @throws IOException	Exception should there be an input/output error during copying
	 */
	public void copy(String srFile, String dtFile) throws FileNotFoundException, IOException 
	{
		// Create the two file handlers
		File f1 = new File(srFile);
		File f2 = new File(dtFile);

		// Read in the updated file
		InputStream in = new FileInputStream(f1);

		// Create an output stream for the updated file
		OutputStream out = new FileOutputStream(f2);

		// Perform the update
		byte[] buf = new byte[1024];
		int len;
		while ((len = in.read(buf)) > 0) {
			out.write(buf, 0, len);
		}
		
		// Close the buffers
		in.close();
		out.close();
	}

	/**
	 * \brief Unzips the downloaded update package to enable the required files to be updated
	 * 
	 * This method unzips the update.zip package that was downloaded from the server. This is unzipped into a folder called update.
	 * Once the update has been completed, this folder is removed.
	 * 
	 * @throws IOException	Exception generated if the zip file cannot be found or the files cannot be extracted onto the hard disk
	 */
	private void unzip() throws IOException 
	{
		// Set up the streams required to perform the unzip
		int BUFFER = 2048;
		BufferedOutputStream dest = null;
		BufferedInputStream is = null;
		ZipEntry entry;
		
		// Name of the zip file downloaded - this is always update.zip
		ZipFile zipfile = new ZipFile("update.zip");
		Enumeration e = zipfile.entries();
		(new File(root)).mkdir();
		
		// Perform the unzip
		while (e.hasMoreElements()) {
			entry = (ZipEntry) e.nextElement();
			outText.setText(outText.getText() + "\nExtracting: " + entry);
			if (entry.isDirectory())
				(new File(root + entry.getName())).mkdir();
			else {
				(new File(root + entry.getName())).createNewFile();
				is = new BufferedInputStream(zipfile.getInputStream(entry));
				int count;
				byte data[] = new byte[BUFFER];
				FileOutputStream fos = new FileOutputStream(root
						+ entry.getName());
				dest = new BufferedOutputStream(fos, BUFFER);
				while ((count = is.read(data, 0, BUFFER)) != -1) {
					dest.write(data, 0, count);
				}

				dest.flush();
				dest.close();
				is.close();
			}
		}
		zipfile.close();

	}

	/**
	 * \brief Downloads the simulation update from the repository
	 * 
	 * This method downloads the latest release of the simulation from the iDynomics repository, storing this as update.zip.
	 * This references an address in the file url.html, that should be contained within the iDynomics repository
	 * 
	 * @param link	The link specified in url.html - fetched by the method that calls this one
	 * @throws MalformedURLException	Exception generated should this link be wrong
	 * @throws IOException	Exception generated should the file not be able to be stored (i/o error on storage)
	 * 
	 */
	private void downloadFile(String link) throws MalformedURLException,IOException 
	{
		// Store the link to the update file as a URL
		//URL url = new URL(link);
		URL url = new URL("http://www.kieranalden.info/update/iDynoMiCSv1_3.zip");
		
		// Open this connection
		URLConnection conn = url.openConnection();
		InputStream is = conn.getInputStream();
		long max = conn.getContentLength();
		
		// Set some status text to show that an update has been found and something is being done
		outText.setText(outText.getText() + "\n"
				+ "Downloding file...\nUpdate Size(compressed): " + max
				+ " Bytes");
		
		// Save the downloaded file as update.zip
		BufferedOutputStream fOut = new BufferedOutputStream(
				new FileOutputStream(new File("update.zip")));
		
		
		byte[] b = new byte[1024];
        int count;

        // Download the file
        while ((count = is.read(b)) >= 0) {
            fOut.write(b, 0, count);
        }
        // Close the buffers
        fOut.close();
        is.close();
		fOut.flush();
		fOut.close();
		is.close();
		
		// Inform the user that the download is complete
		outText.setText(outText.getText() + "\nDownload Complete!");

	}

	/**
	 * \brief Retrieves the address of the update from the iDynomics repository
	 * 
	 * The updater relies on a file within the iDynomics repository - url.html - that contains the full web address to the 
	 * update package zip file. This address is specified within [url][/url] tags. This method processes this file and retrieves the 
	 * URL from these tags
	 * 
	 * @return	A string containing the web address of the update package
	 * @throws MalformedURLException	Exception generated if the URL that specifies where url.html can be located is incorrect
	 * @throws IOException	Exception generated if a stream cannot be opened to read this file
	 */
	private String getDownloadLinkFromHost() throws MalformedURLException,IOException 
	{
		// Generate a URL from the string retrieved from URL.html
		URL url = new URL(pathToURLhtml);

		// Read this string as a buffer, deleting the [url] and [/url] tags, to return the address of the package
		InputStream html = null;

		html = url.openStream();

		int c = 0;
		StringBuilder buffer = new StringBuilder("");

		while (c != -1) {
			c = html.read();
			buffer.append((char) c);

		}
		
		// Return the web address of the update package
		return buffer.substring(buffer.indexOf("[url]") + 5,
				buffer.indexOf("[/url]"));
	}

	/**
	 * \brief Main method - called to begin the process of updating iDynomics
	 * 
	 * This method is launched by the simulator when the user chooses to update iDynomics. The relevant download, unzip, and copy 
	 * methods are then run as required.
	 * 
	 * @param args	Program arguments - though none are required in this case
	 */
	public static void main(String args[]) {
		java.awt.EventQueue.invokeLater(new Runnable() {
			public void run() 
			{
				new Update_IDynomics().setVisible(true);
			}
		});
	}
}
