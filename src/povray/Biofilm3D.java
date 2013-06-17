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

import java.awt.Color;
import java.io.FileWriter;
import java.io.IOException;
import java.io.Serializable;

import java.util.*;

import simulator.geometry.ContinuousVector;
import simulator.geometry.boundaryConditions.*;
import simulator.agent.*;

/**
 * \brief Creates a union with a box (the biofilm carrier) and spheres (the bacteria), creating a 3D scene for the PovRay rendering engine
 * 
 * Creates a union with a box (the biofilm carrier) and spheres (the bacteria), creating a 3D scene for the PovRay rendering engine
 * 
 * @author João Xavier (xavierj@mskcc.org), Memorial Sloan-Kettering Cancer Center (NY, USA)
 * @author Laurent Lardon (lardonl@supagro.inra.fr), INRA, France
 */
public class Biofilm3D implements Serializable 
{
	/**
	 *  Serial version used for the serialisation of the class
	 */
	private static final long serialVersionUID = 1L;

	/**
	 * A POV-Ray 3D scene in which the display will be created
	 */
	private Povray3DScene         _pov;
	
	/**
	 * Translation used to set the agent coordinates to the correct output angle
	 */
	private VectorProperty        translate;
	
	/**
	 * Rotation used to set the agent coordinates to the correct output angle
	 */
	private VectorProperty		rotate;
	
	/**
	 * Linked list of box objects that will be drawn in the output
	 */
	private LinkedList<Box>        _biofilmCarrier = new LinkedList<Box>();
	
	/**
	 * Array of all cells to be displayed, stored as ParticleWithCapsule objects
	 */
	private ParticleWithCapsule[] _cells;
	
	/**
	 * Reference to the next cell to be processed in the _cells array
	 */
	private int                   _next;

	/**
	 * \brief Create a Biofilm3D scene from a specified POV-Ray 3D scene and initialise the output
	 * 
	 * Create a Biofilm3D scene from a specified POV-Ray 3D scene and initialise the output
	 * 
	 * @param pov	A specific POV-Ray 3D scene to generate
	 */
	protected Biofilm3D(Povray3DScene pov) {
		_pov = pov;
		initializeAll();
	}

	/**
	 * \brief Create a Biofilm3D scene from a specified POV-Ray 3D scene and initialise the output and the array of cells to include
	 * 
	 * Create a Biofilm3D scene from a specified POV-Ray 3D scene and initialise the output and the array of cells to include
	 * 
	 * @param pov	A specific POV-Ray 3D scene to generate
	 * @param n	The number of cells that will be included in this scene
	 */
	protected Biofilm3D(Povray3DScene pov, int n) {
		_pov = pov;
		initializeAll();
		_cells = new ParticleWithCapsule[n];
	}

	/**
	 * \brief Initialise the POV-Ray scene, setting up required angles of translation and rotation if 3D, and establishing each of the domain boundaries in the scene
	 * 
	 * Initialise the POV-Ray scene, setting up required angles of translation and rotation if 3D, and establishing each of the domain boundaries in the scene
	 */
	private void initializeAll() {
		if (_pov.getDomain().is3D()) {
			translate = new VectorProperty("translate");
			translate.setValues(-_pov.getX()*0.4f, -_pov.getY()*0.5f, -_pov.getZ()*0.5f);

			rotate = new VectorProperty("rotate");
			rotate.setValues(0, 0, 90);
		} else {
			// just want to center the view
			translate = new VectorProperty("translate");
			translate.setValues(-_pov.getX()*0.5f, -_pov.getY()*0.5f, 0);

			rotate = new VectorProperty("rotate");
			rotate.setValues(0, 0, 90);
		}
		// OLD METHOD where by default the carrier is on the bottom
		// by default the carrier is on the bottom
		//_biofilmCarrier = new Box();
		//_biofilmCarrier.setColor(0.2f, 0.2f, 0.2f);
		//_biofilmCarrier.setCorner1(0, 0, 0);
		// thickness of box
		//_biofilmCarrier.setCorner2(-_pov.getX()*.01f, _pov.getY(), _pov.getZ());

		// NEW METHOD that actually looks at the boundaries to put boxes in
		// but note that this WILL NOT draw the x0y or yNz supports (if there are any)
		// because that blocks the default povray field of view.
		double x1, x2, y1, y2, z1, z2;
		String bdryName;
		for (AllBC aBdry : _pov.mySim.world.domainList.get(0).getAllBoundaries()) {
			if (aBdry.isSupport()
					&& !aBdry.getSide().equals("x0y") && !aBdry.getSide().equals("yNz")) {
				// only draw a box if it is a support AND it won't block the view

				bdryName = aBdry.getSide();

				Box tempbox = new Box();
				tempbox.setColor(0.2f, 0.2f, 0.2f);

				// now set corners by checking which side the boundary is on

				if (bdryName.contains("x") || bdryName.contains("X")) {
					// the bdry is parallel to x, so x varies on bdry
					// box corners are extents in x direction
					x1 = 0;
					x2 = _pov.getX();
				} else {
					// bdry is orthogonal to x, so x will not vary
					// need to get box corner based on if this is 0 or N side
					if (bdryName.contains("0")) {
						x1 = 0;
						x2 = -_pov.getX() * 0.01f;
					} else {
						x1 = _pov.getX();
						x2 = _pov.getX() * 1.01f;
					}
				}
				if (bdryName.contains("y") || bdryName.contains("Y")) {
					// the bdry is parallel to y, so y varies on bdry
					// box corners are extents in y direction
					y1 = 0;
					y2 = _pov.getY();
				} else {
					// bdry is orthogonal to y, so y will not vary
					// need to get box corner based on if this is 0 or N side
					if (bdryName.contains("0")) {
						y1 = 0;
						y2 = -_pov.getY() * 0.01f;
					} else {
						y1 = _pov.getY();
						y2 = _pov.getY() * 1.01f;
					}
				}
				if (bdryName.contains("z") || bdryName.contains("Z")) {
					// the bdry is parallel to z, so z varies on bdry
					// box corners are extents in z direction
					z1 = 0;
					z2 = _pov.getZ();
				} else {
					// bdry is orthogonal to z, so z will not vary
					// need to get box corner based on if this is 0 or N side
					if (bdryName.contains("0")) {
						z1 = 0;
						z2 = -_pov.getZ() * 0.01f;
					} else {
						z1 = _pov.getZ();
						z2 = _pov.getZ() * 1.01f;
					}
				}

				tempbox.setCorner1(x1, y1, z1);
				tempbox.setCorner2(x2, y2, z2);

				// now add this box to the list to be drawn
				_biofilmCarrier.add(tempbox);
			}
		}
	}

	/**
	 * \brief Add a new cell to the POV-Ray 3D output
	 * 
	 * Add a new cell to the POV-Ray 3D output
	 * 
	 * @param x	X coordinate of the cell to add
	 * @param y	Y coordinate of the cell to add
	 * @param z	Z coordinate of the cell to add
	 * @param rad	Radius of the cell to add
	 * @param r	Amount of red colour to use in creating the colour of this cell
	 * @param g	Amount of green colour to use in creating the colour of this cell
	 * @param b	Amount of blue colour to use in creating the colour of this cell
	 */
	protected void addCell(float x, float y, float z, float rad, int r, int g, int b) 
	{
		Color color = new Color(r, g, b);
		_cells[_next] = new ParticleWithCapsule();
		_cells[_next].setCenter(new ContinuousVector(x, y, z));
		_cells[_next].setCoreRadius(rad);
		_cells[_next].setColorCore(color);
		_next++;
	}

	/**
	 * \brief Writes the cell positions to the POV-Ray file output stream
	 * 
	 * Writes the cell positions to the POV-Ray file output stream
	 * 
	 * @param f	POV-Ray output file stream to be written to
	 * @throws IOException	Exception that is thrown if there are problems with this output stream
	 */
	protected void toFile(FileWriter f) throws IOException {
		f.write("union {\n");
		f.write(_biofilmCarrier.toString());
		for (int i = 0; i<_cells.length; i++) {
			f.write(_cells[i].toString());
		}
		f.write("\t"+translate+"\n");
		f.write("\t"+rotate+"\n");
		f.write("}");
	}

	/**
	 * \brief Writes the current cells position in model to the POV-Ray file output stream by not using the _cells array, which is too memory consuming
	 * 
	 * Writes the current cells position in model to the POV-Ray file output stream by not using the _cells array, which is too memory consuming
	 * 
	 * @param f	POV-Ray output file stream to be written to
	 * @throws IOException	Exception that is thrown if there are problems with this output stream
	 */
	protected void modelStateToFile(FileWriter f) throws IOException {
		biofilmHeaderToFile(f);
		particlesToFile(f);
		biofilmFooterToFile(f);
	}

	/**
	 * \brief Writes the union open and the carrier to file
	 * 
	 * Writes the union open and the carrier to file
	 * 
	 * @param f	POV-Ray output file stream to be written to
	 * @throws IOException	Exception that is thrown if there are problems with this output stream
	 */
	protected void biofilmHeaderToFile(FileWriter f) throws IOException {
		f.write("union {\n");
		if (!_biofilmCarrier.isEmpty()) {
			for (Box aBox : _biofilmCarrier) {
				f.write(aBox.toString());
			}
		}
	}

	/**
	 * \brief Writes the union translations and rotations to the POV-Ray file stream and closes the file
	 * 
	 * Writes the union translations and rotations to the POV-Ray file stream and closes the file
	 * 
	 * @param f	POV-Ray output file stream to be written to
	 * @throws IOException	Exception that is thrown if there are problems with this output stream
	 */
	protected void biofilmFooterToFile(FileWriter f) throws IOException {
		f.write("\t"+translate+"\n");
		f.write("\t"+rotate+"\n");
		f.write("}");
	}

	/**
	 * \brief Write the particles only to the POV-Ray file output stream
	 * 
	 * Write the particles only to the POV-Ray file output stream
	 * 
	 * @param f	POV-Ray output file stream to be written to
	 * @throws IOException	Exception that is thrown if there are problems with this output stream
	 */
	protected void particlesToFile(FileWriter f) throws IOException {
		for (Agent anAgent : _pov.mySim.agentGrid.agentList) {
			ParticleWithCapsule s = new ParticleWithCapsule((LocatedAgent)anAgent);
			f.write(s.toString());
		}
	}

}
