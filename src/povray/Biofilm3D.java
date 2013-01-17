/**
 * Project iDynoMiCS (copyright -> see Idynomics.java)
 * ______________________________________________________
 * 
 * Creates a union with a box (the biofilm carrier) and spheres (the bacteria)
 * Create a 3D scene for the PovRay rendering engine
 * 
 */

/**
 * @since Feb 2007
 * @version 1.0
 * @author João Xavier (xavierj@mskcc.org), Memorial Sloan-Kettering Cancer Center (NY, USA)
 * @author Laurent Lardon (lardonl@supagro.inra.fr), INRA, France
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


public class Biofilm3D implements Serializable {
	// Serial version used for the serialisation of the class
	private static final long serialVersionUID = 1L;

	private Povray3DScene         _pov;
	private VectorProperty        translate, rotate;
	private LinkedList<Box>        _biofilmCarrier = new LinkedList<Box>();
	private ParticleWithCapsule[] _cells;
	private int                   _next;

	/* ________________ CONSTRUCTOR _______________________ */
	protected Biofilm3D(Povray3DScene pov) {
		_pov = pov;
		initializeAll();
	}

	protected Biofilm3D(Povray3DScene pov, int n) {
		_pov = pov;
		initializeAll();
		_cells = new ParticleWithCapsule[n];
	}

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
	 * Add a new cell
	 * 
	 * @param x
	 * @param y
	 * @param z
	 * @param rad
	 * @param r
	 * @param g
	 * @param b
	 */
	protected void addCell(float x, float y, float z, float rad, int r, int g, int b) {
		Color color = new Color(r, g, b);
		_cells[_next] = new ParticleWithCapsule();
		_cells[_next].setCenter(new ContinuousVector(x, y, z));
		_cells[_next].setCoreRadius(rad);
		_cells[_next].setColorCore(color);
		_next++;
	}

	/**
	 * Writes the cell positions to pov file
	 * 
	 * @param f
	 * @throws IOException
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
	 * Writes the current cells position in model to file not using the _cells
	 * array, wich is too memory consuming
	 * 
	 * @param f
	 * @throws IOException
	 */
	protected void modelStateToFile(FileWriter f) throws IOException {
		biofilmHeaderToFile(f);
		particlesToFile(f);
		biofilmFooterToFile(f);
	}

	/**
	 * Writes the union open and the carrier to file
	 * 
	 * @param f
	 * @throws IOException
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
	 * writes the union tranlations and rotations and close to file
	 * 
	 * @param f
	 * @throws IOException
	 */
	protected void biofilmFooterToFile(FileWriter f) throws IOException {
		f.write("\t"+translate+"\n");
		f.write("\t"+rotate+"\n");
		f.write("}");
	}

	/**
	 * Write the particles only to file
	 * 
	 * @param f
	 * @throws IOException
	 */
	protected void particlesToFile(FileWriter f) throws IOException {
		for (Agent anAgent : _pov.mySim.agentGrid.agentList) {
			ParticleWithCapsule s = new ParticleWithCapsule((LocatedAgent)anAgent);
			f.write(s.toString());
		}
	}

}
