/**
 * \package utils
 * \brief Package of classes that perform utility functions in the process of running an iDynoMiCS Simulation
 * 
 * Package of classes that perform utility functions in the process of running an iDynoMiCS Simulation. This package is part of iDynoMiCS v1.2, governed by the 
 * CeCILL license under French law and abides by the rules of distribution of free software.  You can use, modify and/ or redistribute 
 * iDynoMiCS under the terms of the CeCILL license as circulated by CEA, CNRS and INRIA at the following URL  "http://www.cecill.info".
 */
package utils;
import java.awt.Color;

/**
 * \brief Static class used to convert units (mass, length, time)
 * 
 * Static class used to convert units (mass, length, time). Internal units: Mass = femtograms (=1e-15 grams), Length = micrometer (=1e-6 meter), 
 * Time = hour
 * 
 * @author Andreas Dötsch (andreas.doetsch@helmholtz-hzi.de), Helmholtz Centre for Infection Research (Germany)
 * @author Laurent Lardon (lardonl@supagro.inra.fr), INRA, France
 */
public class UnitConverter 
{

	/**
     * \brief Takes a unit of time and returns a factor to multiply a parameter by to obtain the correct unit (hour)
     * 
     * Takes a unit of time and returns a factor to multiply a parameter by to obtain the correct unit (hour)
     * 
     * @param timeUnit	Unit of time for which a conversion is required
     * @return a factor to multiply to obtain the correct unit (hour)
     */
	public static Double time(String timeUnit) {
		Double out = 1.0;
		String unit = "";

		if (timeUnit.contains("day")) {
			out = 24.0;
			unit = "day";
		}
		if (timeUnit.contains("hour")) {
			out = 1.0;
			unit = "hour";
		}
		if (timeUnit.contains("minute")) {
			out = 1.0/60.0;
			unit = "minute";
		}
		if (timeUnit.contains("second")) {
			out = 1.0/3600.0;
			unit = "second";
		}

		if (timeUnit.contains(unit+"-1")) {
			out = 1.0/out;
		}

		return out;
	}

	/**
	 * \brief Takes a length unit and returns a double value of length for that unit
	 * 
	 * Takes a length unit and returns a length value of mass for that unit
	 * 
	 * @param lengthUnit	String containing a unit of length (m,cm,etc)
	 * @return	Double value of length for that unit
	 */
	public static Double length(String lengthUnit)
	{
		Double out = 1.0;
		String unit = "";

		if (lengthUnit.contains("m")) {
			out = 1e6;
			unit = "m";
		}
		if (lengthUnit.contains("cm")) {
			out = 1e4;
			unit = "cm";
		}
		if (lengthUnit.contains("mm")) {
			out = 1e3;
			unit = "mm";
		}
		if (lengthUnit.contains("μm")) {
			out = 1.0;
			unit = "microm";
		}
		if (lengthUnit.contains("um")) {
			out = 1.0;
			unit = "microm";
		}
		if (lengthUnit.contains("fm")) {
			out = 1e-9;
			unit = "fm";
		}

		if (lengthUnit.contains(unit+"-")) {
			out = 1/out;
		}

		return out;
	}

	/**
	 * \brief Takes a mass unit and returns a double value of mass for that unit
	 * 
	 * Takes a mass unit and returns a double value of mass for that unit
	 * 
	 * @param massUnit	String containing a unit of mass (g,kg,etc)
	 * @return	Double value of mass for that unit
	 */
	public static double mass(String massUnit) {
		double out = 1;
		String unit = "";

		if (massUnit.contains("g")) {
			out = 1e15;
			unit = "g";
		}

		if (massUnit.contains("kg")) {
			out = 1e18;
			unit = "kg";
		}

		if (massUnit.contains("mg")) {
			out = 1e12;
			unit = "mg";
		}

		if (massUnit.contains("μg")) {
			out = 1e9;
			unit = "μg";
		}
		if (massUnit.contains("ug")) {
			out = 1e9;
			unit = "μg";
		}
		if (massUnit.contains("fg")) {
			out = 1;
			unit = "fg";
		}

		if (massUnit.contains(unit+"-1")) {
			out = 1/out;
		}
		return out;
	}

	/**
	 * \brief Takes a mass unit and returns a double value of volume for that unit
	 * 
	 * Takes a mass unit and returns a double value of volume for that unit
	 * 
	 * @param massUnit	String containing a unit of mass (m,L,etc)
	 * @return	Double value of mass for that unit
	 */
	public static double volume(String massUnit) 
	{
		double out = 1;
		String unit = "";

		if (massUnit.contains("m")) {
			out = 1e18;
			unit = "m";
		}
		if (massUnit.contains("L")) {
			out = 1e15;
			unit = "L";
		}

		if (massUnit.contains("μm")) {
			out = 1;
			unit = "μm";
		}
		if (massUnit.contains("um")) {
			out = 1;
			unit = "μm";
		}

		if (massUnit.contains(unit+"-")) {
			out = 1/out;
		}
		return out;
	}

	/**
	 * \brief For a given colour name, returns a colour object
	 * 
	 * For a given colour name, returns a colour object
	 * 
	 * @param colorName	String name of a colour, for which a colour object is required
	 * @return	A Color object for that specified colour name
	 */
	public static Color getColor(String colorName) 
	{
		if (colorName.equals("red")) return Color.red;
		if (colorName.equals("blue")) return Color.blue;
		if (colorName.equals("cyan")) return Color.cyan;
		if (colorName.equals("darkgray")) return Color.darkGray;
		if (colorName.equals("lightgray")) return Color.lightGray;
		if (colorName.equals("gray")) return Color.gray;
		if (colorName.equals("green")) return Color.green;
		if (colorName.equals("magenta")) return Color.magenta;
		if (colorName.equals("orange")) return Color.orange;
		if (colorName.equals("pink")) return Color.pink;
		if (colorName.equals("yellow")) return Color.yellow;
		if (colorName.equals("white")) return Color.white;
		if (colorName.equals("black")) return Color.black;

		utils.LogFile.writeLog(colorName+" is not recognized as a color");
		return Color.white;
	}
}
