/**
 * \package test
 * \brief Package of classes used to test methods within iDynoMiCS
 * 
 * Package of classes used to test methods within iDynoMiCS. This package is part of iDynoMiCS v1.2, governed by the 
 * CeCILL license under French law and abides by the rules of distribution of free software.  You can use, modify and/ or redistribute 
 * iDynoMiCS under the terms of the CeCILL license as circulated by CEA, CNRS and INRIA at the following URL  "http://www.cecill.info".
 */
package test;

import utils.GridOperations;

public class TestUtils 
{

	/**
	 * \brief Main program script used to call features of iDynoMiCS being tested
	 * 
	 * Main program script used to call features of iDynoMiCS being tested
	 * 
	 * @param args	Any system input arguments required by this test
	 */
	public static void main(String[] args) 
	{
		int[] a = {64,64,64};
		int[] b = {16,16,16};
		int[] c = {32,16,16};
		int[] d = {10,10,10};

		boolean ab = GridOperations.isConsistent(a,b);
		boolean ac = GridOperations.isConsistent(a,c);
		boolean ad = GridOperations.isConsistent(a,d);
		
		System.out.println(ab);
		System.out.println(ac);
		System.out.println(ad);
		System.out.println("End");
	}

}
