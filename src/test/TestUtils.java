/**
 * Project iDynoMiCS (copyright -> see Idynomics.java)
 */


package test;

import utils.GridOperations;


public class TestUtils {


	public static void main(String[] args) {
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
