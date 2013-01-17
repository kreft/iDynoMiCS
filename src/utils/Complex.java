/**
 * Project iDynoMiCS (copyright -> see Idynomics.java)
 * 
 * ______________________________________________________
 * Class for storing complex numbers as two double representing real and imaginary part, respectively.
 * Includes basic math operations (addition, subtraction, multiplication, division)
 */

/**
 * 
 * @since August 2006
 * @version 1.0
 * @author Andreas Dötsch (andreas.doetsch@helmholtz-hzi.de), Helmholtz Centre for Infection Research (Germany)
 */

package utils;


public class Complex {
	

	private double real, imag;
	
	/**
	 * Default constructor, setting number to 0 + 0i
	 *
	 */
	public Complex(){
		this.real = 0d;
		this.imag = 0d;
	}
	
	/**
	 * specific constructor, setting number to given real and imaginary part
	 * 
	 * @param realPart
	 * @param imaginaryPart
	 */
	public Complex(double realPart, double imaginaryPart){
		this.real = realPart;
		this.imag = imaginaryPart;
	}
	
//### operators for complex numbers
	public void add(Complex c){
		this.real += c.getReal();
		this.imag += c.getImag();
	}
	
	public void sub(Complex c){
		this.real -= c.getReal(); 
		this.imag -= c.getImag();
	}
	
	public void mul(Complex c){
		double re = this.real;
		double im = this.imag;
		this.real = re * c.getReal() - im * c.getImag();
		this.imag = re * c.getImag() + im * c.getReal();
	}
	
	public void div(Complex c){
		double re = this.real;
		double im = this.imag;
		double  d = c.getReal() * c.getReal() + c.getImag() * c.getImag();
		this.real = (re * c.getReal() + im * c.getImag())/d;
		this.imag = (im * c.getReal() - re * c.getImag())/d;
	}

//	### operators for doubles
	public void add(double z){
		this.real += z;
	}
	
	public void sub(double z){
		this.real -= z; 
	}
	
	public void mul(double z){
		this.real *= z;
		this.imag *= z;
	}
	
	public void div(double z){
		this.real /= z;
		this.imag /= z;
	}
	
//	### setting/getting
	
	public void setReal(double realPart){
		this.real = realPart;
	}

	public void setImag(double imaginaryPart){
		this.imag = imaginaryPart;
	}
	
	public double getReal(){
		return this.real;
	}
	
	public double getImag(){
		return this.imag;
	}
}
