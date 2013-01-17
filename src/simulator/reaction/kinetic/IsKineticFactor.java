/**
 * Project iDynoMiCS (copyright -> see Idynomics.java)
 */

package simulator.reaction.kinetic;

import java.io.Serializable;

import org.jdom.Element;


@SuppressWarnings("serial")
public abstract class IsKineticFactor implements Serializable{

	public int nParam;
	
	public abstract void init(Element defMarkUp);
	public abstract void initFromAgent(Element aReactionRoot,double[] kineticParam,int paramIndex);
	
	public abstract double kineticValue(double solute);
	public abstract double kineticDiff(double solute);
	public abstract double kineticMax();
	
	public abstract double kineticValue(double solute,double[] paramTable,int index);
	public abstract double kineticDiff(double solute,double[] paramTable,int index);
	public abstract double kineticMax(double[] paramTable,int index);
	
	
}
