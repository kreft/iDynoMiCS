/**
 * Project iDynoMiCS (copyright -> see Idynomics.java)
 * 
 */

package simulator.agent;

public interface HasSpecies {

	/**
	 * Called by a progenitor during initialisation phases to create several
	 * agents Create an Agent copied with random deviation of the parameters and
	 * register it inside the pertinent containers
	 * @param position
	 * @throws CloneNotSupportedException
	 */
	public void createNewAgent()
	throws CloneNotSupportedException;

	public HasSpecies sendNewAgent() throws CloneNotSupportedException;
	public void setSpecies(Species aSpecies);
	public String sendHeader();
	public SpeciesParam getSpeciesParam();

	/**
	 * Mutate the parameters according the speciesParameters	
	 */
	public void mutatePop();
}
