
/**
 * Project iDynoMiCS (copyright -> see Idynomics.java)
 * 
 */

package simulator.agent;

import simulator.reaction.Reaction;

public interface HasReaction {

	public void addActiveReaction(Reaction aReaction, boolean useDefaultParam);

	public void addReaction(Reaction aReaction, boolean useDefaultParam);

	public void removeReaction(Reaction aPathway);

	public void switchOffreaction(Reaction aPathway);

	public void switchOnReaction(Reaction aPathway);

}
