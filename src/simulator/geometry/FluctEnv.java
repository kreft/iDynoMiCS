package simulator.geometry;

import java.io.*;
import java.util.*;
import org.jdom.Element;
import Jama.Matrix;

import idyno.SimTimer;
import simulator.Simulator;
import utils.XMLParser;

/**
 * This class is intended to randomly change the external environment felt by the agents. The transition between these
 * different states is determined by a transition probability matrix whose values are in probability per hour.
 * The correspondent cumulative probability matrix is used to decide the next state, simply because it avoids
 * the problem of having 2 equal values (and thus having to decide randomly between them...). 
 * It has also been implemented the recalculation of the cumulative probability values according to the 
 * current time step (as long as it is lower or equal to 1h... since higher values would produce probabilities 
 * higher than 1...)
 * 
 * It has also been implemented a method to deterministically change the environmental conditions and their duration.
 *  
 * @author SoniaMartins
 *
 */
public class FluctEnv extends Bulk {

	/***--------------------------------------------------------------***/
	//Variables used in the STOCHASTIC version
	
	//counters to use the parameters read from a file only in the first step
	static int counterA = 0;
	static int indcounter =0;
	static int indcounter2 =0;
	static int envCounter =0;
	
	Random num = new Random();
	
	public static ArrayList<String> envName = new ArrayList<String>() ;
	public static double [][] envCumProb;
	public static double [][] transitions;
	
   	// arrayList containing the indexes of the transitions array corresponding to the
   	// transition probability chosen.
	public static ArrayList<Integer> indices = new ArrayList<Integer>();
	
	//arrayList containing the names of the environments involved in the transition, in the correct order...
	public static ArrayList<String> envTrans = new ArrayList<String>();
	
	//declaration and initialisation of important variables
	public static double cumProb =0;
  	public static double transProb = 0;
  	private static int currentPos = 0;
  	private static int startPos = 0;
 
  	
  	/***--------------------------------------------------------------***/
  	
  	//Variables used in the DETERMiNiSTiC version 
  	
	public static HashMap<String, Double> envListTime = new HashMap<String, Double>();
	public static ArrayList<String> envNameList = new ArrayList<String>();
	public static HashMap<String,Integer> envIter = new HashMap<String, Integer>();

	public static String envStatus;
	public static double counter = 0;
	
	
	
	/**public FluctEnv(Simulator aSim, String name) {
		super(aSim, name);
	
	}*/
	
	public FluctEnv (Simulator aSim, XMLParser aBulkRoot) throws FileNotFoundException {
		super (aSim, aBulkRoot);

		LinkedList<Element> envList = aBulkRoot.buildSetMarkUp("envTime");
		XMLParser parser;
		
		for (Element aEnvMarkUp : envList) {
		parser = new XMLParser(aEnvMarkUp);
			envNameList.add(parser.getAttribute("name"));
			String s = parser.getAttribute("name");
			envListTime.put(s, parser.getAttributeDbl("time"));
	}
		for (int i=0; i<envNameList.size(); i++){
			Double val = envListTime.get(envNameList.get(i))/SimTimer.getCurrentTimeStep();
			int iter = val.intValue();
			envIter.put(envNameList.get(i), iter);
		}
		
	}

	
	/** 
     * reads a text file containing the names of the environments involved and returns a list with those names
     * @return ArrayList
     * @throws FileNotFoundException
     */
     public ArrayList<String> envDic()throws FileNotFoundException {
      
            try {
               
                FileReader file = new FileReader("EnvDic.txt");
                BufferedReader br = new BufferedReader (file);
                StreamTokenizer tokens = new StreamTokenizer (br);
                
                  while (tokens.nextToken()!= StreamTokenizer.TT_EOF){
                      
                      switch(tokens.ttype) {
                     
                      case StreamTokenizer.TT_WORD:
                        
                          envName.add(tokens.sval);
                          break;
                      }
                  }    
                 
              } catch (FileNotFoundException e) {
                      e.printStackTrace();
                    } catch (IOException e) {
                      e.printStackTrace();
                            }
           return envName;
        
     }
	

     /** 
      * reads a txt file containing a "table" with the transitions probabilities and returns a 2D array of it
      * @return double 2D array transitions
      * @throws FileNotFoundException
      */
     public double[][] envProb() throws FileNotFoundException {
    
     transitions = new double [envName.size()] [envName.size()];
     
      try {
              
              FileReader file = new FileReader("EnvProbData.txt");
              BufferedReader br = new BufferedReader (file);
              StreamTokenizer tokens = new StreamTokenizer (br);
   
               
              for (int row=0; row<transitions.length; row++) {
                  for(int column=0; column<transitions[row].length; column++){
                     
                      if (tokens.nextToken()!= StreamTokenizer.TT_EOF){
                
                          switch(tokens.ttype) {
                              case StreamTokenizer.TT_NUMBER:
                                  transitions[row][column]=tokens.nval;
                                  break;
                              default:
                                  break;
                                       
                          }
                      }        
                  }
             }          
            
       } catch (FileNotFoundException e) {
           e.printStackTrace();
       } catch (IOException e) {
           e.printStackTrace();
       }

    return transitions;
 }
     
     
     /**
      * reads a txt file containing the cumulative probability matrix describing the transitions between environments.
      * @return double 2D array envCumProb
      * @throws FileNotFoundException
      */
     public double [][] envCumProb() throws FileNotFoundException{

    	    
        envCumProb = new double [envName.size()] [envName.size()];
         
          try {
                  
                  FileReader file = new FileReader("EnvCumProb.txt");
                  BufferedReader br = new BufferedReader (file);
                  StreamTokenizer tokens = new StreamTokenizer (br);
       
                   
                  for (int row=0; row<envCumProb.length; row++) {
                      for(int column=0; column<envCumProb[row].length; column++){
                         
                          if (tokens.nextToken()!= StreamTokenizer.TT_EOF){
                    
                              switch(tokens.ttype) {
                                  case StreamTokenizer.TT_NUMBER:
                                      envCumProb[row][column]=tokens.nval;
                                      break;
                                  default:
                                      break;
                                           
                              }
                          }        
                      }
                 }          
                
           } catch (FileNotFoundException e) {
               e.printStackTrace();
           } catch (IOException e) {
               e.printStackTrace();
           }
           
           // recalculate cum prob values according to timeStep value
           
           double timeStep = SimTimer.getCurrentTimeStep();
           
           if(timeStep<1){
        	   for (int r=0; r<envCumProb.length; r++){
        		   for (int c=0; c<envCumProb.length; c++){
        			   envCumProb[r][c]=envCumProb[r][c]*(1/timeStep);
        		   }
        	   }
           }
           
           //else... think how to use a timestep above 1h...and having probabilities hogher than 1...
     
           
      return envCumProb;
    	 
     }
     
     /** 
      * transforms the 2D array containing the transition probabilities into a matrix whose format can be used to perform 
      * matrix operations
      * @return Matrix
      * @throws FileNotFoundException
      */
     public Matrix transMatrix() throws FileNotFoundException{

    	 double [][] tarr = this.envProb();
         Matrix ma = new Matrix(tarr);
         return ma;
     }
     
    /**
     * this method sets the new environment state according to the condition: 1st cumulative probability>= random number;
     * it returns the cumulative probability value that has been chosen;  
     * This is called in the step() in Simulator class and it is by default commented out.
     * @return double cumpr_val
     * @throws FileNotFoundException
     */
	public double setEnv() throws FileNotFoundException{

        	double rand = num.nextDouble();
 	
        	System.out.println( "random number is" + rand);
        	
        	//TODO read startPos from XMLfile
        	
        	if (counterA == 0){
        	loop1:	
        	for (int i=startPos; i<envCumProb[startPos].length; i++){
        		for ( int j=0; j<envCumProb.length; j++){
        			if( envCumProb[startPos][j] > rand){
        				cumProb = (Double)envCumProb[startPos][j];
        				counterA++;
        				break loop1;		
        				}
        			}
        		}
	
        	}else {	
        	currentPos = indices.get(1);
        	loop2:
            for (int row=currentPos; row<envCumProb[currentPos].length; row++){
            	for ( int column=0; column<envCumProb.length; column++){
            		if( envCumProb[currentPos][column] > rand){
            			cumProb = (Double)envCumProb[currentPos][column];
            			break loop2;		
            		}	
            	}
            }
        	}
        	return cumProb;
		}
   
	
	/**
	 * knowing the cumpr_val and the which row of probabilities is being used for this iteration, we can infer
	 * the indices of the cumulative probability matrix that correspond to the environment transition.
	 * 
	 * @return arrayList indices
	 * @throws FileNotFoundException
	 */
    public ArrayList<Integer> indexCalc() throws FileNotFoundException{
       	
       	//I had to add values to the empty list so that I can set them... instead of keeping adding eaht timestep...
       	if (indcounter == 0){
       	indices.add(1);
       	indices.add(2);
      	indcounter++;
       	}

        
        if (indcounter2 == 0){
         loop3:
         for (int r=startPos; r<envCumProb[startPos].length; r++){
             for (int c=0; c<envCumProb[r].length; c++){
                 if((Double) cumProb == envCumProb[startPos][c]){
                        indices.set(0,startPos);
                        indices.set(1,c);  
                        indcounter2++;             
                        break loop3;
                 	}
             	}
             }
        }else{
        	loop4:
            for (int r=currentPos; r<envCumProb[currentPos].length; r++){
                for (int c=0; c<envCumProb[r].length; c++){
                    if((Double) cumProb == envCumProb[currentPos][c]){
                           indices.set(0,currentPos);
                           indices.set(1,c);    
                           break loop4;
                       }
                	}
                }
        }

    	
       	return indices;
        	
        }
        
    /**
     * we would also like to know to which transition probability value is this transition associated with; 
     * that's what is done here. Using the indices we can retrieve the transition probability value.  
     * @return
     * @throws FileNotFoundException
     */
    public double transProbVal() throws FileNotFoundException{

        	 	loop5:
        	 	for (int k=0; k<transitions.length; k++){
        	 		for(int l=0; l<transitions[k].length; l++){
        	 			ArrayList<Integer> ind = new ArrayList<Integer>();
        	 			ind.add(k);
        	 			ind.add(l);

        	 			if(ind.equals(indices)){
        	 				transProb=transitions[k][l];

        	 				break loop5;
        	 			}	       	 			
        	 		}	
        	 	}
			
			return transProb;
					
		}
	
		
		/**
		 * taking into account the retrieved indices, builds a list with the names of the environments involved in the 
		 * transitions in the correct order
		 * @return ArrayList
		 * @throws FileNotFoundException
		 */
	public String envTransName() throws FileNotFoundException{
			
		 	if (envCounter == 0){
		       	envTrans.add("a");
		       	envTrans.add("b");
		      	envCounter++;
		       	}
			
			for(int i=0; i<indices.size(); i++){
				
				String envN = envName.get(indices.get(i));
				envTrans.set(i, envN);
			}
			
			envStatus = envTrans.get(1); 
		
			
			System.out.println(envTrans);
			return envStatus;
			
		}
	
	
	/***--------------------------------------------------------------***/
	                 // Deterministic Version //
	/**
	 * This method is the deterministic version of environment fluctuations, where different environments have
	 * defined time durations thus cycling between them at defined time points.
	 * @param index
	 * @return the environment Name just set.
	 */
		public static String setEnvCycle(int index){
		
		double timeStep = SimTimer.getCurrentTimeStep();
		
		double timeNow = SimTimer.getCurrentTime();
		int i;
		if(timeNow == 1){
			i = 0;
		}else{
			i = index;
		}
	
		//outer:
			for(String env = envNameList.get(i); i<envNameList.size() ; i++){
				//System.out.println("i is going from..." + i);
				double envTime = envListTime.get(env);
				if (counter < envTime){
					envStatus = envNameList.get(i);
					counter += timeStep;
					break;
					}else {
						counter = 0;
						if(i==envNameList.size()-1){
							i=-1;
						}
						//continue outer;
					}
					}
	
		System.out.println("Environment status is " + envStatus);
			
		return envStatus;
	
	
	}

		
	
}
