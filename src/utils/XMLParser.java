/**
 * \package utils
 * \brief Package of classes that perform utility functions in the process of running an iDynoMiCS Simulation
 * 
 * Package of classes that perform utility functions in the process of running an iDynoMiCS Simulation. This package is part of iDynoMiCS v1.2, governed by the 
 * CeCILL license under French law and abides by the rules of distribution of free software.  You can use, modify and/ or redistribute 
 * iDynoMiCS under the terms of the CeCILL license as circulated by CEA, CNRS and INRIA at the following URL  "http://www.cecill.info".
 */
package utils;

import java.io.File;
import java.io.Serializable;
import java.util.LinkedList;
import java.util.List;
import java.util.ListIterator;
import org.jdom.Element;
import org.jdom.Document;
import org.jdom.input.SAXBuilder;
import simulator.geometry.DiscreteVector;

/**
 * \brief Class providing methods of parsing the XML simulation protocol file
 * 
 * In iDynoMiCS, protocol files are used to create the conditions under which a simulation is run. This is specified in XML. During 
 * initialisation, the simulation will require access to tags within this file. This class provides a means of accessing the data in 
 * these tags 
 */
public class XMLParser implements Serializable 
{

	/**
	 * Serial version used for the serialisation of the class
	 */
	private static final long   serialVersionUID = 1L;
	
	/**
	 * Document Object Model used to process an XML file
	 */
	private Document            document;
	
	/**
	 * The XML head tag of the group of elements currently being processed
	 */
	private Element             _localRoot;
	
	/**
	 * Temporary store of the value read in from the XML file. Used in cases where this value needs to be post-processed before returning 
	 */
	private static double       value;
	
	/**
	 * The unit of the protocol file parameter being processed
	 */
	private static StringBuffer unit;


	/**
	 * \brief Create an XML parser for a given XML local grouping
	 * 
	 * Create an XML parser for a given XML local grouping
	 * 
	 * @param localRoot	An element that is the root of a group of XML tags
	 */
	public XMLParser(Element localRoot) 
	{
		_localRoot = localRoot;
	}

	/**
	 * Create an XML parser object for an XML file of a given name
	 * 
	 * @param fileName	The name of the XML file to be processed
	 */
	public XMLParser(String fileName) 
	{
		openXMLDocument(fileName, false);
		_localRoot = document.getRootElement();
	}
	
	/**
	 * \brief Used to parse an XML file before a simulation log file is generated, when the simulation is being initialised
	 * 
	 * There are some checks the simulation needs to do on initialisation (with the protocol file). This is before the simulation log 
	 * is generated. These errors get stored in a separate log file
	 * 
	 * @param activePath	The path to the protocol file to open
	 * @param protocolFile	The protocol file to read in
	 */
	public XMLParser(String activePath,String protocolFile) 
	{
		openXMLDocument(activePath,protocolFile, false);
		_localRoot = document.getRootElement();
	}

	/**
	 * \brief Creates a DOM document for an XML protocol file of a given name
	 * 
	 * Creates a DOM document for an XML protocol file of a given name
	 * 
	 * @param fileName	The name of the XML protocol file being processed
	 * @param testDTD	Boolean stating whether or not schema validation should be applied on reading the XML file
	 */
	public void openXMLDocument(String fileName, boolean testDTD) {
		try {
			document = (new SAXBuilder(testDTD)).build(new File(fileName));
		} catch (Exception e) 
		{
			LogFile.writeLog("XMLParser.openXMLdocument(): Initialisation of the XML parser failed");
			LogFile.writeLog("Error: "+e);
			System.exit(-1);
		}
	}
	
	/**
	 * \brief Used to open an XML file in simulation initialisation - where the log file has yet to be created
	 * 
	 * Used to open an XML file in simulation initialisation - where the log file has yet to be created. This is the case as the 
	 * simulation log file location is specified in the protocol file. If that can't be read then there is a problem
	 * 
	 * @param activePath	String of path to protocol file being processed
	 * @param protocolFile	The protocol file to be read
	 * @param testDTD	Boolean stating whether or not schema validation should be applied on reading the XML file
	 */
	public void openXMLDocument(String activePath,String protocolFile, boolean testDTD) 
	{
		try 
		{
			document = (new SAXBuilder(testDTD)).build(new File(activePath+protocolFile));
		} catch (Exception e) 
		{
			InitialisationErrorLog.writeLog("XMLParser.openXMLdocument(): Initialisation of the XML parser failed");
			InitialisationErrorLog.writeLog("Check your XML files for Errors");
			InitialisationErrorLog.writeLog("Error: "+e);
			InitialisationErrorLog.closeFile();
			System.exit(-1);
		}
	}

	@SuppressWarnings("unchecked")
	/**
	 * \brief Returns a linked list of all XML elements of a given tag name
	 * 
	 * Returns a linked list of all XML elements within that protocol file that have a given tag name
	 * 
	 * @param childMarkup	The tag name of the elements to put in the list
	 * @return	A linked list of all the elements of the given tag name
	 */
	public LinkedList<Element> buildSetMarkUp(String childMarkup) 
	{
		List<Element> childList = _localRoot.getChildren(childMarkup);
		LinkedList<Element> out = new LinkedList<Element>();
		for (Object anElement : childList) 
		{
			out.add((Element) anElement);
		}
		return out;
	}

	@SuppressWarnings("unchecked")
	/**
	 * \brief Returns a list of XML parsers from combining nodes in the protocol file. Used for boundary conditions
	 * 
	 * Returns a list of XML parsers from combining nodes in the protocol file. Used for boundary conditions
	 * 
	 * @param childMarkup	The XML tag for which many nodes may exist. boundaryCondition is a good example
	 * @return	Linked list of XML parsers for each of the nodes of this name that are found in the protocol file
	 */
	public LinkedList<XMLParser> buildSetParser(String childMarkup) 
	{
		List<Element> childList = _localRoot.getChildren(childMarkup);
		LinkedList<XMLParser> out = new LinkedList<XMLParser>();
		for (Object anElement : childList) {
			out.add(new XMLParser((Element) anElement));
		}
		return out;
	}

	/* ___________________ TO READ PARAMETERS MARK-UPS _____________________ */

	/**
	 * \brief Returns the child node for a given XML tag name
	 * 
	 * Returns the child node for a given XML tag name
	 * 
	 * @param childName	The XML tag name to return
	 * @return	An XML parser containing the tags relating to this XML child tag
	 */
	public XMLParser getChild(String childName) {
		return new XMLParser(getChildElement(childName));
	}

	@SuppressWarnings("unchecked")
	/**
	 * \brief Searches the child nodes of a given tag for a particular parameter name. Returns the assigned string value if present
	 * 
	 * Searches the child nodes of a given tag for a particular parameter name. Returns the assigned string value if present
	 * 
	 * @param paramName	The name of the parameter for which a String value is being sought
	 * @return	String value assigned to that parameter, if present
	 */
	public String getParam(String paramName) 
	{
		List<Element> childList = _localRoot.getChildren("param");
		Element aParam;
		for (Object aChild : childList) 
		{
			aParam = (Element) aChild;
			if (aParam.getAttributeValue("name").equals(paramName)) { return aParam.getText(); }
		}
		return null;
	}
	
	@SuppressWarnings("unchecked")
	/**
	 * \brief Searches the child nodes of a given tag for a particular parameter name. Returns the XML element for that tag if present
	 * 
	 * Searches the child nodes of a given tag for a particular parameter name. Returns the XML element for that tag if present
	 * 
	 * @param paramName	The name of the parameter for which a String value is being sought
	 * @return	XML Element contains the tags relating to that parameter, if present
	 */
	public Element getParamMarkUp(String paramName) 
	{
		List<Element> childList = _localRoot.getChildren("param");
		Element aParam;
		for (Object aChild : childList) {
			aParam = (Element) aChild;
			if (aParam.getAttributeValue("name").equals(paramName)) { return aParam; }
		}
		return null;
	}
	
	@SuppressWarnings("unchecked")
	/**
	 * \brief Searches the child nodes of a given tag for a particular parameter name and a specified unit. Returns the String value of that tag if present
	 * 
	 * Searches the child nodes of a given tag for a particular parameter name and a specified unit. Returns the String value of that tag if present
	 * 
	 * @param paramName	The name of the parameter for which a String value is being sought
	 * @param unit	The unit that this parameter should be
	 * @return	The value assigned to this tag in the protocol file, if present
	 */
	public String getParam(String paramName, StringBuffer unit) {
		List<Element> childList = _localRoot.getChildren("param");
		Element aParam;
		for (Object aChild : childList) {
			aParam = (Element) aChild;
			if (aParam.getAttributeValue("name").equals(paramName)) {
				unit.append(aParam.getAttributeValue("unit"));
				return aParam.getText();
			}
		}
		return null;
	}

	/**
	 * \brief Returns the double value assigned to an XML tag in the protocol file
	 * 
	 * Returns the double value assigned to an XML tag in the protocol file, if the tag is present
	 * 
	 * @param paramName	The XML parameter for which the value is required
	 * @return	Double value assigned to that XML tag
	 */
	public Double getParamDbl(String paramName) 
	{
		if (getParam(paramName)==null) 
		{
			return Double.NaN;
		} 
		else 
		{
			return Double.parseDouble(getParam(paramName));
		}
	}

	/**
	 * \brief For a given XML tag name, returns the value in the specified unit (if that tag exists)
	 * 
	 * For a given XML tag name, returns the value in the specified unit (if that tag exists)
	 * 
	 * @param paramName	The name of the XML tag for which a value is required
	 * @param unit	The unit that this parameter is required to be within
	 * @return	The double value assigned to that parameter, in the required unit
	 */
	public Double getParamDbl(String paramName, StringBuffer unit) 
	{
		if (getParam(paramName)==null) 
		{
			return Double.NaN;
		} 
		else 
		{
			return Double.parseDouble(getParam(paramName, unit));
		}
	}

	/**
	 * \brief For a given XML tag name, returns the value it is assigned (if that tag exists)
	 * 
	 * For a given XML tag name, returns the value it is assigned (if that tag exists)
	 * 
	 * @param paramName	The name of the XML tag for which a value is required
	 * @return	The int value assigned to that parameter
	 */
	public Integer getParamInt(String paramName) 
	{
		if (getParam(paramName)==null) 
		{
			return 0;
		} 
		else 
		{
			return Integer.parseInt(getParam(paramName));
		}
	}

	/**
	 * \brief For a given XML tag name, and a given unit of measurement, returns the value the tag is assigned (if that tag exists)
	 * 
	 * For a given XML tag name, and a given unit of measurement, returns the value the tag is assigned (if that tag exists)
	 * 
	 * @param paramName	The name of the XML tag for which a value is required
	 * @param unit	The unif of measurement specified for this parameter
	 * @return	The int value assigned to that parameter
	 */
	public Integer getParamInt(String paramName, StringBuffer unit) {
		if (getParam(paramName)==null) {
			return 0;
		} else {
			return Integer.parseInt(getParam(paramName, unit));
		}
	}

	/**
	 * \brief Returns a length parameter from the XML, converting to the correct unit as required
	 * 
	 * Returns a length parameter from the XML, converting to the correct unit as required
	 * 
	 * @param paramName	The name of the parameter for which the length value should be returned
	 * @return	The length value assigned to this parameter in the protocol file
	 */
	public Double getParamLength(String paramName)
	{
		unit = new StringBuffer("");
		value = getParamDbl(paramName, unit);
		value *= utils.UnitConverter.length(unit.toString());
		return value;
	}

	/**
	 * \brief Returns the mass of a specified parameter from the XML, calculating this using the unit of measurement for that parameter
	 * 
	 * Returns the mass of a specified parameter from the XML, calculating this using the unit of measurement for that parameter
	 * 
	 * @param paramName	The name of the parameter for which the mass should be returned
	 * @return	The calculated mass of this parameter
	 */
	public Double getParamMass(String paramName) {
		unit = new StringBuffer("");
		value = getParamDbl(paramName, unit);
		value *= utils.UnitConverter.mass(unit.toString());
		return value;
	}

	/**
	 * \brief Gets a given parameter from the protocol file and converts into a double representing time
	 * 
	 * Gets a given parameter from the protocol file and converts into a double representing time
	 * 
	 * @param paramName	The parameter to be retrieved from the XML file
	 * @return	Double of the value of this parameter converted into a unit of time
	 */
	public Double getParamTime(String paramName) {
		unit = new StringBuffer("");
		value = getParamDbl(paramName, unit);
		value *= utils.UnitConverter.time(unit.toString());
		return value;
	}

	/**
	 * \brief Reads a concentration from the protocol file and converts this to the required unit
	 * 
	 * Reads a concentration from the protocol file and converts this to the required unit
	 * 
	 * @param paramName	The concentration of the parameter to be retrieved from the protocol file
	 * @return	The calculated concentration level for this simulation parameter
	 */
	public Double getParamConc(String paramName) {
		unit = new StringBuffer("");
		value = getParamDbl(paramName, unit);
		value *= utils.UnitConverter.mass(unit.toString());
		value *= utils.UnitConverter.volume(unit.toString());
		return value;
	}

	/**
	 * \brief Converts coordinates specified in the protocol file into a discrete vector
	 * 
	 * Converts coordinates specified in the protocol file into a discrete vector. Used for boundary conditions and specifications of 
	 * birth area of agents
	 * 
	 * @param paramName	The protocol file parameter specified as X,Y,and Z coordinated
	 * @return	A discrete vector containing these three coordinates
	 */
	public DiscreteVector getParamXYZ(String paramName){
		//Element param = getParamMarkUp(paramName);
		return new DiscreteVector(getParamMarkUp(paramName));
	}
	
	/**
	 * \brief Read in boolean value of specified parameter from the XML file
	 * 
	 * Reads in boolean value of specified parameter from the XML file
	 * 
	 * @param paramName	The parameter to be retrieved from the XML file
	 * @return	Boolean value assigned to this parameter
	 */
	public Boolean getParamBool(String paramName) {
		return Boolean.parseBoolean(getParam(paramName));
	}

	/**
	 * \brief Read in boolean value and unit of measurement of a specified parameter from the XML file
	 * 
	 * Read in boolean value and unit of measurement of a specified parameter from the XML file
	 * 
	 * @param paramName	The parameter to be retrieved from the XML file
	 * @param unit	The unit of measurement specified for this parameter
	 * @return	Boolean value assigned to this parameter
	 */
	public Boolean getParamBool(String paramName, StringBuffer unit) {
		return Boolean.parseBoolean(getParam(paramName, unit));
	}

	@SuppressWarnings("unchecked")
	/**
	 * \brief Searches through the attributes of the XML tags of a given parameter name to find the String value of a specified detail within that tag
	 * 
	 * Searches through the attributes of the XML tags of a given parameter name to find the value of a specified detail within that tag
	 * 
	 * @param paramName	The parameter name for which the value is required
	 * @param detailName	The name of the detail element which is part of that tag, if present
	 * @return	The text associated with the specified detail of this tag
	 */
	public String getParamSuch(String paramName, String detailName) {

		List<Element> childList = _localRoot.getChildren("param");
		Element aParam;
		for (Object aChild : childList) {
			aParam = (Element) aChild;
			if (aParam.getAttributeValue("name").equals(paramName)) {
				if (aParam.getAttributeValue("detail").equals(detailName)) return aParam.getText();
			} else {
				continue;
			}
		}
		return null;
	}

	@SuppressWarnings("unchecked")
	/**
	 * \brief Searches through the attributes of the XML tags of a given parameter name, in a specified unit to find the String value of a specified detail within that tag
	 * 
	 * Searches through the attributes of the XML tags of a given parameter name to find the value of a specified detail within that tag
	 * 
	 * @param paramName	The parameter name for which the value is required
	 * @param detailName	The name of the detail element which is part of that tag, if present
	 * @param unit	The unit that this parameter should be measured in
	 * @return	The text associated with the specified detail of this tag
	 */
	public String getParamSuch(String paramName, String detailName, StringBuffer unit) {
		List<Element> childList = _localRoot.getChildren("param");
		Element aParam;
		for (Object aChild : childList) {
			aParam = (Element) aChild;
			if (aParam.getAttributeValue("name").equals(paramName)) {
				if (aParam.getAttributeValue("detail").equals(detailName)) {
					unit.append(aParam.getAttributeValue("unit"));
					return aParam.getText();
				}
			} else {
				return null;
			}
		}
		return null;
	}

	/**
	 * \brief Searches through the attributes of the XML tags of a given parameter name to find the Double value of a specified detail within that tag
	 * 
	 * Searches through the attributes of the XML tags of a given parameter name to find the Double value of a specified detail within that tag
	 * 
	 * @param paramName	The parameter name for which the value is required
	 * @param detailName	The name of the detail element which is part of that tag, if present
	 * @return	The double value associated with the specified detail of this tag
	 */
	public Double getParamSuchDbl(String paramName, String detailName) {
		if (getParamSuch(paramName, detailName)==null) {
			return Double.NaN;
		} else {
			return Double.parseDouble(getParamSuch(paramName, detailName));
		}
	}

	/**
	 * \brief Searches through the attributes of the XML tags of a given parameter name to find the boolean value of a specified detail within that tag
	 * 
	 * Searches through the attributes of the XML tags of a given parameter name to find the boolean value of a specified detail within that tag
	 * 
	 * @param paramName	The parameter name for which the value is required
	 * @param detailName	The name of the detail element which is part of that tag, if present
	 * @return	The boolean value associated with the specified detail of this tag
	 */
	public Boolean getParamSuchBool(String paramName, String detailName) {
		if (getParamSuch(paramName, detailName)==null) {
			return null;
		} else {
			return Boolean.parseBoolean(getParamSuch(paramName, detailName));
		}
	}

	/**
	 * \brief Searches through the attributes of the XML tags of a given parameter name, in a specified unit to find the Double value of a specified detail within that tag
	 * 
	 * Searches through the attributes of the XML tags of a given parameter name to find the Double value of a specified detail within that tag
	 * 
	 * @param paramName	The parameter name for which the value is required
	 * @param detailName	The name of the detail element which is part of that tag, if present
	 * @param unit	The unit that this parameter should be measured in
	 * @return	The double value associated with the specified detail of this tag
	 */
	public double getParamSuchDbl(String paramName, String detailName, StringBuffer unit) {
		if (getParamSuch(paramName, detailName)==null) {
			return Double.NaN;
		} else {
			return Double.parseDouble(getParamSuch(paramName, detailName, unit));
		}
	}

	/**
	 * \brief Searches through the attributes of the XML tags of a given parameter name, in a specified unit to find the boolean value of a specified detail within that tag
	 * 
	 * Searches through the attributes of the XML tags of a given parameter name to find the boolean value of a specified detail within that tag
	 * 
	 * @param paramName	The parameter name for which the value is required
	 * @param detailName	The name of the detail element which is part of that tag, if present
	 * @param unit	The unit that this parameter should be measured in
	 * @return	The boolean value associated with the specified detail of this tag
	 */
	public Boolean getParamSuchBool(String paramName, String detailName, StringBuffer unit) {
		if (getParamSuch(paramName, detailName)==null) {
			return null;
		} else {
			return Boolean.parseBoolean(getParamSuch(paramName, detailName, unit));
		}
	}

	
	@SuppressWarnings("unchecked")
	/**
	 * \brief Returns the XML element matching the specified tag name, attribute name, and attribute value. If not find return the local root
	 * 
	 * Returns the XML element matching the specified tag name, attribute name, and attribute value. If not find return the local root
	 * 
	 * @param childName	The XML tag required
	 * @param attrName	The name of the attribute within that tag
	 * @param attrValue	The value of the attribute within that tag
	 * @return	The XML element that matches these three criteria
	 */
	public Element getChildSuchAttribute(String childName, String attrName, String attrValue) {
		List<Element> childList = _localRoot.getChildren(childName);
		for (Object aChild : childList) {
			if (((Element) aChild).getAttributeValue(attrName).equals(attrValue)) { return (Element) aChild; }
		}
		return _localRoot;
	}

	@SuppressWarnings("unchecked")
	
	
	/**
	 * \brief Return the string value of an attribute of the current localroot of the XML file
	 * 
	 * Return the string value of an attribute of the current localroot of the XML file
	 * 
	 * @param attributeName	The attribute name for which the value is being sought
	 * @return	The string value of that attribute, if present
	 */
	public String getAttributeStr(String attributeName) {
		return _localRoot.getAttributeValue(attributeName);
	}

	/**
	 * \brief Return the double value of an attribute of the current localroot of the XML file
	 * 
	 * Return the double value of an attribute of the current localroot of the XML file
	 * 
	 * @param attributeName	The attribute name for which the value is being sought
	 * @return	The double value of that attribute, if present
	 */
	public Double getAttributeDbl(String attributeName) 
	{
		// KA NOV 13 - Added a return of NaN if the attribute was not found. Required when setting up the cell birth area by volume/area
		try
		{
			return Double.parseDouble(getAttributeStr(attributeName));
		}
		catch(Exception e)
		{
			return Double.NaN;
		}
	}

	/**
	 * \brief Return the string value of a child attribute of the current localroot of the XML file
	 * 
	 * Return the string value of a child attribute of the current localroot of the XML file
	 * 
	 * @param childName	Name of the XML child tag
	 * @param attrName	The attribute of that tag for which the value is being sought
	 * @return	The string value of that attribute, if present
	 */
	public String getChildAttrStr(String childName, String attrName) {
		try {
			return _localRoot.getChild(childName).getAttributeValue(attrName);
		} catch (Exception e) {
			return null;
		}
	}

	/**
	 * \brief Returns the value assigned to an attribute within a child node
	 * 
	 * Returns the value assigned to an attribute within a child node. Checking computation domain dimension is one example of its use
	 * 
	 * @param childName	The name of the child node of which the attribute value is being sought
	 * @param attrName	The name of the attribute for which the value is required
	 * @return	Double value assigned to that attribute, if present
	 */
	public Double getChildAttrDbl(String childName, String attrName) 
	{
		String value = getChildAttrStr(childName, attrName);
		if (value==null) return 0.0;
		else return Double.parseDouble(value);
	}
	
	/**
	 * \brief Returns the value assigned to an attribute within a child node
	 * 
	 * Returns the value assigned to an attribute within a child node. Checking computation domain dimension is one example of its use
	 * 
	 * @param childName	The name of the child node of which the attribute value is being sought
	 * @param attrName	The name of the attribute for which the value is required
	 * @return	Double value assigned to that attribute, if present
	 */
	public Integer getChildAttrInt(String childName, String attrName) 
	{
		String value = getChildAttrStr(childName, attrName);
		if (value==null) return 0;
		else return Integer.parseInt(value);
	}
	
	@SuppressWarnings("unchecked")
	/**
	 * \brief Returns a list of the children of a given XML child tag name
	 * 
	 * Returns a list of the children of a given XML child tag name
	 * 
	 * @param childName	The name of a child tag in the XML file for which the list of children is required
	 * @return	List of XML elements containing the children of the specified tag name
	 */
	public List<Element> getChildren(String childName) {
		return (List<Element>) _localRoot.getChildren(childName);
	}

	/**
	 * \brief Returns the XML element of a given child node tag name
	 * 
	 * Returns the XML element of a given child node tag name
	 * 
	 * @param childName	The child tag name for which the XML element is required
	 * @return	An XML element containing the tags for the specified child name
	 */
	public Element getChildElement(String childName) {
		return _localRoot.getChild(childName);
	}

	/**
	 * \brief Return the value assigned to a given XML tag
	 * 
	 * Returns the value assigned to a given XML tag in a protocol file
	 * 
	 * @param attributeName	The XML tag for which the associated value is needed
	 * @return	The string value assigned to that tag in the protocol file
	 */
	public String getAttribute(String attributeName) 
	{
		return _localRoot.getAttributeValue(attributeName);
	}

	/**
	 * \brief Return the current XML local root element
	 * 
	 * Return the current XML local root element
	 * 
	 * @return	XML Element for the current local root
	 */
	public Element getElement() {
		return _localRoot;
	}

	/**
	 * \brief Creates an instance of a class using a string containing that class name
	 * 
	 * Creates an instance of a class using a string containing that class name. Useful for a set of boundary conditions, for example
	 * 
	 * @param prefix	The class for which a new instance is being created
	 * @return	An object of the class stated in the prefix string
	 */
	public Object instanceCreator(String prefix) 
	{
		prefix += "."+this.getAttribute("class");

		try 
		{
			return Class.forName(prefix).newInstance();
		} 
		catch (Exception e) 
		{
			LogFile.writeLog("Unable to create class");
			return null;
		}
	}
		

	/**
	 * \brief Used for Epi-Bac scenarios - retrieving list of environments to which this species is sensitive to and the correspondent probability of dying if under the influence of that environment
	 * 
	 * Used for Epi-Bac scenarios - retrieving list of environments to which this species is sensitive to and the correspondent probability of dying if under the influence of that environment
	 * 
	 * @param childName	The name of the XML tag that contains the information on this parameter
	 * @param attrName	The name of the first attribute within that XML tag
	 * @param attrValue	The value of the first attribute
	 * @param attr2Name	The name of the related attribute that is required to find this value
	 * @return	Double value assigned to this tag
	 */
	public Double getDblAttrOfChildSuchAttribute(String childName, String attrName,
	        String attrValue, String attr2Name) {
		String out = getChildSuchAttribute(childName, attrName, attrValue).getAttributeValue(
		        attr2Name);
		return Double.parseDouble(out);
	}

	

	
}
