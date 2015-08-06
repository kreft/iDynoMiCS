/**
 * \package utils
 * \brief Package of classes that perform utility functions in the process of
 * running an iDynoMiCS Simulation.
 * 
 * This package is part of iDynoMiCS v1.2, governed by the CeCILL license
 * under French law and abides by the rules of distribution of free software.
 * You can use, modify and/ or redistribute iDynoMiCS under the terms of the
 * CeCILL license as circulated by CEA, CNRS and INRIA at the following URL
 * "http://www.cecill.info".
 */
package utils;

import java.io.File;
import java.io.Serializable;
import java.util.LinkedList;
import java.util.List;

import org.jdom.Document;
import org.jdom.Element;
import org.jdom.input.SAXBuilder;

/**
 * \brief Class providing methods of parsing the XML simulation protocol file.
 * 
 * In iDynoMiCS, protocol files are used to create the conditions under which a
 * simulation is run. This is specified in XML. During initialisation, the
 * simulation will require access to tags within this file. This class provides
 * a means of accessing the data in these tags. 
 */
public class XMLParser implements Serializable 
{
	/**
	 * Serial version used for the serialisation of the class.
	 */
	private static final long serialVersionUID = 1L;
	
	/**
	 * Document Object Model used to process an XML file.
	 */
	private Document document;
	
	/**
	 * The XML head tag of the group of elements currently being processed.
	 */
	private Element _localRoot;
	
	/**
	 * Temporary store of the value read in from the XML file. Used in cases
	 * where this value needs to be post-processed before returning. 
	 */
	private static Double value;
	
	/**
	 * The unit of the protocol file parameter being processed.
	 */
	private static StringBuffer unit;
	
	/**
	 * The default value to set integers to, if there's a problem.
	 */
	public static Integer nullInt = null; 
	
	/**
	 * The default value to set doubles to, if there's a problem.
	 */
	public static Double nullDbl = Double.NaN;
	
	/**
	 * The default value to set booleans to, if there's a problem.
	 */
	public static Boolean nullBool = null;
	
	/*************************************************************************
	 * Creating the XMLParser
	 */
	
	/**
	 * \brief Create an XML parser for a given XML local grouping.
	 * 
	 * @param localRoot	An element that is the root of a group of XML tags.
	 */
	public XMLParser(Element localRoot) 
	{
		_localRoot = localRoot;
	}
	
	/**
	 * Create an XML parser object for an XML file of a given name.
	 * 
	 * @param fileName	The name of the XML file to be processed.
	 */
	public XMLParser(String fileName) 
	{
		openXMLDocument(fileName, false);
		_localRoot = document.getRootElement();
	}
	
	/**
	 * \brief Used to parse an XML file before a simulation log file is
	 * generated, when the simulation is being initialised.
	 * 
	 * There are some checks the simulation needs to do on initialisation (with
	 * the protocol file). This is before the simulation log is generated.
	 * These errors get stored in a separate log file.
	 * 
	 * @param activePath	The path to the protocol file to open.
	 * @param protocolFile	The protocol file to read in.
	 */
	public XMLParser(String activePath, String protocolFile) 
	{
		openXMLDocument(activePath, protocolFile, false);
		_localRoot = document.getRootElement();
	}
	
	/*************************************************************************
	 * Opening XML documents
	 */
	
	/**
	 * \brief Creates a DOM document for an XML protocol file of a given name.
	 * 
	 * @param fileName	The name of the XML protocol file being processed.
	 * @param testDTD	Boolean stating whether or not schema validation should
	 * be applied on reading the XML file.
	 */
	public void openXMLDocument(String fileName, Boolean testDTD)
	{
		try
		{
			document = (new SAXBuilder(testDTD)).build(new File(fileName));
		}
		catch (Exception e) 
		{
			LogFile.writeLogAlways("Error while trying to open "+fileName);
			LogFile.writeError(e, "XMLParser.openXMLdocument()");
			System.exit(-1);
		}
	}
	
	/**
	 * \brief Used to open an XML file in simulation initialisation, where the
	 * log file has yet to be created.
	 * 
	 * This is the case as the simulation log file location is specified in 
	 * the protocol file. If that can't be read then there is a problem.
	 * 
	 * @param activePath	String of path to protocol file being processed.
	 * @param protocolFile	The protocol file to be read.
	 * @param testDTD	Boolean stating whether or not schema validation 
	 * should be applied on reading the XML file.
	 */
	public void openXMLDocument(String activePath,
										String protocolFile, Boolean testDTD)
	{
		try 
		{
			document = (new SAXBuilder(testDTD))
									.build(new File(activePath+protocolFile));
		}
		catch (Exception e) 
		{
			InitialisationErrorLog.writeLog("XMLParser.openXMLdocument(): "+
								"Initialisation of the XML parser failed\n"+
								"Check your XML files for Errors\nError: "+e);
			InitialisationErrorLog.closeFile();
			System.exit(-1);
		}
	}
	
	/*************************************************************************
	 * Reading XML Elements: top-level methods
	 */
	
	/**
	 * \brief Return the current XML local root element.
	 * 
	 * @return	XML Element for the current local root.
	 */
	public Element getElement()
	{
		return _localRoot;
	}
	
	/**
	 * @return	String value of the local root.
	 */
	public String getValue()
	{
		return _localRoot.getValue();
	}
	
	/**
	 * \brief Return the string value of an attribute of the current localRoot
	 * of the XML file.
	 * 
	 * @param attributeName	The attribute name for which the value is being
	 * sought.
	 * @return	The string value of that attribute, if present.
	 */
	public String getAttribute(String attributeName)
	{
		return _localRoot.getAttributeValue(attributeName);
	}
	
	/**
	 * \brief Gets the name of this Element.
	 * 
	 * Same as getAttribute("name"). Created as used very often
	 * 
	 * @return The name of the Element, as a String.
	 */
	public String getName()
	{
		return getAttribute("name");
	}
	
	/**
	 * \brief Returns the XML element of a given child node tag name.
	 * 
	 * @param childName	The child tag name for which the XML element is
	 * required.
	 * @return	An XML element containing the tags for the specified child
	 * name.
	 */
	public Element getChildElement(String childName)
	{
		return _localRoot.getChild(childName);
	}
	
	/**
	 * \brief Return the string value of a child attribute of the current
	 * local root of the XML file.
	 * 
	 * @param childName	Name of the XML child tag.
	 * @param attrName	The attribute of that tag for which the value is being
	 * sought.
	 * @return	The string value of that attribute, if present.
	 */
	public String getChildAttrStr(String childName, String attrName)
	{
		try 
		{
			return getChildElement(childName).getAttributeValue(attrName);
		} 
		catch (Exception e) 
		{
			return null;
		}
	}
	
	/**
	 * \brief Returns the child node for a given XML tag name.
	 * 
	 * @param childName	The XML tag name to return.
	 * @return	An XML parser containing the tags relating to this XML child
	 * tag.
	 */
	public XMLParser getChildParser(String childName)
	{
		return new XMLParser(getChildElement(childName));
	}
	
	/**
	 * \brief Returns a list of the children of a given XML child mark up.
	 * 
	 * @param childMarkup	The name of a child tag in the XML file for which the
	 * list of children is required.
	 * @return	List of XML elements containing the children of the specified
	 * tag mark up.
	 */
	@SuppressWarnings("unchecked")
	public LinkedList<Element> getChildrenElements(String childMarkup)
	{
		List<Object> childList = _localRoot.getChildren(childMarkup);
		LinkedList<Element> out = new LinkedList<Element>();
		for (Object child : childList)
			out.add((Element) child);
		return out;
	}
	
	/**
	 * \brief Returns a list of names for the Child Elements with the given
	 * mark up.
	 * 
	 * @param childMarkup
	 * @return
	 */
	public LinkedList<String> getChildrenNames(String childMarkup)
	{
		LinkedList<String> out = new LinkedList<String>();
		for (XMLParser aChild : getChildrenParsers(childMarkup))
			out.add(aChild.getName());
		return out;
	}
	
	/**
	 * \brief Returns a list of XML parsers from combining nodes in the
	 * protocol file.
	 * 
	 * Used for boundary conditions.
	 * 
	 * @param childMarkup	The XML tag for which many nodes may exist.
	 * @return	Linked list of XML parsers for each of the nodes of this name
	 * that are found in the protocol file.
	 */
	public LinkedList<XMLParser> getChildrenParsers(String childMarkup) 
	{
		LinkedList<Element> childList = getChildrenElements(childMarkup);
		LinkedList<XMLParser> out = new LinkedList<XMLParser>();
		for (Element anElement : childList)
			out.add(new XMLParser(anElement));
		return out;
	}
	
	/**
	 * \brief Checks whether an XML element with the specified tag, attribute
	 * and value exists among the children
	 * 
	 * @param childMarkup
	 * @param attrName
	 * @param attrValue
	 * @return
	 */
	@SuppressWarnings("unchecked")
	public Boolean isChildSuchGiven(String childMarkup, 
									String attrName, String attrValue)
	{
		List<Element> childList = _localRoot.getChildren(childMarkup);
		for (Element aChild : childList)
			if ( aChild.getAttributeValue(attrName).equals(attrValue) ) 
				return true;
		return false;
	}
	
	/**
	 * \brief Returns the XML element matching the specified tag name,
	 * attribute name, and attribute value.
	 * 
	 * If not found returns the local root.
	 * 
	 * @param childMarkup	The XML tag required.
	 * @param attrName	The name of the attribute within that tag.
	 * @param attrValue	The value of the attribute within that tag.
	 * @return	The XML element that matches these three criteria.
	 */
	@SuppressWarnings("unchecked")
	public Element getChildSuchAttribute(String childMarkup, 
							String attrName, String attrValue)
	{
		List<Element> childList = _localRoot.getChildren(childMarkup);
		for (Element aChild : childList)
			if ( aChild.getAttributeValue(attrName).equals(attrValue) ) 
				return (Element) aChild;
		return null;
	}
	
	/*************************************************************************
	 * Reading parameters: generic methods
	 */
	
	/**
	 * TODO
	 * @param paramName
	 * @return
	 */
	public Boolean isParamGiven(String paramName)
	{
		return isChildSuchGiven("param", "name", paramName);
	}
	
	/**
	 * \brief Searches the child nodes of a given tag for a particular
	 * parameter name.
	 * 
	 * Returns the XML element for that tag if present.
	 * 
	 * @param paramName	The name of the parameter for which a String value is
	 * being sought.
	 * @return	XML Element contains the tags relating to that parameter, if
	 * present.
	 */
	public Element getParamElement(String paramName) 
	{
		return getChildSuchAttribute("param", "name", paramName);
	}
	
	
	public XMLParser getParamParser(String paramName)
	{
		return new XMLParser(getParamElement(paramName));
	}
	
	/**
	 * \brief Searches the child nodes of a given tag for a particular
	 * parameter name
	 * 
	 * Returns the assigned string value if present.
	 * 
	 * @param paramName	The name of the parameter for which a String value is
	 * being sought.
	 * @return	String value assigned to that parameter, if present
	 */
	public String getParam(String paramName) 
	{
		Element aParam = getParamElement(paramName);
		return ( aParam == null ) ? null : aParam.getText();
	}
	
	/**
	 * \brief Searches the child nodes of a given tag for a particular
	 * parameter name and a specified unit.
	 * 
	 * Returns the String value of that tag if present. If the parameter is 
	 * not specified this should not deal with the unit.
	 * 
	 * @param paramName	The name of the parameter for which a String value is
	 * being sought.
	 * @param unit	The unit that this parameter should be.
	 * @return	The value assigned to this tag in the protocol file, if
	 * present.
	 */
	public String getParam(String paramName, StringBuffer unit)
	{
		Element aParam = getParamElement(paramName);
		if ( aParam == null )
			return null;
		unit.append(aParam.getAttributeValue("unit"));
		return aParam.getText();
	}
	
	/**
	 * \brief Searches through the attributes of the XML tags of a given
	 * parameter name to find the String value of a specified detail within
	 * that tag.
	 * 
	 * @param paramName	The parameter name for which the value is required.
	 * @param detailName	The name of the detail element which is part of
	 * that tag, if present.
	 * @return	The text associated with the specified detail of this tag.
	 */
	public String getParamSuch(String paramName, String detailName)
	{
		Element aParam = getChildSuchAttribute("param", paramName, detailName);
		return ( aParam == null ) ? null : aParam.getText();
	}
	
	/**
	 * \brief Searches through the attributes of the XML tags of a given
	 * parameter name, in a specified unit to find the String value of a
	 * specified detail within that tag.
	 * 
	 * @param paramName	The parameter name for which the value is required.
	 * @param detailName	The name of the detail element which is part of
	 * that tag, if present.
	 * @param unit	The unit that this parameter should be measured in.
	 * @return	The text associated with the specified detail of this tag.
	 */
	public String getParamSuch(String paramName,
										String detailName, StringBuffer unit)
	{
		Element aParam = getChildSuchAttribute("param", paramName, detailName);
		if ( aParam == null )
			return null;
		unit.append(aParam.getAttributeValue("unit"));
		return aParam.getText();
	}
	
	/*************************************************************************
	 * Reading Integer values: if not present, return null
	 */
	
	/**
	 * \brief Converts the given string to an integer, in a clear and 
	 * consistent way.
	 * 
	 * @param in	String to be converted.
	 * @return	Integer value of this string, or nullInt if it cannot be.
	 */
	private Integer stringToInteger(String in)
	{
		if ( in == null || in.equals("") )
			return nullInt;
		return Integer.parseInt(in);
	}
	
	/**
	 * \brief Returns the value assigned to an attribute within a child node.
	 * 
	 * Checking computation domain dimension is one example of its use.
	 * 
	 * @param childName	The name of the child node of which the attribute value
	 * is being sought.
	 * @param attrName	The name of the attribute for which the value is
	 * required.
	 * @return	Integer value assigned to that attribute, if present.
	 */
	public Integer getChildAttrInt(String childName, String attrName) 
	{
		return stringToInteger(getChildAttrStr(childName, attrName));
	}
	
	/**
	 * \brief Return the Integer value of an attribute of the current localRoot
	 * of the XML file.
	 * 
	 * @param attributeName	The attribute name for which the value is being
	 * sought.
	 * @return	The Integer value of that attribute, if present.
	 */
	public Integer getAttributeInt(String attributeName)
	{
		return stringToInteger(getAttribute(attributeName));
	}
	
	/**
	 * \brief For a given XML tag name, returns the value it is assigned
	 * (if that tag exists).
	 * 
	 * @param paramName	The name of the XML tag for which a value is required.
	 * @return	The int value assigned to that parameter.
	 */
	public Integer getParamInt(String paramName) 
	{
		return stringToInteger(getParam(paramName));
	}
	
	/*************************************************************************
	 * Reading Boolean values: if not present, return null
	 */
	
	/**
	 * \brief Converts the given string to an boolean, in a clear and 
	 * consistent way.
	 * 
	 * @param in	String to be converted.
	 * @return	Boolean value of this string, or nullBool if it cannot be.
	 */
	private Boolean stringToBoolean(String in)
	{
		if ( in == null || in.equals("") )
			return nullBool;
		return Boolean.parseBoolean(in);
	}
	
	/**
	 * \brief Read in boolean value of specified parameter from the XML file.
	 * 
	 * @param paramName	The parameter to be retrieved from the XML file.
	 * @return	Boolean value assigned to this parameter.
	 */
	public Boolean getParamBool(String paramName)
	{
		return stringToBoolean(getParam(paramName));
	}
	
	/**
	 * \brief Searches through the attributes of the XML tags of a given
	 * parameter name to find the boolean value of a specified detail within
	 * that tag.
	 * 
	 * @param paramName	The parameter name for which the value is required.
	 * @param detailName	The name of the detail element which is part of
	 * that tag, if present.
	 * @return	The boolean value associated with the specified detail of this
	 * tag.
	 */
	public Boolean getParamSuchBool(String paramName, String detailName)
	{
		return stringToBoolean(getParamSuch(paramName, detailName));
	}
	
	/*************************************************************************
	 * Reading Double values: if not present, return Double.NaN
	 */
	
	/**
	 * \brief Converts the given string to a double, in a clear and consistent
	 * way.
	 * 
	 * @param in	String to be converted.
	 * @return	Double value of this string, or nullDbl if it cannot be.
	 */
	private Double stringToDouble(String in)
	{
		if ( in == null || in.equals("") )
			return nullDbl;
		return Double.parseDouble(in);
	}
	
	/**
	 * \brief Gets the value of this local root as a Double.
	 * 
	 * @return	Double value of the local root.
	 */
	public Double getValueDbl()
	{
		return stringToDouble(getValue());
	}
	
	/**
	 * \brief Return the Double value of an attribute of the current localRoot
	 * of the XML file.
	 * 
	 * @param attributeName	The attribute name for which the value is being
	 * sought.
	 * @return	The Double value of that attribute, if present.
	 */
	public Double getAttributeDbl(String attributeName)
	{
		return stringToDouble(getAttribute(attributeName));
	}
	
	/**
	 * \brief Returns the value assigned to an attribute within a child node.
	 * 
	 * Checking computation domain dimension is one example of its use.
	 * 
	 * @param childName	The name of the child node of which the attribute
	 * value is being sought.
	 * @param attrName	The name of the attribute for which the value is
	 * required.
	 * @return	Double value assigned to that attribute, if present
	 */
	public Double getChildAttrDbl(String childName, String attrName) 
	{
		return stringToDouble(getChildAttrStr(childName, attrName));
	}
	
	/**
	 * \brief Returns the double value assigned to an XML tag in the protocol
	 * file, if the tag is present.
	 * 
	 * @param paramName	The XML parameter for which the value is required
	 * @return	Double value assigned to that XML tag
	 */
	public Double getParamDbl(String paramName) 
	{
		return stringToDouble(getParam(paramName));
	}

	/**
	 * \brief For a given XML tag name, returns the value in the specified
	 * unit (if that tag exists).
	 * 
	 * @param paramName	The name of the XML tag for which a value is required.
	 * @param unit	The unit that this parameter is required to be within.
	 * @return	The double value assigned to that parameter, in the required
	 * unit.
	 */
	public Double getParamDbl(String paramName, StringBuffer unit) 
	{
		return stringToDouble(getParam(paramName, unit));
	}
	
	/**
	 * \brief Searches through the attributes of the XML tags of a given
	 * parameter name to find the Double value of a specified detail within
	 * that tag.
	 * 
	 * @param paramName	The parameter name for which the value is required.
	 * @param detailName	The name of the detail element which is part of
	 * that tag, if present.
	 * @return	The double value associated with the specified detail of this
	 * tag.
	 */
	public Double getParamSuchDbl(String paramName, String detailName)
	{
		return stringToDouble(getParamSuch(paramName, detailName));
	}
	
	/**
	 * \brief Searches through the attributes of the XML tags of a given
	 * parameter name, in a specified unit to find the Double value of a
	 * specified detail within that tag.
	 * 
	 * @param paramName	The parameter name for which the value is required.
	 * @param detailName	The name of the detail element which is part of
	 * that tag, if present.
	 * @param unit	The unit that this parameter should be measured in.
	 * @return	The double value associated with the specified detail of this
	 * tag.
	 */
	public Double getParamSuchDbl(String paramName,
										String detailName, StringBuffer unit)
	{
		return stringToDouble(getParamSuch(paramName, detailName, unit));
	}
	
	/**
	 * \brief Used for Epi-Bac scenarios.
	 * 
	 * Retrieving list of environments to which this species is
	 * sensitive to and the correspondent probability of dying if under
	 * the influence of that environment.
	 * 
	 * @param childName	The name of the XML tag that contains the information
	 * on this parameter.
	 * @param attrName	The name of the first attribute within that XML tag.
	 * @param attrValue	The value of the first attribute
	 * @param attr2Name	The name of the related attribute that is required to
	 * find this value.
	 * @return	Double value assigned to this tag.
	 */
	public Double getDblAttrOfChildSuchAttribute(String childName,
						String attrName, String attrValue, String attr2Name)
	{
		return stringToDouble(getChildSuchAttribute(childName, 
				attrName, attrValue).getAttributeValue(attr2Name));
	}
	
	/*************************************************************************
	 * Reading values with units (all Doubles)
	 */
	
	/**
	 * \brief Returns a length parameter from the XML, converting to the
	 * correct unit as required.
	 * 
	 * @param paramName	The name of the parameter for which the length value
	 * should be returned.
	 * @return	The length value assigned to this parameter in the protocol
	 * file.
	 */
	public Double getParamLength(String paramName)
	{
		unit = new StringBuffer("");
		value = getParamDbl(paramName, unit);
		value *= utils.UnitConverter.length(unit.toString());
		return value;
	}
	
	/**
	 * \brief Returns an area parameter from the XML, converting to the
	 * correct unit as required.
	 * 
	 * @param paramName	The name of the parameter for which the area value
	 * should be returned.
	 * @return	The area value assigned to this parameter in the protocol
	 * file.
	 */
	public Double getParamArea(String paramName)
	{
		unit = new StringBuffer("");
		value = getParamDbl(paramName, unit);
		value *= utils.UnitConverter.area(unit.toString());
		return value;
	}

	/**
	 * \brief Returns the mass of a specified parameter from the XML,
	 * calculating this using the unit of measurement for that parameter.
	 * 
	 * @param paramName	The name of the parameter for which the mass should be
	 * returned.
	 * @return	The calculated mass of this parameter.
	 */
	public Double getParamMass(String paramName)
	{
		unit = new StringBuffer("");
		value = getParamDbl(paramName, unit);
		value *= utils.UnitConverter.mass(unit.toString());
		return value;
	}

	/**
	 * \brief Gets a given parameter from the protocol file and converts into
	 * a double representing time.
	 * 
	 * @param paramName	The parameter to be retrieved from the XML file.
	 * @return	Double of the value of this parameter converted into a unit
	 * of time.
	 */
	public Double getParamTime(String paramName)
	{
		unit = new StringBuffer("");
		value = getParamDbl(paramName, unit);
		value *= utils.UnitConverter.time(unit.toString());
		return value;
	}

	/**
	 * \brief Reads a concentration from the protocol file and converts this
	 * to the required unit.
	 * 
	 * @param paramName	The name of the parameter to be retrieved
	 * from the protocol file.
	 * @return	The calculated concentration level for this simulation
	 * parameter.
	 */
	public Double getParamConcn(String paramName)
	{
		unit = new StringBuffer("");
		value = getParamDbl(paramName, unit);
		value *= utils.UnitConverter.mass(unit.toString());
		value *= utils.UnitConverter.volume(unit.toString());
		return value;
	}
	
	/**
	 * \brief Reads a speed from the protocol file and converts this
	 * to the required unit.
	 * 
	 * @param paramName	The name of the parameter to be retrieved
	 * from the protocol file.
	 * @return	The calculated speed for this simulation parameter.
	 */
	public Double getParamSpeed(String paramName)
	{
		unit = new StringBuffer("");
		value = getParamDbl(paramName, unit);
		value *= utils.UnitConverter.length(unit.toString());
		value *= utils.UnitConverter.time(unit.toString());
		return value;
	}
	
	
	public Double getParamDiffusivity(String paramName)
	{
		unit = new StringBuffer("");
		value = getParamDbl(paramName, unit);
		value *= utils.UnitConverter.time(unit.toString());
		value *= UnitConverter.length(unit.toString());
		value *= UnitConverter.length(unit.toString());
		return value;
	}
	
	/*************************************************************************
	 * Creating instances
	 */
		
	/**
	 * \brief Creates an instance of a class using a string containing that
	 * class name.
	 * 
	 * Useful for a set of boundary conditions, for example.
	 * 
	 * @param prefix	The class for which a new instance is being created.
	 * @return	An object of the class stated in the prefix string.
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
			LogFile.writeLogAlways("Unable to create class: "+prefix);
			return null;
		}
	}
}
