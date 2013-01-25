/**
 * Project iDynoMiCS (copyright -> see Idynomics.java) 
 *__________________________________________________________________________
 * XML parser
 * 
 */

/**
 * @since June 2006
 * @version 1.0
 * @author Andreas D�tsch (andreas.doetsch@helmholtz-hzi.de), Helmholtz Centre for Infection Research (Germany)
 * @author Laurent Lardon (lardonl@supagro.inra.fr), INRA, France
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

public class XMLParser implements Serializable {

	// Serial version used for the serialisation of the class
	private static final long   serialVersionUID = 1L;
	private Document            document;
	private Element             _localRoot;
	private static double       value;
	private static StringBuffer unit;

	// protected LinkedList<Element> markUpList;
	// private Element currentElement;

	/* ______________________ CONSTRUCTOR ___________________________________ */
	public XMLParser(Element localRoot) {
		_localRoot = localRoot;
	}

	public XMLParser(String fileName) {
		openXMLDocument(fileName, false);
		_localRoot = document.getRootElement();
	}

	/**
	 * Now create the DOM document associated to the XML file
	 * @param fileName
	 * @param testDTD
	 */
	public void openXMLDocument(String fileName, boolean testDTD) {
		try {
			document = (new SAXBuilder(testDTD)).build(new File(fileName));
		} catch (Exception e) {
			LogFile.writeLog("XMLParser.openXMLdocument(): Initialisation of the XML parser failed");
			LogFile.writeLog("Error: "+e);
			//LogFile.writeLog("File does not exist: "+fileName);
			System.exit(-1);
		}
	}

	@SuppressWarnings("unchecked")
	public LinkedList<Element> buildSetMarkUp(String childMarkup) {
		List<Element> childList = _localRoot.getChildren(childMarkup);
		LinkedList<Element> out = new LinkedList<Element>();
		for (Object anElement : childList) {
			out.add((Element) anElement);
		}
		return out;
	}

	@SuppressWarnings("unchecked")
	public LinkedList<XMLParser> buildSetParser(String childMarkup) {
		List<Element> childList = _localRoot.getChildren(childMarkup);
		LinkedList<XMLParser> out = new LinkedList<XMLParser>();
		for (Object anElement : childList) {
			out.add(new XMLParser((Element) anElement));
		}
		return out;
	}

	/* ___________________ TO READ PARAMETERS MARK-UPS _____________________ */

	public XMLParser getChild(String childName) {
		return new XMLParser(getChildElement(childName));
	}

	@SuppressWarnings("unchecked")
	public String getParam(String paramName) {
		List<Element> childList = _localRoot.getChildren("param");
		Element aParam;
		for (Object aChild : childList) {
			aParam = (Element) aChild;
			if (aParam.getAttributeValue("name").equals(paramName)) { return aParam.getText(); }
		}
		return null;
	}
	
	@SuppressWarnings("unchecked")
	public Element getParamMarkUp(String paramName) {
		List<Element> childList = _localRoot.getChildren("param");
		Element aParam;
		for (Object aChild : childList) {
			aParam = (Element) aChild;
			if (aParam.getAttributeValue("name").equals(paramName)) { return aParam; }
		}
		return null;
	}
	
	@SuppressWarnings("unchecked")
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

	public double getParamDbl(String paramName) {
		if (getParam(paramName)==null) {
			return Double.NaN;
		} else {
			return Double.parseDouble(getParam(paramName));
		}
	}

	public double getParamDbl(String paramName, StringBuffer unit) {
		if (getParam(paramName)==null) {
			return Double.NaN;
		} else {
			return Double.parseDouble(getParam(paramName, unit));
		}
	}

	public int getParamInt(String paramName) {
		if (getParam(paramName)==null) {
			return 0;
		} else {
			return Integer.parseInt(getParam(paramName));
		}
	}

	public int getParamInt(String paramName, StringBuffer unit) {
		if (getParam(paramName)==null) {
			return 0;
		} else {
			return Integer.parseInt(getParam(paramName, unit));
		}
	}

	public double getParamLength(String paramName) {
		unit = new StringBuffer("");
		value = getParamDbl(paramName, unit);
		value *= utils.UnitConverter.length(unit.toString());
		return value;
	}

	public double getParamMass(String paramName) {
		unit = new StringBuffer("");
		value = getParamDbl(paramName, unit);
		value *= utils.UnitConverter.mass(unit.toString());
		return value;
	}

	public double getParamTime(String paramName) {
		unit = new StringBuffer("");
		value = getParamDbl(paramName, unit);
		value *= utils.UnitConverter.time(unit.toString());
		return value;
	}

	/**
	 * read a concentration and converts the unit
	 * @param paramName
	 * @return
	 */
	public double getParamConc(String paramName) {
		unit = new StringBuffer("");
		value = getParamDbl(paramName, unit);
		value *= utils.UnitConverter.mass(unit.toString());
		value *= utils.UnitConverter.volume(unit.toString());
		return value;
	}

	public DiscreteVector getParamXYZ(String paramName){
		//Element param = getParamMarkUp(paramName);
		return new DiscreteVector(getParamMarkUp(paramName));
	}
	
	public Boolean getParamBool(String paramName) {
		return Boolean.parseBoolean(getParam(paramName));
	}

	public Boolean getParamBool(String paramName, StringBuffer unit) {
		return Boolean.parseBoolean(getParam(paramName, unit));
	}

	@SuppressWarnings("unchecked")
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

	public double getParamSuchDbl(String paramName, String detailName) {
		if (getParamSuch(paramName, detailName)==null) {
			return Double.NaN;
		} else {
			return Double.parseDouble(getParamSuch(paramName, detailName));
		}
	}

	public Boolean getParamSuchBool(String paramName, String detailName) {
		if (getParamSuch(paramName, detailName)==null) {
			return null;
		} else {
			return Boolean.parseBoolean(getParamSuch(paramName, detailName));
		}
	}

	public double getParamSuchDbl(String paramName, String detailName, StringBuffer unit) {
		if (getParamSuch(paramName, detailName)==null) {
			return Double.NaN;
		} else {
			return Double.parseDouble(getParamSuch(paramName, detailName, unit));
		}
	}

	public Boolean getParamSuchBool(String paramName, String detailName, StringBuffer unit) {
		if (getParamSuch(paramName, detailName)==null) {
			return null;
		} else {
			return Boolean.parseBoolean(getParamSuch(paramName, detailName, unit));
		}
	}

	/* ______________________________ _______________________________ */

	/**
	 * Search a child mark-up called childName and return its first child whose
	 * the attribute attrName equals attrValue If no such child is found,
	 * returns the localRoot
	 * @param childName
	 * @param attrName
	 * @param attrValue
	 * @return
	 */
	@SuppressWarnings("unchecked")
	public Element getChildSuchAttribute(String childName, String attrName, String attrValue) {
		List<Element> childList = _localRoot.getChildren(childName);
		for (Object aChild : childList) {
			if (((Element) aChild).getAttributeValue(attrName).equals(attrValue)) { return (Element) aChild; }
		}
		return _localRoot;
	}

	@SuppressWarnings("unchecked")
	public Element getChildSuchAttribute(String childName, String attrName[], String attrValue[]) {
		List<Element> childList = _localRoot.getChildren(childName);
		Element myChild;
		boolean isOk;

		for (Object aChild : childList) {
			myChild = (Element) aChild;
			isOk = true;
			for (int iAtt = 0; iAtt<attrName.length; iAtt++) {
				if (!(myChild.getAttributeValue(attrName[iAtt]).equals(attrValue[iAtt]))) {
					isOk &= false;
				}
			}
			if (isOk) { return myChild; }
		}
		return _localRoot;
	}

	public double getDblAttrOfChildSuchAttribute(String childName, String[] attrName,
	        String[] attrValue, String attr2Name) {
		String out = getChildSuchAttribute(childName, attrName, attrValue).getAttributeValue(
		        attr2Name);
		return (double) Double.parseDouble(out);
	}

	public String getStrAttrOfChildSuchAttribute(String childName, String attrName,
	        String attrValue, String attr2Name) {
		String out = getChildSuchAttribute(childName, attrName, attrValue).getAttributeValue(
		        attr2Name);
		return out;
	}

	/* Read Attribute of the root____________________________________________ */
	public String getAttributeStr(String attributeName) {
		return _localRoot.getAttributeValue(attributeName);
	}

	public double getAttributeDbl(String attributeName) {
		return Double.parseDouble(getAttributeStr(attributeName));
	}

	public String getChildAttrStr(String childName, String attrName) {
		try {
			return _localRoot.getChild(childName).getAttributeValue(attrName);
		} catch (Exception e) {
			return null;
		}
	}

	public double getChildAttrDbl(String childName, String attrName) {
		String value = getChildAttrStr(childName, attrName);
		if (value==null) return 0;
		else return Double.parseDouble(value);
	}

	@SuppressWarnings("unchecked")
	public List<Element> getChildren(String childName) {
		return (List<Element>) _localRoot.getChildren(childName);
	}

	/* */
	public String getChildValueStr(String childName) {
		return _localRoot.getChild(childName).getText();
	}

	public double getChildValueDbl(String childName) {
		return Double.parseDouble(getChildValueStr(childName));
	}

	public String getValue(String attributeName) {
		return _localRoot.getText();
	}

	public Element getChildElement(String childName) {
		return _localRoot.getChild(childName);
	}

	public String getAttribute(String attributeName) {
		return _localRoot.getAttributeValue(attributeName);
	}

	public Element getElement() {
		return _localRoot;
	}

	public Object instanceCreator(String prefix) {
		prefix += "."+this.getAttribute("class");
		try {
			return Class.forName(prefix).newInstance();
		} catch (Exception e) {
			LogFile.writeLog("Unable to create class");
			return null;
		}
	}

	/* ______________ NOT USED ANYMORE ____________________________ */

	@SuppressWarnings("unchecked")
	public double[] getChildListValueDbl(String childName) {
		LinkedList<Element> childList = (LinkedList<Element>) _localRoot.getChildren(childName);
		ListIterator<Element> iter = childList.listIterator();
		double out[] = new double[childList.size()];
		while (iter.hasNext()) {
			out[iter.nextIndex()] = Double.parseDouble(((Element) iter.next()).getText());
		}
		return out;
	}

	/* An array of children value (several child markup)______________________ */
	@SuppressWarnings("unchecked")
	public String[] getChildListValueStr(String childName) {
		LinkedList<Element> childList = (LinkedList<Element>) _localRoot.getChildren(childName);
		LinkedList<String> out = new LinkedList<String>();
		for (Element anElement : childList) {
			out.add(anElement.getText());
		}
		return (String[]) out.toArray(new String[0]);
	}

	public double getDblAttrOfChildSuchAttribute(String childName, String attrName,
	        String attrValue, String attr2Name) {
		String out = getChildSuchAttribute(childName, attrName, attrValue).getAttributeValue(
		        attr2Name);
		return (double) Double.parseDouble(out);
	}

	/* An array of attributes (several child markup, same attribute name______ */
	@SuppressWarnings("unchecked")
	public String[] getChildListAttrStr(String childName, String attrName) {
		List<Element> childList = _localRoot.getChildren(childName);
		LinkedList<String> out = new LinkedList<String>();

		for (Object anElement : childList) {
			out.add(((Element) anElement).getAttributeValue(attrName));
		}
		return (String[]) out.toArray(new String[0]);
	}

	@SuppressWarnings("unchecked")
	public double[] getChildListAttrDbl(String childName, String attrName) {
		List<Element> childList = _localRoot.getChildren(childName);
		ListIterator<Element> iter = childList.listIterator();
		double out[] = new double[childList.size()];
		while (iter.hasNext()) {
			out[iter.nextIndex()] = Double.parseDouble(((Element) iter.next())
			        .getAttributeValue(attrName));
		}
		return out;
	}
}
