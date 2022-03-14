#!/usr/bin/python3

import argparse
import configparser
import os
from string import ascii_letters
import subprocess
import xml.etree.ElementTree as ET
from xml.etree.ElementTree import Element


def convertUnits(myLine):
	# Check if delete
	if ("DEL!" in myLine):
		# Get the setting
		settingString = myLine.split("=")[0]
		return (settingString + "= DEL!")

	# Check if just a comment
	elif ("[" not in myLine):
		return myLine.split("#")[0]

	# Get units
	unitString = myLine.split("#")[1]
	
	# Set multiplier for units
	multiplier = 0
	
	if ("[sec]" in unitString):
		multiplier = 1

	elif ("[min]" in unitString):
		multiplier = 60

	elif ("[hrs]" in unitString):
		multiplier = 3600

	else:
		# Find what is in unit string an tell the user it is incorrect
		incorrectUnit = unitString.split("[")[1]
		incorrectUnit = incorrectUnit.split("]")[0]
		print("ERROR: Incorrect unit type for conversion: " + incorrectUnit)
		exit()

	# Get the setting
	settingString = myLine.split("=")[0]

	# Get the number
	numberList = []
	for t in myLine.split():
		try:
			numberList.append(float(t))
		except:
			pass

	# Multiply
	actualNumber = numberList[0] * multiplier

	# Set this to the number string
	numberString = str(actualNumber)

	# Put the whole thing back together and return
	return (settingString + "= " + numberString)


def initialize(settingsFile):
    """Initialize using the settings file."""

    # Set up whitelist for section names (letters only)
    whiteList = set(ascii_letters)

    # Open up settings file
    settings = settingsFile # First argument is input file name

    # Create a sub-folder called "Parsed Settings", or something
    os.system("mkdir .Settings")

    # Go through the settings file and check for section flags
    with open(settings, 'r') as file:
        content = file.read()

    # Move to subfilder
    os.chdir(".Settings")

    # Put everything between section flags in it's own subfile in the created directory.
    # The name should correspond to the section.
    contentSplit = content.splitlines()

    # Create a temporary string to hold the temporary settings files
    temporary = ""
    name = ''.join(l for l in contentSplit[0] if l in whiteList)
    contentSplit.pop(0)

    pos = -1

    # Loop through the split file
    for line in contentSplit:
        pos += 1
        
        if (len(line) < 1):
            temporary = temporary + "\n"
        # If there is a # at the beginning of the string
        elif (line[0] == '#'):
            # Then write file, change name, and reset temporary
            tempFile = open(name + ".ini", "w")
            tempFile.write(temporary)
            tempFile.close()
            name = ''.join(l for l in line if l in whiteList)
            temporary = ""
        elif (pos == (len(contentSplit) - 1)):
            # Last one, add line then write everything out!
            # Check for converting!
            if ("#" in line):
                temporary = temporary + convertUnits(line) + "\n"

            else:
                temporary = temporary + line + "\n"

            tempFile = open(name + ".ini", "w")
            tempFile.write(temporary)
            tempFile.close()
        # Check for conversion
        elif ("#" in line):
            temporary = temporary + convertUnits(line) + "\n"
        else:
            # Add line to temporary
            temporary = temporary + line + "\n"


# Found on Stack Overflow.
# This indents everything in the elem node properly since apparently etree doesn't do that on it's own...
def indent(elem, level=0):
    i = "\n" + level*"  "
    if len(elem):
        if not elem.text or not elem.text.strip():
            elem.text = i + "  "
        if not elem.tail or not elem.tail.strip():
            elem.tail = i
        for elem in elem:
            indent(elem, level+1)
        if not elem.tail or not elem.tail.strip():
            elem.tail = i
    else:
        if level and (not elem.tail or not elem.tail.strip()):
            elem.tail = i


def create_from_xml_template(ini_file, xml_file, template):
    """Convert the .ini file to a .xml file using a template."""

    # Get settings file name
    settings = ini_file

    # Get output file name
    output = xml_file

    # Run Initialization for Settings
    initialize(settings)

    os.chdir('..')

    # Check if initialize failed/didn't produce output
    settingsFolder = os.listdir(".Settings")
    if (len(settingsFolder) == 0):
        print("Initialization failed. Aborting.")
        # Cleanup Settigns folder
        subprocess.Popen("rm -rf .Settings/", shell=True)
        exit(1)

    # Read in default settings for this template
    tree = ET.parse(template)
    templateRoot = tree.getroot()

    # Make a Parser
    user = configparser.RawConfigParser()

    # Try to turn on case sensitivity
    user.optionxform = lambda option: option

    os.chdir(".Settings")

    inputDicts = {}

    # Iterate through each file there and make the root node the key for the resultant tree in a dictionary
    for filename in os.listdir():
        user = configparser.RawConfigParser()
        user.optionxform = lambda option: option
        user.read(filename)
        temp = filename.split('.')
        newFileName = temp[0]
        inputDicts[newFileName] = user

    #print(inputDicts)

    # Make bigger ETree by adding all roots as sub elements
    top = Element('Kaiju')

    for key in inputDicts.keys():
        temp = Element(key)
        # ET.dump(temp)
        for section in inputDicts[key].sections():
            deeperTemp = Element(section)
            for option in inputDicts[key].options(section):
                deeperTemp.set(option, inputDicts[key].get(section, option))
            ET.SubElement(temp, deeperTemp.tag, deeperTemp.attrib)
        top.append(temp)

    # ET.dump(top)

    # Go through the new settings and see if they match elements in the default
    for child in top:
        #print(child.tag)
        # Try to find that child tag in the default tree
        if (templateRoot.find(child.tag) is not None):
            # If it exists, go one level down and iterate through those nodes
            firstLevel = templateRoot.find(child.tag)
            for lower in child:
                # Find the corresponding tag in the default
                nextLevel = firstLevel.find(lower.tag)
                # Check if that tag exists. If not, just add it
                if  (nextLevel is not None):
                    for item in lower.keys():
                        # Check for the delete flag
                        if ("DEL!" in lower.get(item)):
                            # Check if that option exists.
                            if (nextLevel.get(item) is None):
                                # Just print Debug statements
                                #print(nextLevel)
                                #print(item)
                                #print(lower.get(item))
                                #print(nextLevel.attrib)
                                continue
                            
                            else:
                                del nextLevel.attrib[item]
                        
                        else:
                            nextLevel.set(item, lower.get(item))
                else:
                    # Check for the delete flag
                    if ("DEL!" in lower.attrib):
                        # print(lower.tag)
                        # print(lower.attrib)
                        # Don't add anything
                        continue

                    else:
                        ET.SubElement(firstLevel, lower.tag, lower.attrib)
        else:
            # Else, just add that element to the root
            templateRoot.insert(0, child)
    #	for key in child.keys():
    #		# If tag appears, check sub-entries
    #		if (templateRoot.find(key) is not None):
    #			subelement = templateRoot.find(key)
    #			print("I found " + key + " in the default tree")
    #			# For each option in element, add that to the ETree element
    #			for item in key:
    #				subelement.set(item[0], item[1])
    #		# If tag does not appear, append new one to ETree
    #		else:
    #			print("I did not find " + key + " in the default tree")
    #			tempElement = ET.Element(key)
    #			# Go through the options and add them to a new element
    #			for item in child:
    #				tempElement.set(item[0],item[1])
    #		
    #			# Insert this new element at the end of the current section	
    #			templateRoot.insert(len(list(templateRoot)),tempElement)

    # Run root through the indentation function
    indent(templateRoot)

    os.chdir('..')

    # Write the XML file
    tree.write(output)

    # Cleanup Settigns folder
    subprocess.Popen("rm -rf .Settings/", shell=True)

    print("\n\nXML generation complete!\n\n")


def create_xml_template(ini_file, xml_file):
    """Convert the .ini file to a template .xml file."""

    # Get settings file name
    settings = ini_file

    # Get output file name
    output = xml_file

    # Run Initialization for Settings
    initialize(settings)

    # Make a Parser
    user = configparser.RawConfigParser()

    # Try to turn on case sensitivity
    user.optionxform = lambda option: option

    #print(os.getcwd())

    inputDicts = {}

    # Iterate through each file there and make the root node the key for the resultant tree in a dictionary
    for filename in os.listdir():
        user = configparser.RawConfigParser()
        user.optionxform = lambda option: option
        user.read(filename)
        temp = filename.split('.')
        newFileName = temp[0]
        inputDicts[newFileName] = user

    # Make bigger ETree by adding all roots as sub elements
    top = Element('Kaiju')

    for key in inputDicts.keys():
        temp = Element(key)
        # ET.dump(temp)
        for section in inputDicts[key].sections():
            deeperTemp = Element(section)
            for option in inputDicts[key].options(section):
                if ("DEL!" in inputDicts[key].get(section, option)):
                    continue
                else:
                    deeperTemp.set(option, inputDicts[key].get(section, option))
            ET.SubElement(temp, deeperTemp.tag, deeperTemp.attrib)
        top.append(temp)

    # ET.dump(top)

    # Run root through the indentation function
    indent(top)

    os.chdir('..')

    # Create Etree with the root
    myTree = ET.ElementTree(top)

    # Write the XML file
    myTree.write(output)

    # Cleanup Settigns folder
    subprocess.Popen("rm -rf .Settings/", shell=True)

    print("\n\nTemplate creation complete!\n\n")


if __name__ == "__main__":
    """Convert a .ini file to a .xml file."""

    # Create the command-line argument parser.
    parser = argparse.ArgumentParser(description="Convert Kaiju .ini files to XML configuration or XML configuration template files.")
    parser.add_argument("-t", dest="template", help="Path to .xml file as source template")
    parser.add_argument("ini_file", help="Path to .ini file to convert")
    parser.add_argument("xml_file", help="Path to .xml file output")
    parser.add_argument("-d", "--debug", help="Activate debug mode", action="store_true")
    parser.add_argument("-v", "--verbose", help="Activate verbose execution mode", action="store_true")

    # Parse the commmand-line arguments.
    args = parser.parse_args()

    # If a template file was specified, use it. Otherwise, convert directly
    # from .ini format to .xml format.
    if args.template:
        if args.verbose:
            print("Converting %s to %s using template %s." % (args.ini_file, args.xml_file, args.template))
        create_from_xml_template(args.ini_file, args.xml_file, args.template)
    else:
        if args.verbose:
            print("Converting %s to XML output %s." % (args.ini_file, args.xml_file))
        create_xml_template(args.ini_file, args.xml_file)

