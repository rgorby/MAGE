#!/bin/env python

"""Convert a MAGE .ini file to a MAGE .XML file.

Convert a MAGE .ini file to a MAGE .XML file.

Authors
-------
Brent Smith
Eric Winter
"""

# Include standard modules.
import argparse
import configparser
import os
from string import ascii_letters
import subprocess
import xml.etree.ElementTree as ET
from xml.etree.ElementTree import Element

# Include 3rd-party modules.

# Include project modules.


# Program constants.

# Program description string.
description = """Convert a MAGE .ini file to a MAGE .XML file."""


def create_command_line_parser():
    """Create the command-line argument parser.

    Create the command-line parser.

    Parameters
    ----------
    None

    Returns
    -------
    parser : argparse.ArgumentParser
        Parser for command-line arguments.
    """
    parser = argparse.ArgumentParser(
        description=description, formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument(
        "-d", "--debug", action="store_true", default=False,
        help="Print debugging output (default: %(default)s)."
    )
    parser.add_argument(
        "-v", "--verbose", action="store_true", default=False,
        help="Print verbose output (default: %(default)s)."
    )
    parser.add_argument(
        "ini_path",
        help="Path to .ini file to convert."
    )
    parser.add_argument(
        "xml_path",
        help="Path to .xml file to create."
    )
    return parser


def convertUnits(myLine):
    """Convert units."""
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
    with open(settings, "r") as file:
        content = file.read()

    # Move to subfilder
    os.chdir(".Settings")

    # Put everything between section flags in it's own subfile in the created directory.
    # The name should correspond to the section.
    contentSplit = content.splitlines()

    # Create a temporary string to hold the temporary settings files
    temporary = ""
    name = "".join(l for l in contentSplit[0] if l in whiteList)
    contentSplit.pop(0)  # Remove leading comment line.

    # Loop through the split file
    pos = -1
    for line in contentSplit:
        pos += 1
        if len(line) < 1:
            # Blank lines
            temporary = temporary + "\n"
        elif (line[0] == "#"):
            # If there is a # at the beginning of the string
            # Then write file, change name, and reset temporary
            tempFile = open(name + ".ini", "w")
            tempFile.write(temporary)
            tempFile.close()
            name = "".join(l for l in line if l in whiteList)
            temporary = ""
        elif pos == len(contentSplit) - 1:
            # Last one, add line then write everything out!
            # Check for converting!
            if ("#" in line):
                temporary = temporary + convertUnits(line) + "\n"
            else:
                temporary = temporary + line + "\n"
            tempFile = open(name + ".ini", "w")
            tempFile.write(temporary)
            tempFile.close()
        elif ("#" in line):
            # Check for conversion
            temporary = temporary + convertUnits(line) + "\n"
        else:
            # Add line to temporary
            temporary = temporary + line + "\n"


def indent(elem, level=0):
    """Found on Stack Overflow.

    This indents everything in the elem node properly since apparently etree
    doesn't do that on it's own...
    """
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


def create_xml_template(ini_path, xml_path):
    """Convert a MAGE .ini file to .xml format.

    Convert a MAGE .ini file to .xml format.

    Parameters
    ----------
    ini_path : str
        Path to .ini file to convert.
    xml_path : str
        Path to .xml file to create.
    
    Returns
    -------
    None
    """
    # Run Initialization for Settings
    settings = ini_path

    # Get output file name
    output = xml_path

    # Run Initialization for Settings
    initialize(settings)

    # Make a Parser
    user = configparser.RawConfigParser()

    # Try to turn on case sensitivity
    user.optionxform = lambda option: option

    inputDicts = {}

    # Iterate through each file there and make the root node the key for the
    # resultant tree in a dictionary.
    for filename in os.listdir():
        user = configparser.RawConfigParser()
        user.optionxform = lambda option: option
        user.read(filename)
        temp = filename.split(".")
        newFileName = temp[0]
        inputDicts[newFileName] = user

    # Make bigger ETree by adding all roots as sub elements
    top = Element('Kaiju')

    for key in inputDicts.keys():
        temp = Element(key)
        for section in inputDicts[key].sections():
            deeperTemp = Element(section)
            for option in inputDicts[key].options(section):
                if ("DEL!" in inputDicts[key].get(section, option)):
                    continue
                else:
                    deeperTemp.set(option, inputDicts[key].get(section, option))
            ET.SubElement(temp, deeperTemp.tag, deeperTemp.attrib)
        top.append(temp)

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


# def create_from_xml_template(ini_file, xml_file, template):
#     """Convert the .ini file to a .xml file using a template."""

#     # Get settings file name
#     settings = ini_file

#     # Get output file name
#     output = xml_file

#     # Run Initialization for Settings
#     initialize(settings)

#     os.chdir('..')

#     # Check if initialize failed/didn't produce output
#     settingsFolder = os.listdir(".Settings")
#     if (len(settingsFolder) == 0):
#         print("Initialization failed. Aborting.")
#         # Cleanup Settigns folder
#         subprocess.Popen("rm -rf .Settings/", shell=True)
#         exit(1)

#     # Read in default settings for this template
#     tree = ET.parse(template)
#     templateRoot = tree.getroot()

#     # Make a Parser
#     user = configparser.RawConfigParser()

#     # Try to turn on case sensitivity
#     user.optionxform = lambda option: option

#     os.chdir(".Settings")

#     inputDicts = {}

#     # Iterate through each file there and make the root node the key for the resultant tree in a dictionary
#     for filename in os.listdir():
#         user = configparser.RawConfigParser()
#         user.optionxform = lambda option: option
#         user.read(filename)
#         temp = filename.split('.')
#         newFileName = temp[0]
#         inputDicts[newFileName] = user

#     #print(inputDicts)

#     # Make bigger ETree by adding all roots as sub elements
#     top = Element('Kaiju')

#     for key in inputDicts.keys():
#         temp = Element(key)
#         # ET.dump(temp)
#         for section in inputDicts[key].sections():
#             deeperTemp = Element(section)
#             for option in inputDicts[key].options(section):
#                 deeperTemp.set(option, inputDicts[key].get(section, option))
#             ET.SubElement(temp, deeperTemp.tag, deeperTemp.attrib)
#         top.append(temp)

#     # ET.dump(top)

#     # Go through the new settings and see if they match elements in the default
#     for child in top:
#         #print(child.tag)
#         # Try to find that child tag in the default tree
#         if (templateRoot.find(child.tag) is not None):
#             # If it exists, go one level down and iterate through those nodes
#             firstLevel = templateRoot.find(child.tag)
#             for lower in child:
#                 # Find the corresponding tag in the default
#                 nextLevel = firstLevel.find(lower.tag)
#                 # Check if that tag exists. If not, just add it
#                 if  (nextLevel is not None):
#                     for item in lower.keys():
#                         # Check for the delete flag
#                         if ("DEL!" in lower.get(item)):
#                             # Check if that option exists.
#                             if (nextLevel.get(item) is None):
#                                 # Just print Debug statements
#                                 #print(nextLevel)
#                                 #print(item)
#                                 #print(lower.get(item))
#                                 #print(nextLevel.attrib)
#                                 continue
                            
#                             else:
#                                 del nextLevel.attrib[item]
                        
#                         else:
#                             nextLevel.set(item, lower.get(item))
#                 else:
#                     # Check for the delete flag
#                     if ("DEL!" in lower.attrib):
#                         # print(lower.tag)
#                         # print(lower.attrib)
#                         # Don't add anything
#                         continue

#                     else:
#                         ET.SubElement(firstLevel, lower.tag, lower.attrib)
#         else:
#             # Else, just add that element to the root
#             templateRoot.insert(0, child)
#     #	for key in child.keys():
#     #		# If tag appears, check sub-entries
#     #		if (templateRoot.find(key) is not None):
#     #			subelement = templateRoot.find(key)
#     #			print("I found " + key + " in the default tree")
#     #			# For each option in element, add that to the ETree element
#     #			for item in key:
#     #				subelement.set(item[0], item[1])
#     #		# If tag does not appear, append new one to ETree
#     #		else:
#     #			print("I did not find " + key + " in the default tree")
#     #			tempElement = ET.Element(key)
#     #			# Go through the options and add them to a new element
#     #			for item in child:
#     #				tempElement.set(item[0],item[1])
#     #		
#     #			# Insert this new element at the end of the current section	
#     #			templateRoot.insert(len(list(templateRoot)),tempElement)

#     # Run root through the indentation function
#     indent(templateRoot)

#     os.chdir('..')

#     # Write the XML file
#     tree.write(output)

#     # Cleanup Settigns folder
#     subprocess.Popen("rm -rf .Settings/", shell=True)



def main():
    """Main program.

    This is the main program routine.

    Parameters
    ----------
    None

    Returns
    -------
    None
    """

    # Set up the command-line parser.
    parser = create_command_line_parser()

    # Parse the commmand-line arguments.
    args = parser.parse_args()
    debug = args.debug
    verbose = args.verbose
    ini_path = args.ini_path
    xml_path = args.xml_path
    if debug:
        print("args = %s" % args)

    # Convert the .ini to .xml format.
    if verbose:
        print("Converting %s to %s." % (ini_path, xml_path))
    create_xml_template(args.ini_path, args.xml_path)


if __name__ == "__main__":
    """Begin main program."""
    main()
