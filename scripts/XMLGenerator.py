import configparser
import xml.etree.ElementTree as ET

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

# Read in default settings
tree = ET.parse('test.xml')
root = tree.getroot()

# Initilization
user = configparser.RawConfigParser()

# Try to turn on case sensitivity
user.optionxform = lambda option: option

# Read user-defined options in
user.read('config.ini')

# Get all sections of user config file
sections = user.sections()

# Go through the new settings and see if they match elements in the default
for element in sections:
	# If tag appears, check sub-entries
	if (root.find(element) is not None):
		subelement = root.find(element)
		print("I found " + element + " in the default tree")
		# For each option in element, add that to the ETree element
		for item in user.items(element):
			subelement.set(item[0], item[1])
	# If tag does not appear, append new one to ETree
	else:
		print("I did not find " + element + " in the default tree")
		tempElement = ET.Element(element)
		# Go through the options and add them to a new element
		item = user.items(element)
		for item in user.items(element):
			tempElement.set(item[0],item[1])
		
		# Insert this new element at the end of the current section	
		root.insert(len(list(root)),tempElement)

# Run root through the indentation function
indent(root)

# Write the XML file
tree.write('output.xml')
