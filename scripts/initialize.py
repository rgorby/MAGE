import sys
import os
from string import ascii_letters

def convertUnits(myLine):
	# Check if just a comment
	if ("[" not in myLine):
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
		print("ERROR: Incorrect unit type for conversion")
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

# Check number of arguments
if (len(sys.argv) < 2):
	print("ERROR: Too few arguments")
	exit()
elif (len(sys.argv) > 2):
	print("ERROR: Too many arguments")
	exit()

# Set up whitelist for section names (letters only)
whiteList = set(ascii_letters)

# Open up settings file
settings = sys.argv[1] # First argument is input file name

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
	
print(contentSplit)

