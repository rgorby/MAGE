import matplotlib as mpl
import numpy as np
import os

#Pull colormap from text file
__location__ = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))

fIn = os.path.join(__location__,"cmDDiv.txt")
Q = np.loadtxt(fIn,skiprows=1)
cmDiv = mpl.colors.ListedColormap(Q/255.0)

fIn = os.path.join(__location__,"cmMLT.txt")
Q = np.loadtxt(fIn,skiprows=1)
cmMLT = mpl.colors.ListedColormap(Q/255.0)

