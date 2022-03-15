import matplotlib as mpl
import numpy as np
import os

#Pull colormap from text file
__location__ = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))
fIn = os.path.join(__location__,"cmDDiv.txt")
Q = np.loadtxt(fIn)
cmDiv = mpl.colors.ListedColormap(Q/255.0)

# Q = np.loadtxt("cmMLT.txt")
# cmMLT = mpl.colors.ListedColormap(Q/255.0)

