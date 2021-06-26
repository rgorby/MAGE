#!/usr/bin/env python
# Make XMF files from an MPI-decomposed gamera run

import argparse
from argparse import RawTextHelpFormatter
import numpy as np
import kaipy.gamera.block_gampp as gampp
import xml.etree.ElementTree as et
import xml.dom.minidom
import kaipy.kaiH5 as kh5

import os
# Generate name of restart file

def VectorOutput(Vname,Vecs,VDims,n,h5f):
    vAtt = et.SubElement(Grid, "Attribute")
    vAtt.set("Name", Vname)
    vAtt.set("AttributeType", "Vector")
    vAtt.set("Center", "Cell")
    FDI = et.SubElement(vAtt, "DataItem")
    FDI.set("Dimensions", VDims)
    FDI.set("Function", "JOIN( $0, $1, $2)" )
    FDI.set("ItemType", "Function")
    XDI = et.SubElement(FDI, "DataItem")
    XDI.set("Dimensions", cDims)
    XDI.set("NumberType", "Float")
    XDI.set("Precision", "4")
    XDI.set("Format", "HDF")
    XDI.text = "%s:/Step#%d/%s" % (h5F, n, Vecs[0])
    YDI = et.SubElement(FDI, "DataItem")
    YDI.set("Dimensions", cDims)
    YDI.set("NumberType", "Float")
    YDI.set("Precision", "4")
    YDI.set("Format", "HDF")
    YDI.text = "%s:/Step#%d/%s" % (h5F, n, Vecs[1])
    ZDI = et.SubElement(FDI, "DataItem")
    ZDI.set("Dimensions", cDims)
    ZDI.set("NumberType", "Float")
    ZDI.set("Precision", "4")
    ZDI.set("Format", "HDF")
    ZDI.text = "%s:/Step#%d/%s" % (h5F, n, Vecs[2])


def genName(bStr, i, j, k, Ri, Rj, Rk,isOld):
    n = j + i*Rj + k*Ri*Rj
    if (isOld):
        fID = bStr + \
        "_%04d_%04d_%04d_%04d_%04d_%04d_%012d.h5" % (Ri, Rj, Rk, i, j, k,n)
    else:
        fID = bStr + \
        "_%04d_%04d_%04d_%04d_%04d_%04d.gam.h5" % (Ri, Rj, Rk, i, j, k)
    return fID


if __name__ == "__main__":
    # Defaults
    fdir = os.getcwd()
    ftag = "msphere"
    outid = "sim"
    MainS = """Creates series of XMF files from MPI-decomposed Gamera run
	"""

    parser = argparse.ArgumentParser(
        description=MainS, formatter_class=RawTextHelpFormatter)
    parser.add_argument('-d', type=str, metavar="directory", default=fdir,
                        help="Directory to read from (default: %(default)s)")
    parser.add_argument('-outd', type=str, metavar="directory", default=fdir,
                        help="Directory to output (default: %(default)s)")
    parser.add_argument('-id', type=str, metavar="runid",
                        default=ftag, help="RunID of data (default: %(default)s)")
    parser.add_argument('-outid', type=str, metavar="outid", default=outid,
                        help="RunID of output XMF files (default: %(default)s)")

    # Finalize parsing
    args = parser.parse_args()
    fdir = args.d
    fodir = args.outd
    ftag = args.id
    outid = args.outid

    # ---------------------
    # Init data
    print(fdir,"ftag",ftag)
    gamData = gampp.GameraPipe(fdir, ftag)
    print("isOld  = ", gamData.isOld)

    print("number of variable",gamData.Nv)
    print("Vars ",gamData.vIDs)

    # ---------------------
    # Do work
    Ri = gamData.Ri
    Rj = gamData.Rj
    Rk = gamData.Rk

    tOut = 0

    topoStr = "3DSMesh"
    geoStr = "X_Y_Z"

#    vIDs = ["D", "Vx", "Vy", "Vz", "P", "Bx", "By", "Bz"]  # ,"Jx","Jy","Jz"]
    vIDs = ["D", "P"]
    VVec =["Vx", "Vy", "Vz"] 
    BVec =["Bx", "By", "Bz"] 
    JVec =["Jx", "Jy", "Jz"]
    EVec =["Ex", "Ey", "Ez"]

    haveV = False
    haveB = False
    haveJ = False
    haveE = False
    for i in gamData.vIDs:
        if ( i == "Vx" ): haveV = True
        if ( i == "Bx" ): haveB = True
        if ( i == "Jx" ): haveJ = True
        if ( i == "Ex" ): haveE = True

    Nv = len(vIDs)
    for n in range(gamData.s0, gamData.sFin+1):
        nslc = n-gamData.s0
        print(n,gamData.T[nslc])

        # Create XMF tree
        Xdmf = et.Element("Xdmf")
        Xdmf.set("Version", "2.0")

        Dom = et.SubElement(Xdmf, "Domain")

        # Spatial collection
        meshName = "g%02d_mesh"%(nslc)
        gCol = et.SubElement(Dom, "Grid")
        gCol.set("Name", meshName)
        gCol.set("GridType", "Collection")
        gCol.set("CollectionType", "Spatial")

        Time = et.SubElement(gCol, "Time")
        Time.set("Value", "%s" % (str(gamData.T[nslc])))

        Cyc = et.SubElement(gCol, "Cycle")
        Cyc.set("Value", "%d" % (n))

        # Now loop over MPI decomposition
        for i in range(Ri):
            for j in range(Rj):
                for k in range(Rk):
                    nMPI = j + i*Rj + k*Ri*Rj
                    h5F = fdir + '/' + genName(ftag, i, j, k, Ri, Rj, Rk,gamData.isOld)


                    Ni = gamData.dNi[i]
                    Nj = gamData.dNj[j]
                    Nk = gamData.dNk[k]
                    Ndim = 3

                    iDims = "%d %d %d" % (Nk+1, Nj+1, Ni+1)
                    VDims = "%d %d %d %d" % (Nk+0, Nj+0, Ni+0, Ndim )
                    cDims = "%d %d %d" % (Nk+0, Nj+0, Ni+0)
                    iDimA = "%d %d %d" % (Nk*Rk+1, Nj*Rj+1, Ni*Ri+1)
                    
                    # Create new subgrid
                    gName = meshName+"%d" % (nMPI)
                    Grid = et.SubElement(gCol, "Grid")
                    Grid.set("GridType", "Uniform")
                    Grid.set("Name", gName)
                    print(Grid)

                    # Time = et.SubElement(Grid,"Time")
                    # Time.set("TimeType","Single")
                    # Time.set("Value","%s"%(str(gamData.T[nslc])))

                    Topo = et.SubElement(Grid, "Topology")
                    Topo.set("TopologyType", topoStr)
                    Topo.set("NumberOfElements", iDims)
                    Geom = et.SubElement(Grid, "Geometry")
                    Geom.set("GeometryType", geoStr)

                    xC = et.SubElement(Geom, "DataItem")
                    xC.set("Dimensions", iDims)
                    xC.set("NumberType", "Float")
                    xC.set("Precision", "4")
                    xC.set("Format", "HDF")
                    xC.text = "%s:/X" % (h5F)

                    yC = et.SubElement(Geom, "DataItem")
                    yC.set("Dimensions", iDims)
                    yC.set("NumberType", "Float")
                    yC.set("Precision", "4")
                    yC.set("Format", "HDF")
                    yC.text = "%s:/Y" % (h5F)

                    zC = et.SubElement(Geom, "DataItem")
                    zC.set("Dimensions", iDims)
                    zC.set("NumberType", "Float")
                    zC.set("Precision", "4")
                    zC.set("Format", "HDF")
                    zC.text = "%s:/Z" % (h5F)

                    # Create variables
                    for v in range(Nv):
                        vID = vIDs[v]
                        vAtt = et.SubElement(Grid, "Attribute")
                        vAtt.set("Name", vID)
                        vAtt.set("AttributeType", "Scalar")
                        vAtt.set("Center", "Cell")
                        aDI = et.SubElement(vAtt, "DataItem")
                        aDI.set("Dimensions", cDims)
                        aDI.set("NumberType", "Float")
                        aDI.set("Precision", "4")
                        aDI.set("Format", "HDF")
                        aDI.text = "%s:/Step#%d/%s" % ( h5F, n, vID)
                        
                     # create vectors
                    if (haveV): VectorOutput("V",VVec,VDims,n,h5F)
                    if (haveB): VectorOutput("B",BVec,VDims,n,h5F)
                    if (haveJ): VectorOutput("J",JVec,VDims,n,h5F)
                    if (haveE): VectorOutput("E",EVec,VDims,n,h5F)

        # Write output
        fOut = "%s/%s.%06d.xmf" % (fodir, outid, tOut)
        print("Writing %s" % (fOut))
 #       xTree = et.ElementTree(Xdmf)
 #       xTree.write(fOut, pretty_print=True,
 #                   xml_declaration=True, encoding='UTF-8')
        xmlStr = xml.dom.minidom.parseString(et.tostring(Xdmf)).toprettyxml(indent="    ")
        with open(fOut,"w") as f:
            f.write(xmlStr)


        # Prep for next step
        tOut = tOut+1
