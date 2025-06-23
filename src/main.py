import pandas as pd
import re
import os
import time
import tempfile

from . import privateer_core as pvt


def privateer_validation_wrapper(session,OutputFolderPath,m,i,display=False):
    timestr = time.strftime("%Y%m%d-%H%M%S")
    tempdirpath = tempfile.gettempdir()
    InputStructureFilePath = os.path.join(tempdirpath,f"temp_{timestr}.pdb")
    from chimerax.pdb.pdb import save_pdb
    save_pdb(session,InputStructureFilePath,models=[m])
    dpath = os.path.dirname(os.path.abspath(__file__))
    zscorefilepath = os.path.join(dpath,"data","linkage_torsions","privateer_torsions_z_score_database.json")
    glycomicsfilepath = os.path.join(dpath,"data","glycomics","privateer_glycomics_database.json")
    AllGlycans = pvt.validate(InputStructureFilePath,zscorefilepath,glycomicsfilepath)
    os.remove(InputStructureFilePath)
    if display:
        return AllGlycans
    else:
        for glycan in AllGlycans:
            svgstring = glycan["svg"]
            rootID = glycan["RootID"]
            rootID = rootID.replace("/","-")
            svgfile = open(os.path.join(OutputFolderPath,f"model-{i}_{rootID}.svg"),"w")
            svgfile.write(svgstring)
            svgfile.close()
        df = pd.DataFrame.from_dict(AllGlycans)
        df = df.drop(["svg"], axis=1)
        csv_out = os.path.join(OutputFolderPath,f"model-{i}_privateer-report.csv")
        df.to_csv(csv_out)  

_ATOM_FMT = ("ATOM  %5d %-4s%1s"                # serial, atom name, altloc
             "%-3s %1s%4s%1s   "                # res name, chain, seq, insert
             "%8.3f%8.3f%8.3f%6.2f%6.2f      "  # xyz, occupancy, bfactor
             "%4s%2s%2s")                       # segment, element, charge

def draw_glycoblocks(session,Glycans, modelID):
    from chimerax.surface.shapes import cylinder_geometry, box_geometry, sphere_geometry
    from chimerax.core.models import Drawing, Model
    from chimerax.geometry import Place, vector_rotation, scale
    import numpy as np
    dm = Model(f'model {modelID} glycoblocks',session)
    blue = [0,144,188,255]
    green = [0,166,81,255]
    red = [237,28,36,255]
    orange = [255,127,0,255]
    yellow = [255,221,0,255]
    grey = [128,128,128,255]
    cyan = [143,204,233,255]
    purple = [165,67,153,255]
    # Glc       = GLC,BGC       = blue          circle
    # Gal       = GAL,GLA       = yellow        circle
    # Man       = MAN,BMA       = green         circle
    # Fuc       = FUC,FCB,FUL   = red           triangle
    # Xyl       = XYS,XYP       = orange        star
    # GlcNAc    = NAG,NDG       = blue          square
    # GalNAc    = NGA,A2G       = yellow        square
    # ManNAc    = BM3,BM7       = green         square
    # GlcN      = GCS,PA1       = half-blue     square
    # GalN      =               = half-yellow   square
    # ManN      =               = half-green    square
    # GlcA      = BDP,GCU       = blue-up       diamond
    # GalA      = GTR,ADA       = yellow-left   diamond
    # ManA      = MAV,BEM       = green-right   diamond
    # Neu5Gc    = NGC, NGE      = cyan          diamond
    # Neu5Ac    = SIA,SLB       = purple        diamond
    # IdoA      = IDR           = tan-down      diamond
    # KDN       = KDM,KDN       = green         diamond
    bluesugars = ["GLC","BGC","NAG","NDG"]
    greensugars = ["MAN","BMA","BM3","BM7","KDM","KDN"]
    redsugars = ["FUC","FCB","FUL"]
    yellowsugars = ["GAL","GLA","NGA","AG2"]
    orangesugars = ["XYS","XYP"]
    cyansugars = ["NGC","NGE"]
    purplesugars = ["SIA","SLB"]

    circlesugars = ["GLC","BGC","GAL","GLA","MAN","BMA"]
    squaresugars = ["NAG","NDG","NGA","A2G","BM3","BM7"]
    trianglesugars = ["FUC","FCB","FUL"]
    starsugars = ["XYL","XYP"]
    diamondsugars = ["BDP","GCU","GTR","ADA","MAV","BEM","NGC","NGE","SIA","SLB","IDR","KDM","KDN"]
    for glycan in Glycans:
        for sugar in glycan["Sugars"]:
            sugarname = sugar["sugarname"]
            sugar_centre_x = sugar["sugar_centre_x"]
            sugar_centre_y = sugar["sugar_centre_y"]
            sugar_centre_z = sugar["sugar_centre_z"]
            sugar_plane_i = sugar["sugar_plane_i"]
            sugar_plane_j = sugar["sugar_plane_j"]
            sugar_plane_k = sugar["sugar_plane_k"]
            sugar_plane_i = sugar_plane_i / (sugar_plane_i**2 + sugar_plane_j**2 + sugar_plane_k**2)
            sugar_plane_j = sugar_plane_j / (sugar_plane_i**2 + sugar_plane_j**2 + sugar_plane_k**2)
            sugar_plane_k = sugar_plane_k / (sugar_plane_i**2 + sugar_plane_j**2 + sugar_plane_k**2)
            C1 = [sugar["C1_x"],sugar["C1_y"],sugar["C1_z"]]
            C2 = [sugar["C2_x"],sugar["C2_y"],sugar["C2_z"]]
            C3 = [sugar["C3_x"],sugar["C3_y"],sugar["C3_z"]]
            C4 = [sugar["C4_x"],sugar["C4_y"],sugar["C4_z"]]
            C5 = [sugar["C5_x"],sugar["C5_y"],sugar["C5_z"]]
            d = Drawing('sugar')
            if sugarname in circlesugars:
                # Create the shape
                nc = 25
                v, n, t = cylinder_geometry(radius = 1.4, height = 1.0, nc=nc)
                # Rotate so we know which points are closest to C1 and C4
                v1 = v[0,:]
                v4 = v[int(nc/2),:]
                tr = vector_rotation((v1[0]-v4[0],v1[1]-v4[1],v1[2]-v4[2]),(C1[0]-C4[0],C1[1]-C4[1],C1[2]-C4[2]))
                v = tr.transform_points(v, in_place=True)
                n = tr.transform_vectors(n, in_place=True)
                # Rotate to match the plane of the sugar ring
                tr = vector_rotation((n[-1,0],n[-1,1],n[-1,2]),(sugar_plane_i, sugar_plane_j, sugar_plane_k))
                v = tr.transform_points(v, in_place=True)
                n = tr.transform_vectors(n, in_place=True)
                # Rotate so we know which points are closest to C1 and C4 again now we're in plane
                tr = vector_rotation((v1[0]-v4[0],v1[1]-v4[1],v1[2]-v4[2]),(C1[0]-C4[0],C1[1]-C4[1],C1[2]-C4[2]))
                v = tr.transform_points(v, in_place=True)
                n = tr.transform_vectors(n, in_place=True)
                # Translate to the centre of the sugar ring
                vp = Place(origin = (sugar_centre_x, sugar_centre_y, sugar_centre_z))
                v = vp.transform_points(v)
                d.set_geometry(v, n, t)
            elif sugarname in squaresugars:
                # Create the shape
                v, n, t = box_geometry((-1.4,-1.4,-0.5),(1.4,1.4,0.5))
                # Rotate to match linkage positions
                tl = vector_rotation((v[0,0]-v[2,0],v[0,1]-v[2,1],v[0,2]-v[2,2]),(C1[0]-C4[0],C1[1]-C4[1],C1[2]-C4[2]))
                v = tl.transform_points(v, in_place=True)
                n = tl.transform_vectors(n, in_place=True)
                # Rotate to match the plane of the sugar ring
                tr = vector_rotation((n[20,0],n[20,1],n[20,2]),(sugar_plane_i, sugar_plane_j, sugar_plane_k))
                v = tr.transform_points(v, in_place=True)
                n = tr.transform_vectors(n, in_place=True)
                # Rotate to match linkage positions again now we're in plane
                tl = vector_rotation((v[0,0]-v[2,0],v[0,1]-v[2,1],v[0,2]-v[2,2]),(C1[0]-C4[0],C1[1]-C4[1],C1[2]-C4[2]))
                v = tl.transform_points(v, in_place=True)
                n = tl.transform_vectors(n, in_place=True)
                # Translate to the centre of the sugar ring
                vp = Place(origin = (sugar_centre_x, sugar_centre_y, sugar_centre_z))
                v = vp.transform_points(v)
                d.set_geometry(v, n, t)
            elif sugarname in trianglesugars:
                # Create the shape
                v, n, t = triangular_prism_geometry(3.2,1.0) 
                # if sugar["num_bonds"] == 1:
                #     # Rotate to match linkage positions of 2D SNFG based on number of bonds
                #     tl = vector_rotation(((v[2,0]+v[6,0])/2.0-v[0,0],(v[2,1]+v[6,1])/2.0-v[0,1],(v[2,2]+v[6,2])/2.0)-v[0,2],(C1[0]-C4[0],C1[1]-C4[1],C1[2]-C4[2]))
                #     v = tl.transform_points(v, in_place=True)
                #     n = tl.transform_vectors(n, in_place=True)  
                #     # Rotate to match the plane of the sugar ring
                #     tr = vector_rotation((n[15,0],n[15,1],n[15,2]),(sugar_plane_i, sugar_plane_j, sugar_plane_k))
                #     v = tr.transform_points(v, in_place=True)
                #     n = tr.transform_vectors(n, in_place=True)
                #     # Rotate to match linkage positions of 2D SNFG based on number of bonds again now we're in plane
                #     tl = vector_rotation(((v[2,0]+v[6,0])/2.0-v[0,0],(v[2,1]+v[6,1])/2.0-v[0,1],(v[2,2]+v[6,2])/2.0)-v[0,2],(C1[0]-C4[0],C1[1]-C4[1],C1[2]-C4[2]))
                #     v = tl.transform_points(v, in_place=True)
                #     n = tl.transform_vectors(n, in_place=True)  
                #     # Translate to correct location based on vertex of triangle
                #     vp = Place(origin = (C4[0], C4[1], C4[2]))
                #     v = vp.transform_points(v)
                # elif sugar["num_bonds"] == 2: 
                #     # Rotate to match linkage positions of 2D SNFG based on number of bonds 
                #     tl = vector_rotation((v[8,0]-v[10,0],v[8,1]-v[10,1],v[8,2]-v[10,2]),(C1[0]-C3[0],C1[1]-C3[1],C1[2]-C3[2]))
                #     v = tl.transform_points(v, in_place=True)
                #     n = tl.transform_vectors(n, in_place=True)
                #     # Rotate to match the plane of the sugar ring
                #     tr = vector_rotation((n[15,0],n[15,1],n[15,2]),(sugar_plane_i, sugar_plane_j, sugar_plane_k))
                #     v = tr.transform_points(v, in_place=True)
                #     n = tr.transform_vectors(n, in_place=True)  
                #     # Rotate to match linkage positions of 2D SNFG based on number of bonds again now we're in plane
                #     tl = vector_rotation((v[8,0]-v[10,0],v[8,1]-v[10,1],v[8,2]-v[10,2]),(C1[0]-C3[0],C1[1]-C3[1],C1[2]-C3[2]))
                #     v = tl.transform_points(v, in_place=True)
                #     n = tl.transform_vectors(n, in_place=True)
                #     # Translate to correct location based on vertex of triangle
                #     vp = Place(origin = (C2[0], C2[1], C2[2]))
                #     v = vp.transform_points(v)
                # else: 
                #     # Rotate to match linkage positions of 2D SNFG based on number of bonds
                #     tl = vector_rotation((v[0,0]-(v[2,0]+v[6,0])/2.0,v[0,1]-(v[2,1]+v[6,1])/2.0,v[0,2]-(v[2,2]+v[6,2])/2.0),(C1[0]-C4[0],C1[1]-C4[1],C1[2]-C4[2]))
                #     v = tl.transform_points(v, in_place=True)
                #     n = tl.transform_vectors(n, in_place=True)  
                #     # Rotate to match the plane of the sugar ring
                #     tr = vector_rotation((n[15,0],n[15,1],n[15,2]),(sugar_plane_i, sugar_plane_j, sugar_plane_k))
                #     v = tr.transform_points(v, in_place=True)
                #     n = tr.transform_vectors(n, in_place=True)
                #     # Rotate to match linkage positions of 2D SNFG based on number of bonds again now we're in plane
                #     tl = vector_rotation(((v[2,0]+v[6,0])/2.0-v[0,0],(v[2,1]+v[6,1])/2.0-v[0,1],(v[2,2]+v[6,2])/2.0)-v[0,2],(C1[0]-C4[0],C1[1]-C4[1],C1[2]-C4[2]))
                #     v = tl.transform_points(v, in_place=True)
                #     n = tl.transform_vectors(n, in_place=True)  
                #     # Translate to correct location based on vertex of triangle
                #     vp = Place(origin = (C4[0], C4[1], C4[2]))
                #     v = vp.transform_points(v)
                # Rotate to match linkage positions of 2D SNFG based on number of bonds
                tl = vector_rotation((v[0,0]-(v[2,0]+v[6,0])/2.0,v[0,1]-(v[2,1]+v[6,1])/2.0,v[0,2]-(v[2,2]+v[6,2])/2.0),(C1[0]-C4[0],C1[1]-C4[1],C1[2]-C4[2]))
                v = tl.transform_points(v, in_place=True)
                n = tl.transform_vectors(n, in_place=True)  
                # Rotate to match the plane of the sugar ring
                tr = vector_rotation((n[15,0],n[15,1],n[15,2]),(sugar_plane_i, sugar_plane_j, sugar_plane_k))
                v = tr.transform_points(v, in_place=True)
                n = tr.transform_vectors(n, in_place=True)
                # Rotate to match linkage positions of 2D SNFG based on number of bonds again now we're in plane
                tl = vector_rotation((v[0,0]-(v[2,0]+v[6,0])/2.0,v[0,1]-(v[2,1]+v[6,1])/2.0,v[0,2]-(v[2,2]+v[6,2])/2.0),(C1[0]-C4[0],C1[1]-C4[1],C1[2]-C4[2]))
                v = tl.transform_points(v, in_place=True)
                n = tl.transform_vectors(n, in_place=True)  
                #     # Translate to correct location based on vertex of triangle
                vp = Place(origin = (C1[0], C1[1], C1[2]))
                v = vp.transform_points(v)
                d.set_geometry(v, n, t)
            elif sugarname in diamondsugars:
                # Create the shape
                v, n, t = box_geometry((-1.4,-1.4,-0.5),(1.4,1.4,0.5))
                # Rotate to match linkage positions
                tl = vector_rotation((v[9,0]-v[10,0],v[9,1]-v[10,1],v[9,2]-v[10,2]),(C2[0]-C5[0],C2[1]-C5[1],C2[2]-C5[2]))
                v = tl.transform_points(v, in_place=True)
                n = tl.transform_vectors(n, in_place=True)
                # Rotate to match the plane of the sugar ring
                tr = vector_rotation((n[20,0],n[20,1],n[20,2]),(sugar_plane_i, sugar_plane_j, sugar_plane_k))
                v = tr.transform_points(v, in_place=True)
                n = tr.transform_vectors(n, in_place=True)
                # Rotate to match linkage positions again now we're in plane
                tl = vector_rotation((v[9,0]-v[10,0],v[9,1]-v[10,1],v[9,2]-v[10,2]),(C2[0]-C5[0],C2[1]-C5[1],C2[2]-C5[2]))
                v = tl.transform_points(v, in_place=True)
                n = tl.transform_vectors(n, in_place=True)
                # Translate to the centre of the sugar ring
                vp = Place(origin = (sugar_centre_x, sugar_centre_y, sugar_centre_z))
                v = vp.transform_points(v)
                d.set_geometry(v, n, t)
            elif sugarname in starsugars:
                v, n, t = star_geometry(2.0,1)
                # Rotate to match linkage positions
                tl = vector_rotation((v[5,0]-v[0,0],v[5,1]-v[0,1],v[5,2]-v[0,2]),(C1[0]-C4[0],C1[1]-C4[1],C1[2]-C4[2]))
                v = tl.transform_points(v, in_place=True)
                n = tl.transform_vectors(n, in_place=True)
                # Rotate to match the plane of the sugar ring
                tr = vector_rotation((n[0,0],n[0,1],n[0,2]),(sugar_plane_i, sugar_plane_j, sugar_plane_k))
                v = tr.transform_points(v, in_place=True)
                n = tr.transform_vectors(n, in_place=True)
                # Rotate to match linkage positions again now we're in plane
                tl = vector_rotation((v[5,0]-v[0,0],v[5,1]-v[0,1],v[5,2]-v[0,2]),(C1[0]-C4[0],C1[1]-C4[1],C1[2]-C4[2]))
                v = tl.transform_points(v, in_place=True)
                n = tl.transform_vectors(n, in_place=True)
                # Translate to the centre of the sugar ring
                vp = Place(origin = (sugar_centre_x, sugar_centre_y, sugar_centre_z))
                v = vp.transform_points(v)
                d.set_geometry(v, n, t)
            else:
                v, n, t = hexagon_geometry(1.5, 1.0)
                # Rotate to match linkage positions
                tl = vector_rotation(((v[6,0]+v[7,0])/2-(v[18,0]+v[19,0])/2,(v[6,1]+v[7,1])/2-(v[18,1]+v[19,1])/2,(v[6,2]+v[7,2])/2-(v[18,2]+v[19,2])/2),(C1[0]-C4[0],C1[1]-C4[1],C1[2]-C4[2]))
                v = tl.transform_points(v, in_place=True)
                n = tl.transform_vectors(n, in_place=True)
                # Rotate to match the plane of the sugar ring
                tr = vector_rotation((n[0,0],n[0,1],n[0,2]),(sugar_plane_i, sugar_plane_j, sugar_plane_k))
                v = tr.transform_points(v, in_place=True)
                n = tr.transform_vectors(n, in_place=True)
                # Rotate to match linkage positions
                tl = vector_rotation(((v[6,0]+v[7,0])/2-(v[18,0]+v[19,0])/2,(v[6,1]+v[7,1])/2-(v[18,1]+v[19,1])/2,(v[6,2]+v[7,2])/2-(v[18,2]+v[19,2])/2),(C1[0]-C4[0],C1[1]-C4[1],C1[2]-C4[2]))
                v = tl.transform_points(v, in_place=True)
                n = tl.transform_vectors(n, in_place=True)
                # Translate to the centre of the sugar ring
                vp = Place(origin = (sugar_centre_x, sugar_centre_y, sugar_centre_z))
                v = vp.transform_points(v)
                d.set_geometry(v, n, t)
            # Colour according to sugar type
            if sugarname in bluesugars:
                d.color = blue 
            elif sugarname in greensugars:
                d.color = green
            elif sugarname in redsugars:
                d.color = red
            elif sugarname in yellowsugars:
                d.color = yellow
            elif sugarname in orangesugars:
                d.color = orange
            elif sugarname in cyansugars:
                d.color = cyan
            elif sugarname in purplesugars:
                d.color = purple
            else:
                d.color = grey
            dm.add_drawing(d)
            for i, link in enumerate(glycan["Torsions"]):
                if i == 0:
                    pos1 = [link["x1"],link["y1"],link["z1"]]
                if link["chainID"] == sugar["chainID"] and link["sugar_1_resID"] == sugar["resID"] and link["sugar_1"] == sugarname:
                    if sugarname in circlesugars:
                        if link["atom_number_1"] == "1":
                            pos1 = (v[0,:] + v[nc,:])/2.0
                        elif link["atom_number_1"] == "2":
                            pos1 = (v[int(5*nc/6),:]+v[int(5*nc/6)+nc,:])/2.0
                        elif link["atom_number_1"] == "3":
                            pos1 = (v[int(2*nc/3),:]+v[int(2*nc/3)+nc,:])/2.0
                        elif link["atom_number_1"] == "4":
                            pos1 = (v[int(nc/2),:]+v[int(nc/2)+nc,:])/2.0
                        elif link["atom_number_1"] == "5":
                            pos1 = (v[int(nc/3),:]+v[int(nc/3)+nc,:])/2.0
                        elif link["atom_number_1"] == "6":
                            pos1 = (v[int(nc/3),:]+v[int(nc/3)+nc,:])/2.0
                        else: 
                            pos1 = [link["x1"],link["y1"],link["z1"]]
                    elif sugarname in squaresugars:
                        if link["atom_number_1"] == "1":
                            pos1 = (v[4,:] + v[5,:] + v[6,:] + v[7,:])/4.0
                        elif link["atom_number_1"] == "2":
                            if "3" in sugar["link_atoms"]:
                                pos1 = (v[0,:]*2 + v[1,:]*2 + v[2,:] + v[3,:])/6.0
                            else:
                                pos1 = (v[0,:] + v[1,:] + v[2,:] + v[3,:])/4.0
                        elif link["atom_number_1"] == "3":
                            if "2" in sugar["link_atoms"]:
                                pos1 = (v[0,:] + v[1,:] + v[2,:]*2 + v[3,:]*2)/6.0
                            else:
                                pos1 = (v[0,:] + v[1,:] + v[2,:] + v[3,:])/4.0
                        elif link["atom_number_1"] == "4":
                            pos1 = (v[16,:] + v[17,:] + v[18,:] + v[19,:])/4.0
                        elif link["atom_number_1"] == "5":
                            pos1 = (v[12,:] + v[13,:] + v[14,:] + v[15,:])/4.0
                        elif link["atom_number_1"] == "6":
                            pos1 = (v[12,:] + v[13,:] + v[14,:] + v[15,:])/4.0
                        else: 
                            pos1 = [link["x1"],link["y1"],link["z1"]]
                    elif sugarname in trianglesugars:
                        # if sugar["num_bonds"] == 1:
                        #     if link["atom_number_1"] == "1":
                        #         pos1 = (v[8,:] + v[9,:] + v[10,:] + v[11,:])/4.0
                        #     else:
                        #         pos1 = [link["x1"],link["y1"],link["z1"]]
                        # elif sugar["num_bonds"] == 2:
                        #     if link["atom_number_1"] == "1":
                        #         pos1 = (v[4,:]*2 + v[5,:]*2 + v[6,:] + v[7,:])/6.0
                        #     elif link["atom_number_1"] == "3":
                        #         pos1 = (v[0,:]*2 + v[1,:]*2 + v[2,:] + v[3,:])/6.0
                        #     else:
                        #         pos1 = [link["x1"],link["y1"],link["z1"]]
                        # else:
                        #     pos1 = [link["x1"],link["y1"],link["z1"]]
                        if link["atom_number_1"] == "1":
                            pos1 = (v[0,:]+v[1,:])/2.0
                        elif link["atom_number_1"] == "2":
                            pos1 = (v[0,:] + v[1,:] + v[2,:] + v[3,:])/4.0
                        elif link["atom_number_1"] == "3":
                            pos1 = (v[2,:]+v[3,:])/2.0
                        elif link["atom_number_1"] == "4":
                            pos1 = (v[8,:] + v[9,:] + v[10,:] + v[11,:])/4.0
                        elif link["atom_number_1"] == "5":
                            pos1 = (v[8,:]+v[9,:])/2.0
                        elif link["atom_number_1"] == "6":
                            pos1 = (v[8,:]+v[9,:])/2.0
                        #else:
                        #    pos1 = [link["x1"],link["y1"],link["z1"]]
                    elif sugarname in diamondsugars:
                        if link["atom_number_1"] == "1":
                            pos1 = (v[0,:] + v[1,:])/2.0
                        elif link["atom_number_1"] == "2":
                            pos1 = (v[0,:] + v[1,:])/2.0
                        elif link["atom_number_1"] == "3":
                            if "4" in sugar["link_atoms"]:
                                pos1 = (v[0,:] + v[1,:] + v[2,:] + v[3,:])/4.0
                            else:
                                pos1 = (v[2,:] + v[3,:])/2.0
                        elif link["atom_number_1"] == "4":
                            if "3" in sugar["link_atoms"]:
                                pos1 = (v[16,:] + v[17,:] + v[18,:] + v[19,:])/4.0
                            else:
                                pos1 = (v[2,:] + v[3,:])/2.0
                        elif link["atom_number_1"] == "5":
                            pos1 = (v[14,:] + v[15,:])/2.0
                        elif link["atom_number_1"] == "6":
                            pos1 = (v[12,:] + v[13,:])/2.0
                        elif link["atom_number_1"] == "7":
                            pos1 = (v[12,:] + v[13,:])/2.0
                        elif link["atom_number_1"] == "8":
                            pos1 = (v[12,:] + v[13,:])/2.0
                        else: 
                            pos1 = [link["x1"],link["y1"],link["z1"]]
                    elif sugarname in starsugars:
                        if link["atom_number_1"] == "1":
                            pos1 = (v[5,:] + v[55,:])/2.0
                        elif link["atom_number_1"] == "2":
                            pos1 = (v[3,:] + v[53,:])/2.0
                        elif link["atom_number_1"] == "3":
                            pos1 = (v[1,:] + v[51,:])/2.0
                        elif link["atom_number_1"] == "4":
                            pos1 = (v[0,:] + v[50,:])/2.0
                        elif link["atom_number_1"] == "5":
                            pos1 = (v[9,:] + v[59,:])/2.0
                        elif link["atom_number_1"] == "6":
                            pos1 = (v[9,:] + v[59,:])/2.0
                        else:
                            pos1 = [link["x1"],link["y1"],link["z1"]]
                    else:
                        if link["atom_number_1"] == "1":
                            pos1 = (v[0,:] + v[1,:] + v[30,:] + v[31,:])/4.0
                        elif link["atom_number_1"] == "2":
                            pos1 = (v[1,:] + v[2,:] + v[31,:] + v[32,:])/4.0
                        elif link["atom_number_1"] == "3":
                            pos1 = (v[2,:] + v[3,:] + v[32,:] + v[33,:])/4.0
                        elif link["atom_number_1"] == "4":
                            pos1 = (v[3,:] + v[4,:] + v[33,:] + v[34,:])/4.0
                        elif link["atom_number_1"] == "5":
                            pos1 = (v[4,:] + v[5,:] + v[34,:] + v[35,:])/4.0
                        elif link["atom_number_1"] == "6":
                            pos1 = (v[4,:] + v[5,:] + v[34,:] + v[35,:])/4.0
                        else:
                            pos1 = [link["x1"],link["y1"],link["z1"]]
                    glycan["Torsions"][i]["x1"] = pos1[0]
                    glycan["Torsions"][i]["y1"] = pos1[1]
                    glycan["Torsions"][i]["z1"] = pos1[2]
                if link["chainID"] == sugar["chainID"] and link["sugar_2_resID"] == sugar["resID"] and link["sugar_2"] == sugarname:
                    if sugarname in circlesugars:
                        if link["atom_number_2"] == "1":
                            pos2 = (v[0,:] + v[nc,:])/2.0
                        elif link["atom_number_2"] == "2":
                            pos2 = (v[int(5*nc/6),:]+v[int(5*nc/6)+nc,:])/2.0
                        elif link["atom_number_2"] == "3":
                            pos2 = (v[int(2*nc/3),:]+v[int(2*nc/3)+nc,:])/2.0
                        elif link["atom_number_2"] == "4":
                            pos2 = (v[int(nc/2),:]+v[int(nc/2)+nc,:])/2.0
                        elif link["atom_number_2"] == "5":
                            pos2 = (v[int(nc/3),:]+v[int(nc/3)+nc,:])/2.0
                        elif link["atom_number_2"] == "6":
                            pos2 = (v[int(nc/3),:]+v[int(nc/3)+nc,:])/2.0
                        else: 
                            pos2 = [link["x2"],link["y2"],link["z2"]] 
                    elif sugarname in squaresugars:
                        if link["atom_number_2"] == "1":
                            pos2 = (v[4,:] + v[5,:] + v[6,:] + v[7,:])/4.0
                        elif link["atom_number_2"] == "2":
                            if "3" in sugar["link_atoms"]:
                                pos2 = (v[0,:]*2 + v[1,:]*2 + v[2,:] + v[3,:])/6.0
                            else:
                                pos2 = (v[0,:] + v[1,:] + v[2,:] + v[3,:])/4.0
                        elif link["atom_number_2"] == "3":
                            if "2" in sugar["link_atoms"]:
                                pos2 = (v[0,:] + v[1,:] + v[2,:]*2 + v[3,:]*2)/6.0
                            else:
                                pos2 = (v[0,:] + v[1,:] + v[2,:] + v[3,:])/4.0
                        elif link["atom_number_2"] == "4":
                            pos2 = (v[16,:] + v[17,:] + v[18,:] + v[19,:])/4.0
                        elif link["atom_number_2"] == "5":
                            pos2 = (v[12,:] + v[13,:] + v[14,:] + v[15,:])/4.0
                        elif link["atom_number_2"] == "6":
                            pos2 = (v[12,:] + v[13,:] + v[14,:] + v[15,:])/4.0
                        else: 
                            pos2 = [link["x2"],link["y2"],link["z2"]]
                    elif sugarname in trianglesugars:
                        # if sugar["num_bonds"] == 1:
                        #     if link["atom_number_2"] == "1":
                        #         pos2 = (v[8,:] + v[9,:] + v[10,:] + v[11,:])/4.0
                        #     else:
                        #         pos2 = [link["x2"],link["y2"],link["z2"]]
                        # elif sugar["num_bonds"] == 2:
                        #     if link["atom_number_2"] == "1":
                        #         pos2 = (v[4,:]*2 + v[5,:]*2 + v[6,:] + v[7,:])/6.0
                        #     elif link["atom_number_2"] == "3":
                        #         pos2 = (v[0,:]*2 + v[1,:]*2 + v[2,:] + v[3,:])/6.0
                        #     else:
                        #         pos2 = [link["x2"],link["y2"],link["z2"]]
                        # else:
                        #     pos2 = [link["x2"],link["y2"],link["z2"]]
                        if link["atom_number_2"] == "1":
                            pos2 = (v[0,:]+v[1,:])/2.0
                        elif link["atom_number_2"] == "2":
                            pos2 = (v[0,:] + v[1,:] + v[2,:] + v[3,:])/4.0
                        elif link["atom_number_2"] == "3":
                            pos2 = (v[2,:]+v[3,:])/2.0
                        elif link["atom_number_2"] == "4":
                            pos2 = (v[8,:] + v[9,:] + v[10,:] + v[11,:])/4.0
                        elif link["atom_number_2"] == "5":
                            pos2 = (v[8,:]+v[9,:])/2.0
                        elif link["atom_number_2"] == "6":
                            pos2 = (v[8,:]+v[9,:])/2.0
                        #else:
                        #    pos2 = [link["x2"],link["y2"],link["z2"]]
                    elif sugarname in diamondsugars:
                        if link["atom_number_2"] == "1":
                            pos2 = (v[0,:] + v[1,:])/2.0
                        elif link["atom_number_2"] == "2":
                            pos2 = (v[0,:] + v[1,:])/2.0
                        elif link["atom_number_2"] == "3":
                            if "4" in sugar["link_atoms"]:
                                pos2 = (v[0,:] + v[1,:] + v[2,:] + v[3,:])/4.0
                            else:
                                pos2 = (v[2,:] + v[3,:])/2.0
                        elif link["atom_number_2"] == "4":
                            if "3" in sugar["link_atoms"]:
                                pos2 = (v[16,:] + v[17,:] + v[18,:] + v[19,:])/4.0
                            else:
                                pos2 = (v[2,:] + v[3,:])/2.0
                        elif link["atom_number_2"] == "5":
                            pos2 = (v[14,:] + v[15,:])/2.0
                        elif link["atom_number_2"] == "6":
                            pos2 = (v[12,:] + v[13,:])/2.0
                        elif link["atom_number_2"] == "7":
                            pos2 = (v[12,:] + v[13,:])/2.0
                        elif link["atom_number_2"] == "8":
                            pos2 = (v[12,:] + v[13,:])/2.0
                        else: 
                            pos2 = [link["x2"],link["y2"],link["z2"]]
                    elif sugarname in starsugars:
                        if link["atom_number_2"] == "1":
                            pos2 = (v[5,:] + v[55,:])/2.0
                        elif link["atom_number_2"] == "2":
                            pos2 = (v[3,:] + v[53,:])/2.0
                        elif link["atom_number_2"] == "3":
                            pos2 = (v[1,:] + v[51,:])/2.0
                        elif link["atom_number_2"] == "4":
                            pos2 = (v[0,:] + v[50,:])/2.0
                        elif link["atom_number_2"] == "5":
                            pos2 = (v[9,:] + v[59,:])/2.0
                        elif link["atom_number_2"] == "6":
                            pos2 = (v[9,:] + v[59,:])/2.0
                        else:
                            pos2 = [link["x2"],link["y2"],link["z2"]]
                    else:
                        if link["atom_number_2"] == "1":
                            pos2 = (v[0,:] + v[1,:] + v[30,:] + v[31,:])/4.0
                        elif link["atom_number_2"] == "2":
                            pos2 = (v[1,:] + v[2,:] + v[31,:] + v[32,:])/4.0
                        elif link["atom_number_2"] == "3":
                            pos2 = (v[2,:] + v[3,:] + v[32,:] + v[33,:])/4.0
                        elif link["atom_number_2"] == "4":
                            pos2 = (v[3,:] + v[4,:] + v[33,:] + v[34,:])/4.0
                        elif link["atom_number_2"] == "5":
                            pos2 = (v[4,:] + v[5,:] + v[34,:] + v[35,:])/4.0
                        elif link["atom_number_2"] == "6":
                            pos2 = (v[4,:] + v[5,:] + v[34,:] + v[35,:])/4.0
                        else:
                            pos2 = [link["x2"],link["y2"],link["z2"]]
                    glycan["Torsions"][i]["x2"] = pos2[0]
                    glycan["Torsions"][i]["y2"] = pos2[1]
                    glycan["Torsions"][i]["z2"] = pos2[2]
        for link in glycan["Torsions"]:
            pos1 = [link["x1"],link["y1"],link["z1"]]
            pos2 = [link["x2"],link["y2"],link["z2"]]
            if pos1[0] == 0 and pos1[1] == 0 and pos1[2] == 0 and pos2[0] == 0 and pos2[1] == 0 and pos2[2] == 0:
                continue
            # Draw the link between glycans
            dl = Drawing('link')
            linklength = np.sqrt((pos1[0] - pos2[0])**2 + (pos1[1] - pos2[1])**2  + (pos1[2] - pos2[2])**2)
            linkmidpoint = ((pos1[0] + pos2[0])/2 , (pos1[1] + pos2[1])/2  , (pos1[2] + pos2[2])/2)
            linkvector = (pos1[0]-pos2[0],pos1[1]-pos2[1],pos1[2]-pos2[2])
            vlink, nlink, tlink = cylinder_geometry(radius = 0.25, height = linklength, nc=25)
            # Rotate to match link orientation
            trlink = vector_rotation((0,0,1),linkvector)
            vlink = trlink.transform_points(vlink, in_place=True)
            nlink = trlink.transform_vectors(nlink, in_place=True)
            # Translate to the centre of the sugar ring
            vplink = Place(origin = linkmidpoint)
            vlink = vplink.transform_points(vlink)
            dl.set_geometry(vlink, nlink, tlink)
            dm.add_drawing(dl)
            # Add a sphere at either end of bond to smooth the join
            da1 = Drawing('atom1')
            va1, na1, ta1 = sphere_geometry(80)
            vpa1 = Place(origin = (pos1[0],pos1[1],pos1[2]))*scale(0.4)
            va1 = vpa1.transform_points(va1)
            da1.set_geometry(va1, na1, ta1)
            dm.add_drawing(da1)
            da2 = Drawing('atom2')
            va2, na2, ta2 = sphere_geometry(80)
            vpa2 = Place(origin = (pos2[0],pos2[1],pos2[2]))*scale(0.4)
            va2 = vpa2.transform_points(va2)
            da2.set_geometry(va2, na2, ta2)
            dm.add_drawing(da2)
    session.models.add([dm])

def triangular_prism_geometry(l,h):
    # Return vertex, normal vector, and triangle arrays for triangular prism
    from numpy import array, sqrt, float32, int32
    
    #
    #          v4___v5         y  z
    #          /\  / \         | /
    #         /  \/   \        |/______ x
    #        /   /\    \     
    #       /   /  \    \     
    #      /   /    \    \     
    #     /  v1------\----v3  
    #  v0 ------------ v2  
    #
    x = l
    y = sqrt(3.0)*l/2.0
    z = h/2.0
    vertices = array([
        # -x, v0 - v1 - v4 - v5 
        [ 0,    0, -z],
        [ 0,    0,  z],
        [ x/2,  y, -z],
        [ x/2,  y,  z],

        # -y, v0 - v1 - v2 - v3
        [ 0,    0, -z],
        [ 0,    0,  z],
        [ x,    0, -z],
        [ x,    0,  z],

        # x,  v2 - v3 - v4 - v5
        [ x,    0, -z],
        [ x,    0,  z],
        [ x/2,  y, -z],
        [ x/2,  y,  z],

        # -z, v0 - v2 - v4
        [ 0,    0, -z],
        [ x,    0, -z],
        [ x/2,  y, -z],

        # z, v1 - v3 - v5
        [ 0,    0,  z],
        [ x,    0,  z],
        [ x/2,  y,  z],
    ],dtype=float32)

    normals = array([ # FLAG: Need to fix these by trial and error but the basic shape is there now.
        # -x, v0 - v1 - v4 - v5 
        [ -1, 1,  0],
        [ -1, 1,  0],
        [ -1, 1,  0],
        [ -1, 1,  0],

        # -y, v0 - v1 - v2 - v3
        [ 0,  1,  0],
        [ 0,  1,  0],
        [ 0,  1,  0],
        [ 0,  1,  0],

        # x,  v2 - v3 - v4 - v5
        [-1,  1,  0],
        [-1,  1,  0],
        [-1,  1,  0],
        [-1,  1,  0],

        # -z, v0 - v2 - v4
        [ 0,  0,  1],
        [ 0,  0,  1],
        [ 0,  0,  1],

        # z, v1 - v3 - v5
        [ 0,  0,  1],
        [ 0,  0,  1],
        [ 0,  0,  1],
    ],dtype=float32)

    triangles = array([
        # -x, v0 - v1 - v4 - v5 
        [ 0, 1, 2],
        [ 2, 1, 3],
        # -y, v0 - v1 - v2 - v3
        [ 4, 5, 6],
        [ 6, 5, 7],
        # x,  v2 - v3 - v4 - v5
        [ 8, 9,10],
        [10, 9,11],
        # -z, v0 - v2 - v4
        [12,13,14],
        # z, v1 - v3 - v5
        [15,16,17],
    ],dtype=int32)
    return vertices, normals, triangles
def star_geometry(rout, h):
    from numpy import array, zeros, sqrt, cos, sin, pi, float32, int32
    xo = zeros(5)
    yo = zeros(5)
    xi = zeros(5)
    yi = zeros(5)
    # Set inner radius of star
    ratio = 2/(3.0 + sqrt(5))
    rin = rout * ratio
    # Define the vertices of the star by
    for i in range(5):
        # 5 points are equally spaced on outer circle
        xo[i] = rout*cos(2*pi*i/5.0 + pi/2.0)
        yo[i] = rout*sin(2*pi*i/5.0 + pi/2.0)
        # Other 5 vertices are equally spaced on inner circle, offset by 36 deg
        xi[i] = rin*cos(2*pi*i/5.0 + pi/2.0 + pi/5.0)
        yi[i] = rin*sin(2*pi*i/5.0 + pi/2.0 + pi/5.0)
    z = h/2.0
    mx = (xo[0]+xi[4])/2.0
    my = (yo[0]+yi[4])/2
    vertices = array([
        [xo[0],yo[0],z],    # 0
        [xi[0],yi[0],z],    # 1 
        [xo[1],yo[1],z],    # 2
        [xi[1],yi[1],z],    # 3
        [xo[2],yo[2],z],    # 4
        [xi[2],yi[2],z],    # 5
        [xo[3],yo[3],z],    # 6
        [xi[3],yi[3],z],    # 7
        [xo[4],yo[4],z],    # 8
        [xi[4],yi[4],z],    # 9

        [xo[0],yo[0],z],    # 10
        [xi[0],yi[0],z],    # 11
        [xo[0],yo[0],-z],   # 12
        [xi[0],yi[0],-z],   # 13

        [xi[0],yi[0],z],    # 14
        [xo[1],yo[1],z],    # 15
        [xi[0],yi[0],-z],   # 16
        [xo[1],yo[1],-z],   # 17

        [xo[1],yo[1],z],    # 18
        [xi[1],yi[1],z],    # 19
        [xo[1],yo[1],-z],   # 20
        [xi[1],yi[1],-z],   # 21

        [xi[1],yi[1],z],    # 22
        [xo[2],yo[2],z],    # 23
        [xi[1],yi[1],-z],   # 24
        [xo[2],yo[2],-z],   # 25

        [xo[2],yo[2],z],    # 26
        [xi[2],yi[2],z],    # 27
        [xo[2],yo[2],-z],   # 28
        [xi[2],yi[2],-z],   # 29

        [xi[2],yi[2],z],    # 30
        [xo[3],yo[3],z],    # 31
        [xi[2],yi[2],-z],   # 32
        [xo[3],yo[3],-z],   # 33

        [xo[3],yo[3],z],    # 34
        [xi[3],yi[3],z],    # 35
        [xo[3],yo[3],-z],   # 36
        [xi[3],yi[3],-z],   # 37

        [xi[3],yi[3],z],    # 38
        [xo[4],yo[4],z],    # 39
        [xi[3],yi[3],-z],   # 40
        [xo[4],yo[4],-z],   # 41

        [xo[4],yo[4],z],    # 42
        [xi[4],yi[4],z],    # 43
        [xo[4],yo[4],-z],   # 44
        [xi[4],yi[4],-z],   # 45

        [xi[4],yi[4],z],    # 46
        [xo[0],yo[0],z],    # 47
        #[mx,my,0],
        [xi[4],yi[4],-z],   # 48
        [xo[0],yo[0],-z],   # 49

        [xo[0],yo[0],-z],   # 50
        [xi[0],yi[0],-z],   # 51
        [xo[1],yo[1],-z],   # 52
        [xi[1],yi[1],-z],   # 53
        [xo[2],yo[2],-z],   # 54
        [xi[2],yi[2],-z],   # 55
        [xo[3],yo[3],-z],   # 56
        [xi[3],yi[3],-z],   # 57
        [xo[4],yo[4],-z],   # 58
        [xi[4],yi[4],-z],   # 59
    ],dtype=float32)

    normals = array([
        [0,0,z],            # 0
        [0,0,z],            # 1
        [0,0,z],            # 2
        [0,0,z],            # 3
        [0,0,z],            # 4
        [0,0,z],            # 5
        [0,0,z],            # 6
        [0,0,z],            # 7
        [0,0,z],            # 8
        [0,0,z],            # 9

        [-xi[3],-yi[3],0],  # 10
        [-xi[3],-yi[3],0],  # 11
        [-xi[3],-yi[3],0],  # 12
        [-xi[3],-yi[3],0],  # 13

        [-xi[2],-yi[2],0],  # 14
        [-xi[2],-yi[2],0],  # 15
        [-xi[2],-yi[2],0],  # 16
        [-xi[2],-yi[2],0],  # 17

        [-xi[4],-yi[4],0],  # 18
        [-xi[4],-yi[4],0],  # 19
        [-xi[4],-yi[4],0],  # 20
        [-xi[4],-yi[4],0],  # 21

        [-xi[3],-yi[3],0],  # 22
        [-xi[3],-yi[3],0],  # 23
        [-xi[3],-yi[3],0],  # 24
        [-xi[3],-yi[3],0],  # 25

        [-xi[0],-yi[0],0],  # 26
        [-xi[0],-yi[0],0],  # 27
        [-xi[0],-yi[0],0],  # 28
        [-xi[0],-yi[0],0],  # 29

        [-xi[4],-yi[4],0],  # 30
        [-xi[4],-yi[4],0],  # 31
        [-xi[4],-yi[4],0],  # 32
        [-xi[4],-yi[4],0],  # 33

        [-xi[1],-yi[1],0],  # 34
        [-xi[1],-yi[1],0],  # 35
        [-xi[1],-yi[1],0],  # 36
        [-xi[1],-yi[1],0],  # 37

        [-xi[0],-yi[0],0],  # 38
        [-xi[0],-yi[0],0],  # 39
        [-xi[0],-yi[0],0],  # 40
        [-xi[0],-yi[0],0],  # 41

        [-xi[2],-yi[2],0],  # 42
        [-xi[2],-yi[2],0],  # 43
        [-xi[2],-yi[2],0],  # 44
        [-xi[2],-yi[2],0],  # 45

        [-xi[1],-yi[1],0],  # 46
        [-xi[1],-yi[1],0],  # 47
        [-xi[1],-yi[1],0],  # 48
        [-xi[1],-yi[1],0],  # 49

        [0,0,z],            # 50
        [0,0,z],            # 51
        [0,0,z],            # 52
        [0,0,z],            # 53
        [0,0,z],            # 54
        [0,0,z],            # 55
        [0,0,z],            # 56
        [0,0,z],            # 57
        [0,0,z],            # 58
        [0,0,z],            # 59
    ],dtype=float32)

    triangles = array([
        # Top face
        [0,4,7],
        [2,5,8],
        [0,3,6],

        # Sides
        [12,11,10],
        [11,12,13],

        [16,15,14],
        [15,16,17],

        [20,19,18],
        [19,20,21],

        [24,23,22],
        [23,24,25],

        [28,27,26],
        [27,28,29],

        [32,31,30],
        [31,32,33],

        [36,35,34],
        [35,36,37],

        [40,39,38],
        [39,40,41],

        [44,43,42],
        [43,44,45],

        [48,47,46],
        [47,48,49],

        # Bottom face
        [50,54,57],
        [52,55,58],
        [50,53,56],
    ],dtype=int32)
    return vertices, normals, triangles

def hexagon_geometry(l,h):
    from numpy import array, zeros, sqrt, cos, sin, pi, float32, int32
    vertices = array([
        [-l,     0,             h/2],   # 0
        [-l/2,   sqrt(3)*l/2,   h/2],   # 1
        [ l/2,   sqrt(3)*l/2,   h/2],   # 2
        [ l,     0,             h/2],   # 3
        [ l/2,  -sqrt(3)*l/2,   h/2],   # 4
        [-l/2,  -sqrt(3)*l/2,   h/2],   # 5

        [-l,    0,              h/2],   # 6
        [-l/2,  sqrt(3)*l/2,    h/2],   # 7
        [-l,    0,             -h/2],   # 8
        [-l/2,  sqrt(3)*l/2,   -h/2],   # 9

        [-l/2,  sqrt(3)*l/2,     h/2],  # 10
        [ l/2,  sqrt(3)*l/2,     h/2],  # 11
        [-l/2,  sqrt(3)*l/2,    -h/2],  # 12
        [ l/2,  sqrt(3)*l/2,    -h/2],  # 13

        [ l/2,   sqrt(3)*l/2,     h/2], # 14
        [ l,     0,               h/2], # 15
        [ l/2,   sqrt(3)*l/2,    -h/2], # 16
        [ l,     0,              -h/2], # 17

        [ l,     0,              h/2],  # 18
        [ l/2,  -sqrt(3)*l/2,    h/2],  # 19
        [ l,     0,             -h/2],  # 20
        [ l/2,  -sqrt(3)*l/2,   -h/2],  # 21

        [ l/2,  -sqrt(3)*l/2,    h/2],  # 22
        [-l/2,  -sqrt(3)*l/2,    h/2],  # 23
        [ l/2,  -sqrt(3)*l/2,   -h/2],  # 24
        [-l/2,  -sqrt(3)*l/2,   -h/2],  # 25

        [-l/2,  -sqrt(3)*l/2,    h/2],  # 26
        [-l,     0,              h/2],  # 27
        [-l/2,  -sqrt(3)*l/2,   -h/2],  # 28
        [-l,     0,             -h/2],  # 29

        [-l,    0,              -h/2],  # 30
        [-l/2,  sqrt(3)*l/2,    -h/2],  # 31
        [l/2,   sqrt(3)*l/2,    -h/2],  # 32
        [l,     0,              -h/2],  # 33
        [l/2,   -sqrt(3)*l/2,   -h/2],  # 34
        [-l/2,  -sqrt(3)*l/2,   -h/2],  # 35
    ],dtype=float32)
    normals = array([
        [0,0,1],    # 0
        [0,0,1],    # 1
        [0,0,1],    # 2
        [0,0,1],    # 3
        [0,0,1],    # 4
        [0,0,1],    # 5

        [0,1,0],    # 6
        [0,1,0],    # 7
        [0,1,0],    # 8
        [0,1,0],    # 9

        [1,1,0],    # 10
        [1,1,0],    # 11
        [1,1,0],    # 12
        [1,1,0],    # 13

        [1,-1,0],   # 14
        [1,-1,0],   # 15
        [1,-1,0],   # 16
        [1,-1,0],   # 17

        [0,-1,0],   # 18
        [0,-1,0],   # 19
        [0,-1,0],   # 20
        [0,-1,0],   # 21

        [-1,-1,0],  # 22
        [-1,-1,0],  # 23
        [-1,-1,0],  # 24
        [-1,-1,0],  # 25

        [-1,1,0],   # 26
        [-1,1,0],   # 27
        [-1,1,0],   # 28
        [-1,1,0],   # 29

        [0,0,-1],   # 30
        [0,0,-1],   # 31
        [0,0,-1],   # 32
        [0,0,-1],   # 33
        [0,0,-1],   # 34
        [0,0,-1],   # 35
        ],dtype=float32)
    triangles = array([
        # Top face
        [0,3,1],
        [1,4,2],
        [2,5,3],
        [3,0,4],
        [4,1,5],
        [5,2,0],

        # Sides
        [6,7,8],
        [9,8,7],

        [10,11,12],
        [13,12,11],

        [14,15,16],
        [17,16,15],

        [18,19,20],
        [21,20,19],

        [22,23,24],
        [25,24,23],

        [26,27,28],
        [29,28,27],

        # Bottom face
        [31,33,30],
        [32,34,31],
        [33,35,32],
        [34,30,33],
        [35,31,34],
        [30,32,35],
    ],dtype=int32)
    return vertices, normals, triangles

   

