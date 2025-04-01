import pandas as pd
import re
import os
import time
import tempfile

from . import privateer_core as pvt


def privateer_validation_wrapper(OutputFolderPath,m,i,display=False):
    timestr = time.strftime("%Y%m%d-%H%M%S")
    tempdirpath = tempfile.gettempdir()
    InputStructureFilePath = os.path.join(tempdirpath,f"temp_{timestr}.pdb")
    _write_pdb(m,InputStructureFilePath)
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
def _write_pdb(m, filename):
    # Function for writing a pdb from chimeraX model
    # Taken from an example chimeraX bundle
    # Writing then deleting temp.pdb is kinda hacky -- look into proper way to do this converting model too a clipper::mmol object
    atoms = m.atoms
    coords = atoms.coords
    atom_names = atoms.names
    element_names = atoms.element_names
    residues = atoms.residues
    chain_ids = residues.chain_ids
    residue_names = residues.names
    residue_numbers = residues.numbers
    with open(filename, "w") as f:
        serial = 0
        for i in range(len(coords)):
            serial += 1
            atom_name = atom_names[i]
            alt_loc = ' '
            res_name = residue_names[i]
            chain_id = chain_ids[i]
            res_seq = residue_numbers[i]
            res_insert_code = ' '
            x, y, z = coords[i]
            occupancy = 1.0
            b_factor = 0.0
            segment = ' '
            element = element_names[i]
            if len(element) == 1:
                atom_name = ' ' + atom_name
            charge = ' '
            print(_ATOM_FMT % (serial, atom_name, alt_loc,
                               res_name, chain_id, res_seq, res_insert_code,
                               x, y, z, occupancy, b_factor,
                               segment, element, charge),
                               file=f)

def draw_glycoblocks(session,Glycans):
    from chimerax.surface.shapes import cylinder_geometry, box_geometry, sphere_geometry
    from chimerax.core.models import Drawing, Model
    from chimerax.geometry import Place, vector_rotation, scale
    import numpy as np
    dm = Model('privateer glycoblocks',session)
    blue = [0,144,188,255]
    green = [0,166,81,255]
    red = [237,28,36,255]
    orange = [255,127,0,255]
    yellow = [255,221,0,255]
    grey = [128,128,128,255]
    # Glc       = GLC,BGC       = blue      circle
    # Gal       = GAL,GLA       = yellow    circle
    # Man       = MAN,BMA       = green     circle
    # Fuc       = FUC,FCB,FUL   = red       triangle
    # Xyl       = XYS,XYP       = orange    star
    # GlcNAc    = NAG,NDG       = blue      square
    # GalNAc    = NGA,A2G       = yellow    square
    # ManNAc    = BM3,BM7       = green     square
    # GlcN      = GCS,PA1
    # GalN      = 
    # ManN      =
    # GlcA      = BDP,GCU
    # GalA      = GTR,ADA
    # ManA      = MAV,BEM
    # Neu5Gc    =
    # Neu5Ac    = SIA,SLB
    # IdoA      = IDR
    # KDN       = KDM,KDN
    bluesugars = ["GLC","BGC","NAG","NDG"]
    greensugars = ["MAN","BMA","BM3","BM7"]
    redsugars = ["FUC","FCB","FUL"]
    yellowsugars = ["GAL","GLA","NGA","AG2"]
    orangesugars = ["XYS","XYP"]

    circlesugars = ["GLC","BGC","GAL","GLA","MAN","BMA"]
    squaresugars = ["NAG","NDG","NGA","A2G","BM3","BM7"]
    trianglesugars = ["FUC","FCB","FUL"]
    starsugars = []
    diamondsugars = []
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
            d = Drawing('sugar')
            if sugarname in circlesugars:
                # Create the shape
                v, n, t = cylinder_geometry(radius = 1.4, height = 1.0, nc=25)
                # Rotate to match the plane of the sugar ring
                tr = vector_rotation((0,0,1),(sugar_plane_i, sugar_plane_j, sugar_plane_k))
                v = tr.transform_points(v, in_place=True)
                n = tr.transform_vectors(n, in_place=True)
                # Translate to the centre of the sugar ring
                vp = Place(origin = (sugar_centre_x, sugar_centre_y, sugar_centre_z))
                v = vp.transform_points(v)
                d.set_geometry(v, n, t)
            elif sugarname in squaresugars:
                # Create the shape
                v, n, t = box_geometry((-1.4,-1.4,-0.5),(1.4,1.4,0.5))
                # Rotate to match the plane of the sugar ring
                tr = vector_rotation((0,0,1),(sugar_plane_i, sugar_plane_j, sugar_plane_k))
                v = tr.transform_points(v, in_place=True)
                n = tr.transform_vectors(n, in_place=True)
                # Rotate to match linkage positions
                tl = vector_rotation((v[0,0]-v[2,0],v[0,1]-v[2,1],v[0,2]-v[2,2]),(C1[0]-C4[0],C1[1]-C4[1],C1[2]-C4[2]))
                v = tl.transform_points(v, in_place=True)
                n = tl.transform_vectors(n, in_place=True)
                # Translate to the centre of the sugar ring
                vp = Place(origin = (sugar_centre_x, sugar_centre_y, sugar_centre_z))
                v = vp.transform_points(v)
                d.set_geometry(v, n, t)
            elif sugarname in trianglesugars:
                # Create the shape
                v, n, t = triangular_prism_geometry(3.5,1.0) 
                # Rotate to match linkage positions of 2D SNFG based on number of bonds
                if sugar["num_bonds"] == 1:
                    tl = vector_rotation(((v[13,0]+v[14,0])/2.0-v[12,0],(v[13,1]+v[14,1])/2.0-v[12,1],(v[13,2]+v[14,2])/2.0)-v[12,2],(C1[0]-C4[0],C1[1]-C4[1],C1[2]-C4[2]))
                    v = tl.transform_points(v, in_place=True)
                    n = tl.transform_vectors(n, in_place=True)
                    # Rotate to match the plane of the sugar ring
                    tr = vector_rotation((n[15,0],n[15,1],n[15,2]),(sugar_plane_i, sugar_plane_j, sugar_plane_k))
                    v = tr.transform_points(v, in_place=True)
                    n = tr.transform_vectors(n, in_place=True)  
                    # Translate to correct location based on vertex of triangle
                    vp = Place(origin = (C4[0], C4[1], C4[2]))
                    v = vp.transform_points(v)
                if sugar["num_bonds"] == 2:  
                    tl = vector_rotation((v[8,0]-v[10,0],v[8,1]-v[10,1],v[8,2]-v[10,2]),(C1[0]-C3[0],C1[1]-C3[1],C1[2]-C3[2]))
                    v = tl.transform_points(v, in_place=True)
                    n = tl.transform_vectors(n, in_place=True)
                    # Rotate to match the plane of the sugar ring
                    tr = vector_rotation((n[15,0],n[15,1],n[15,2]),(sugar_plane_i, sugar_plane_j, sugar_plane_k))
                    v = tr.transform_points(v, in_place=True)
                    n = tr.transform_vectors(n, in_place=True)  
                    # Translate to correct location based on vertex of triangle
                    vp = Place(origin = (C2[0], C2[1], C2[2]))
                    v = vp.transform_points(v)
                # Translate to the centre of the sugar ring
                #vp = Place(origin = (sugar_centre_x, sugar_centre_y, sugar_centre_z))
                #v = vp.transform_points(v)
                d.set_geometry(v, n, t)
            else:
                v, n, t = cylinder_geometry(radius = 1.5, height = 0.9, nc=25)
                # Rotate to match the plane of the sugar ring
                tr = vector_rotation((0,0,1),(sugar_plane_i, sugar_plane_j, sugar_plane_k))
                v = tr.transform_points(v, in_place=True)
                n = tr.transform_vectors(n, in_place=True)
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
            else:
                d.color = grey
            dm.add_drawing(d)
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
    # x = l/2
    # y = sqrt(3.0)*l/4.0
    # z = h/2.0
    # vertices = array([
    #     # -x, v0 - v1 - v4 - v5 
    #     [-x, -y, -z],
    #     [-x, -y,  z],
    #     [ 0,  y, -z],
    #     [ 0,  y,  z],

    #     # -y, v0 - v1 - v2 - v3
    #     [-x, -y, -z],
    #     [-x, -y,  z],
    #     [ x, -y, -z],
    #     [ x, -y,  z],

    #     # x,  v2 - v3 - v4 - v5
    #     [ x, -y, -z],
    #     [ x, -y,  z],
    #     [ 0,  y, -z],
    #     [ 0,  y,  z],

    #     # -z, v0 - v2 - v4
    #     [-x, -y, -z],
    #     [ x, -y, -z],
    #     [ 0,  y, -z],

    #     # z, v1 - v3 - v5
    #     [-x, -y,  z],
    #     [ x, -y,  z],
    #     [ 0,  y,  z],
    # ],dtype=float32)
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

   

