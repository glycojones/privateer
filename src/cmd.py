from chimerax.core.commands import CmdDesc, SaveFolderNameArg, IntArg, BoolArg, FloatArg, StringArg
from chimerax import atomic
from .main import *



def privateer_validation(session,OutputFolderPath,modelID=None,display=False):
    models = atomic.all_structures(session)
    if display:
        ValidationReportAllGlycans = privateer_validation_wrapper(OutputFolderPath,models[modelID-1],modelID,display)
        return ValidationReportAllGlycans
    else:
        if modelID == None:
            for i,m in enumerate(models):
                privateer_validation_wrapper(OutputFolderPath,m,i+1)
            session.logger.info(f"Privateer has run carbohydrate validation on the structure models currently loaded in the session. The validation report is saved in {OutputFolderPath} with filename model-i_privateer-report.csv where i is the model index in the current session.")
        else:
            if modelID-1 > len(models):
                session.logger.info(f"Error in running Privateer... Chosen modelID exceeds the number of models loaded in the session.")
            privateer_validation_wrapper(OutputFolderPath,models[modelID-1],modelID)
            session.logger.info(f"Privateer has run carbohydrate validation on the structure model {modelID}. The validation report is saved in {OutputFolderPath} with filename model-{modelID}_privateer-report.csv.")

privateer_validation_desc = CmdDesc(
    required=[
        ("OutputFolderPath", SaveFolderNameArg),
    ],
    optional=[
        ("modelID", IntArg),
        ("display", BoolArg)
        ],
    synopsis="Provide a validation report on the model specified, or if no model is specified, on all the structure models currently loaded in the session. The validation report is saved in the specified location with filename model-i_privateer-report.csv where i is the model index in the current session."
)

def privateer_torsion_plot(session,sugar_1,atom_number_1,sugar_2,atom_number_2,phi,psi):
    import os
    dpath = os.path.dirname(os.path.abspath(__file__))
    torsiondatabasefilepath = os.path.join(dpath,"data","linkage_torsions","privateer_torsion_database.json")

    import json
    with open(torsiondatabasefilepath) as json_file:
        torsions = json.load(json_file)
    torsions = torsions["data"]
    phis = []
    psis = []
    for firsts in torsions:
        if firsts["first"] == sugar_1:
            for seconds in firsts["second"]:
                if seconds["sugar"] == sugar_2 and str(seconds["donor_position"]) == str(atom_number_2) and str(seconds["acceptor_position"]) == str(atom_number_1): #Do I need to swap these two positions?
                    for torsion_pairs in seconds["torsions"]:
                        phis.append(torsion_pairs["Phi"])
                        psis.append(torsion_pairs["Psi"])
    if sugar_1 == "ASN" and sugar_2 == "NAG":
        linkageString = f"{sugar_1}-{atom_number_2},{atom_number_1}-{sugar_2}" 
        phimin = -180
        phimax = 180
        psimin = 0
        psimax = 360
    else:
        linkageString = f"{sugar_2}-{atom_number_2},{atom_number_1}-{sugar_1}" 
        phimin = -180
        phimax = 180
        psimin = -180
        psimax = 180
    if len(phis) < 10:
        session.logger.info(f"Not enough data for linkage {linkageString} to calculate statistics or produce torsion plot.")
        return
    from matplotlib.figure import Figure
    from matplotlib.backends.backend_qtagg import (FigureCanvasQTAgg as Canvas,)
    from Qt.QtWidgets import QVBoxLayout
    from chimerax.privateer.tool import FancierPrivateerTool
    tool = session.tools.find_by_class(FancierPrivateerTool)[0]
    torsion_tool_window = tool.tool_window.create_child_window("Torsion Plot")
    ui_area = torsion_tool_window.ui_area
    ui_area.setMinimumHeight(1)
    layout = QVBoxLayout()
    fig = Figure()
    axs = fig.add_subplot(111)
    axs.hist2d(phis,psis,bins = 90,range=[[phimin,phimax],[psimin,psimax]], cmin = 1)
    axs.plot(phi,psi,"rx",label=f"$\phi$ = {round(phi,2)}$^\circ$, $\psi$ = {round(psi,2)}$^\circ$")
    axs.legend()
    axs.set_title(linkageString)
    axs.set_ylabel("$\psi^\circ$")
    axs.set_xlabel("$\phi^\circ$")
    #annotation = axs.annotate(f"({round(phi,2)}$^\circ$,{round(psi,2)}$^\circ$)",xy=(phi,psi),xytext=(10,10),textcoords="offset points",bbox=dict(boxstyle="round", fc="r",alpha=0.4))
    anno = axs.annotate("",xy=(0,0),xytext=(2,2),textcoords="offset points",bbox=dict(boxstyle="round", fc="w",alpha=0.5))
    anno.set_visible(False)

    def on_move(event):
        if event.inaxes:
            anno.xy = (event.xdata,event.ydata)
            anno.set_text(f"({round(event.xdata,2)}$^\circ$,{round(event.ydata,2)}$^\circ$)")
            anno.set_visible(True)
            fig.canvas.draw_idle()
        else:
            anno.set_visible(False)
            fig.canvas.draw_idle()
    fig.tight_layout()
    fig.canvas.mpl_connect("motion_notify_event",on_move)
    canvas = Canvas(fig)
    canvas.setParent(ui_area)
    layout.addWidget(canvas)
    canvas.draw()
    torsion_tool_window.ui_area.setLayout(layout)
    torsion_tool_window.manage('side')

privateer_torsion_plot_desc = CmdDesc(
    required=[
        ("sugar_1", StringArg),
        ("atom_number_1", IntArg),
        ("sugar_2", StringArg),
        ("atom_number_2", IntArg),
        ("phi",FloatArg),
        ("psi",FloatArg),
    ],
    synopsis="Producs a plot of known torsion angles for the specified linkage."
)

def privateer_glycoblocks(session, modelID):
    models = atomic.all_structures(session)
    model = models[modelID-1]
    Glycans = privateer_validation_wrapper(None,model,modelID,True)
    from chimerax.surface.shapes import cylinder_geometry
    from chimerax.core.models import Drawing, Model
    from chimerax.geometry import Place, vector_rotation
    dm = Model('glycoblock test',session)
    blue = [0,0,255,255]
    green = [0,255,0,255]
    red = [255,0,0,255]
    orange = [255,165,0,255]
    yellow = [255,255,0,255]
    grey = [128,128,128,255]
    # Glc = blue circle
    # Gal = yellow circle
    # Man = green circle
    # Fuc = red triangle
    # Xyl = orange star
    # GlcNAc = blue square
    # GalNAc = yellow square
    # ManNAc = green square
    # GlcN
    # GalN
    # ManN
    # GlcA
    # GalA
    # ManA
    # Neu5Gc
    # Neu5Ac
    # IdoA
    # KDN
    bluesugars = ["Glc","GlcNAc"]
    greensugars = ["Man","ManNAc"]
    redsugars = ["Fuc"]
    yellowsugars = ["GalNAc"]
    orangesugars = ["Xyls"]
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
            d = Drawing('privateer glycoblocks')
            # Create the shape
            v, n, t = cylinder_geometry(radius = 1.5, height = 0.25, nc=25)
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
    session.models.add([dm])

privateer_glycoblocks_desc = CmdDesc(
    required=[
        ("modelID", IntArg)
    ],
    synopsis="Display glycans in glycoblock represenstation in the ChimeraX view."
)