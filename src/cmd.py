from chimerax.core.commands import CmdDesc, SaveFolderNameArg, IntArg, BoolArg, FloatArg, StringArg, ModelIdArg
from chimerax import atomic
from .main import *



def privateer_validation(session,OutputFolderPath,modelID=None,display=False):
    models = atomic.all_structures(session)
    if modelID != None:
        model = None
        for m in models:
            if m.id_string == modelID:
                model = m
        if model == None:
            session.logger.info(f"Error in running Privateer... Chosen modelID does not correspond to model loaded in the session.")
    if display:
        ValidationReportAllGlycans = privateer_validation_wrapper(session,OutputFolderPath,model,modelID,display)
        return ValidationReportAllGlycans
    else:
        if modelID == None:
            for i,m in enumerate(models):
                privateer_validation_wrapper(session,OutputFolderPath,m,i+1)
            session.logger.info(f"Privateer has run carbohydrate validation on the structure models currently loaded in the session. The validation report is saved in {OutputFolderPath} with filename model-i_privateer-report.csv where i is the model index in the current session.")
        else:
            privateer_validation_wrapper(session,OutputFolderPath,model,modelID)
            session.logger.info(f"Privateer has run carbohydrate validation on the structure model {modelID}. The validation report is saved in {OutputFolderPath} with filename model-{modelID}_privateer-report.csv.")

privateer_validation_desc = CmdDesc(
    required=[
        ("OutputFolderPath", SaveFolderNameArg),
    ],
    optional=[
        ("modelID", StringArg),
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
    from Qt.QtCore import Qt 
    from Qt.QtWidgets import QVBoxLayout
    from chimerax.privateer.tool import PrivateerTool
    tool = session.tools.find_by_class(PrivateerTool)[0]
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
    torsion_tool_window.manage(placement=None, allowed_areas=Qt.LeftDockWidgetArea|Qt.RightDockWidgetArea)

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

# def privateer_glycoblocks(session, modelID):
#     models = atomic.all_structures(session)
#     model = None
#     for m in models:
#         if m.id_string == modelID:
#             model = m
#     if model == None:
#         session.logger.info(f"Error in running Privateer... Chosen modelID does not correspond to model loaded in the session.")
#     Glycans = privateer_validation_wrapper(session,None,model,modelID,True)
#     draw_glycoblocks(session,Glycans,modelID)

def privateer_glycoblocks(session, modelID, auto_update, scroll_resize):
    models = atomic.all_structures(session)
    model = None
    for m in models:
        if m.id_string == modelID:
            model = m
        if model == None:
            session.logger.info(f"Error in running Privateer... Chosen modelID does not correspond to model loaded in the session.")
    from .tool import Glycoblocks
    return Glycoblocks(model,auto_update,scroll_resize)
        
    
privateer_glycoblocks_desc = CmdDesc(
    required=[
        ("modelID", StringArg),
        ("auto_update", BoolArg),
        ("scroll_resize", BoolArg)
    ],
    synopsis="Display 3D symbols for glycans on model"
)