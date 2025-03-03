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
    linkageString = f"{sugar_1}-{atom_number_2},{atom_number_1}-{sugar_2}" 
    
    if len(phis) < 10:
        session.logger.info(f"Not enough data for linkage {linkageString} to calculate statistics or produce torsion plot.")
        return

    if sugar_1 == "ASN" and sugar_2 == "NAG":
        phimin = -180
        phimax = 180
        psimin = 0
        psimax = 360
    else:
        phimin = -180
        phimax = 180
        psimin = -180
        psimax = 180
    from matplotlib.figure import Figure
    from matplotlib.backends.backend_qtagg import (FigureCanvasQTAgg as Canvas,)
    from Qt.QtWidgets import QVBoxLayout
    from chimerax.privateer.tool import FancierPrivateerTool
    tool = session.tools.find_by_class(FancierPrivateerTool)[0]
    torsion_tool_window = tool.tool_window.create_child_window("Torsion Plot")
    layout = QVBoxLayout()
    fig = Figure()
    axs = fig.add_subplot(111)
    axs.plot(phi,psi,"rx")
    axs.hist2d(phis,psis,bins = 180,range=[[phimin,phimax],[psimin,psimax]], cmin = 1)
    axs.set_title(linkageString)
    axs.set_ylabel("$\psi$")
    axs.set_xlabel("$\phi$")
    fig.tight_layout()
    canvas = Canvas(fig)
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