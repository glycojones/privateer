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
    torsion_plot(sugar_1,atom_number_1,sugar_2,atom_number_2,phi,psi)

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