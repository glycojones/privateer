from chimerax.core.commands import register, Command, CmdDesc, OpenFileNameArg, SaveFolderNameArg
from chimerax import atomic
from .main import *


def privateer_validation(session,OutputFolderPath):
    # FLAG: Change this command to run on specified models. I think I can do this using the model ID.
    logger = session.logger
    models = atomic.all_structures(session)
    for i,m in enumerate(models):
        privateer_validation_wrapper(OutputFolderPath,m,i)
    session.logger.info(f"Privateer has run carbohydrate validation on the structure models currently loaded in the session. The validation report is saved in {OutputFolderPath} with filename model-i_privateer-report.csv where i is the model index in the current session.")

privateer_validation_desc = CmdDesc(
    required=[
        ("OutputFolderPath", SaveFolderNameArg),
    ],
    synopsis="Provide a validation report on the structure models currently loaded in the session. The validation report is saved in the specified location with filename privateer-report-i.csv where i is the model index in the current session."
)

