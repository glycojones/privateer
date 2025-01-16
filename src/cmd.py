from chimerax.core.commands import register, Command, CmdDesc, OpenFileNameArg, SaveFolderNameArg
from chimerax import atomic
from .main import *


def privateer_validation(session,OutputFolderPath):
    logger = session.logger
    models = atomic.all_structures(session)
    for i,m in enumerate(models):
        privateer_validation_wrapper(OutputFolderPath,m,i)
    session.logger.info(f"Privateer has run carbohydrate validation on the structure. The validation report is saved in {OutputFolderPath}.")

privateer_validation_desc = CmdDesc(
    required=[
        ('OutputFolderPath', SaveFolderNameArg),
    ],
    synopsis='Provide a validation report on the input structure'
)

