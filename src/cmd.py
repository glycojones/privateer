from chimerax.core.commands import register, Command, CmdDesc, OpenFileNameArg, SaveFolderNameArg
from .main import *


def privateer_validation(session,InputStructureFilePath,OutputFilePath):
    logger = session.logger
    privateer_validation_wrapper(InputStructureFilePath,OutputFilePath)
    session.logger.info(f"Privateer has run carbohydrate validation on the structure saved at {InputStructureFilePath}. The validation report is saved in {OutputFilePath}.")

privateer_validation_desc = CmdDesc(
    required=[
        ('InputStructureFilePath', OpenFileNameArg),
        ('OutputFolderPath', SaveFolderNameArg),
    ],
    synopsis='Provide a validation report on the input structure'
)

