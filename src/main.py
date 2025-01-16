import pandas as pd
import re
import os
from . import privateer_core as pvt


def privateer_validation_wrapper(InputStructureFilePath,OutputFolderPath):
    structurefilename = os.path.basename(InputStructureFilePath)
    dpath = os.path.dirname(os.path.abspath(__file__))
    zscorefilepath = os.path.join(dpath,"data","linkage_torsions","privateer_torsions_z_score_database.json")
    AllSugars = pvt.validate(InputStructureFilePath,zscorefilepath)
    df = pd.DataFrame.from_dict(AllSugars)
    csv_out = os.path.join(OutputFolderPath,f"{structurefilename}-privateer-report.csv")
    df.to_csv(csv_out)



