import pandas as pd
import re
import os
from . import privateer_core as pvt


def privateer_validation_wrapper(InputStructureFilePath,OutputFilePath):
    structurefilename = os.path.basename(InputStructureFilePath)
    AllSugars = pvt.validate(InputStructureFilePath)
    df = pd.DataFrame.from_dict(AllSugars)
    csv_out = os.path.join(OutputFilePath,f"{structurefilename}-privateer-report.csv")
    df.to_csv(csv_out)



