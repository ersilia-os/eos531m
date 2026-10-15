# imports
import os
import sys

import numpy as np
from ersilia_pack_utils.core import read_smiles, write_out

# make the wrapper importable
root = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, root)

from predhlm_predict import COLUMNS, predict  # noqa: E402

# parse arguments
input_file = sys.argv[1]
output_file = sys.argv[2]


# my model: SMILES -> PredHLM human liver microsome half-life (log10 min), SD, minutes, SDC
def my_model(smiles_list):
    return predict(smiles_list)


# read SMILES from .csv file, assuming one column with header
_, smiles_list = read_smiles(input_file)

# run model
outputs = my_model(smiles_list)

# check input and output have the same length
assert len(smiles_list) == len(outputs)

# write output in a .csv file
write_out(outputs, COLUMNS, output_file, np.float32)
