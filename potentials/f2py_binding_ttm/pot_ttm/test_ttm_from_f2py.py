import os
os.environ['LD_LIBRARY_PATH'] = os.getcwd()
# this sets the library path to looks for the ttm dll in the current directory

import numpy as np
import ttm
print(ttm.ttm_from_f2py.__doc__)

pos = np.array(([1.260474, -1.006574, -0.094009],
                [1.222678, -0.032592, -0.011174],
                [1.91971,  -1.290587,  0.549949]))
print(pos.T)
derivatives = np.zeros(pos.shape)
final_energy = 0.0
ttm_model = 21

print(ttm.ttm_from_f2py(ttm_model, pos.T, 1))
