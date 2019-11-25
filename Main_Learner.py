import Analyze_structure
import sys
import numpy as np

data = sys.argv[1]
test = sys.argv[2]
e,tau = Analyze_structure.analyze(data, test)

np.savetxt("emission.txt",e)
np.savetxt("transition.txt",tau)