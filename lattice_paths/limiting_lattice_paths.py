import numpy as np
import matplotlib.pyplot as plt
import random
import time
import tqdm
from collections import Counter

from lattice_utils import *


# SLOW VERSIONS
# startingState = maxState(300)
# animateRandomWalk(startingState, getUniformLegalContractedStateEnumerative, 100000, ms_per_draw=1, max_updates_per_draw=None)

# startingState = minState(300)
# animateRandomWalk(startingState, getUniformLegalExpandedStateEnumerative, 100000, ms_per_draw=1, max_updates_per_draw=None)

# startingState = minState(100) # NOTE: not actually fully uniform right now.
# animateRandomWalk(startingState, getUniformLegalNextStateEnumerative, 100000, ms_per_draw=1, max_updates_per_draw=None)

# FAST VERSIONS
# startingState = maxState(300)
# animateRandomWalk(startingState, getUniformLegalContractedStateFast, 100000, ms_per_draw=1, max_updates_per_draw=None)

# startingState = minState(300)
# animateRandomWalk(startingState, getUniformLegalExpandedStateFast, 100000, ms_per_draw=1, max_updates_per_draw=None)

startingState = minState(100) # NOTE: not actually fully uniform right now.
animateRandomWalk(startingState, getUniformLegalNextStateFast, 100000, ms_per_draw=1, max_updates_per_draw=None)

# add weights, prob of going up neq going down





# state = maxState(5)
# legalContractions = getLegalContractions(state)
# counts = Counter([getBinaryRepr(getContractionNextState(state, *c)) for c in legalContractions])
# for i in range(1000000):
#     fast_sample = getUniformLegalContractedStateFast(state)
#     # slow_sample = getUniformLegalContractedStateEnumerative(state)
#     counts[getBinaryRepr(fast_sample)] += 1
#     # counts[getBinaryRepr(slow_sample)] += 1

# for key, value in counts.items():
#     print(key, value)
