import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import random
import time
import tqdm

'''
'index' = which of the n balls 
'slot' = the position of a ball or space in the 2*n-1 slots of the x-ray vector
'''
def getBinaryRepr(state, asString=True):
    '''
    state: an array of length n of the indices of the n balls
    asString: if True, return a string of 0s and 1s, else return a list of integers
    '''
    n = len(state)
    binary_repr = np.zeros(2*n-1, dtype=int)
    for i in state:
        binary_repr[i-1] = 1
    if asString:
        return ''.join([str(x) for x in binary_repr])
    return binary_repr

def maxState(n):
    '''
    state of maximum moment of inertia (all balls spread out / anti-diagonal matrix)
    No legal expansions moves from this state.
    '''
    return [2*i+1 for i in range(n)]

def minState(n):
    '''
    state of minimum moment (all balls contracted)
    No legal contraction moves from this state.
    '''
    if n % 2 == 0:
        # for n even, one 0 in the middle, all ones around it
        return [int(n/2) + i for i in range(int(n/2))] + [n + 1 + i for i in range(int(n/2))]
    else:
        # for n odd, all 1s are in the middle
        return [int(n/2) + 1 + i for i in range(n)]

def spaceToRight(state, i_ind):
    '''
    Given a state and ball with index i, check if there is a space immediately to the right.
    i.e. ball can plausibly be the left of a contraction or right of an expansion.
    '''
    n = len(state)
    i_slot = state[i_ind]
    if i_slot == 2*n - 1: # last slot
        return False 
    if i_ind < n-1: # not the last ball
        if state[i_ind+1] != i_slot+1: # space right
            return True
        return False # ball to right
    return True # last ball but not last slot

def spaceToLeft(state, i_ind):
    '''
    Given a state and ball with index i, check if there is a space immediately to the left.
    i.e. ball can plausibly be the right of a contraction or left of an expansion.
    '''
    n = len(state)
    i_slot = state[i_ind]
    if i_slot == 1: # first slot
        return False 
    if i_ind > 0: # not the first ball
        if state[i_ind-1] != i_slot-1: # space left
            return True
        return False # ball to left
    return True # first ball but not first slot

def dominanceCondition(state):
    '''
    state: an array of length n of the slots of the n balls

    check that the sum of the slot numbers of the first k balls is at least k^2,
    i.e. the signed area above the diagonal is always positive moving left to right.
    NOT EFFICIENT. Enumerating possible next states is rarely necessary.
    '''
    running_sum = 0
    for ball_ind, ball_slot in enumerate(state):
        running_sum += ball_slot
        if running_sum < (ball_ind+1)**2:
            return False
    return True

def getAreasAtBallIndices(state):
    '''
    return an array of the signed area above the diagonal at each ball index (including the ball). 
    e.g. if the area above the diagonal after the ball 3 is 2, then the area[3] = 2.
    If this area is positive at index i_ind, then the ball i_ind can be moved left without violating dominance.
    Areas should be all nonnegative for a legal state.
    '''
    areas = []
    running_sum = 0
    for ball_ind, ball_slot in enumerate(state):
        running_sum += ball_slot
        areas.append(running_sum - (ball_ind+1)**2)
    return areas

def getAreasAtSlots(state):
    '''
    the signed area above the diagonal at every slot, not just ball indices.
    Include 0 at the '0th slot', so return a list of length n+1.
    '''
    binary_repr = getBinaryRepr(state, False)
    recorded_areas = [0]
    curr_area = 0
    curr_y = 0
    diag_y = 0
    for digit in binary_repr:
        if digit == 1:
            curr_area += curr_y - diag_y
            recorded_areas.append(curr_area)
        curr_y += 1 - digit
        diag_y += digit
    return recorded_areas

def isLegalState(state):
    '''
    return true if a state is legal. 
    check number of balls, center of mass, and dominance condition.
    '''
    n = len(state)
    if state[-1] > 2*n-1:
        return False
    COM = sum(state) / n
    if COM != n:
        return False
    return dominanceCondition(state)


def preconditionDominanceCheck(state):
    '''
    return a function which takes two ball indices and returns True if the dominance condition is satisfied after expanding at those indices.
    More efficient than checking the dominance condition after every expansion.
    Use a union-find data struct to keep track of the segments of nonzero areas.
    '''
    n = len(state)
    areas = getAreasAtBallIndices(state) # calculate the areas above the diagonal at each ball index once. Use later to check dominance for expansion moves.
    # iterate over area values and label contiguous segments of nonzero areas. 
    segments = {0: 0} # map each ball index to its group leader index, or None if it is a zero area.
    for i_ind in range(1, n):
        if areas[i_ind] > 0:
            if areas[i_ind-1] > 0: # continue the segment.
                segments[i_ind] = segments[i_ind-1]
            else: # start a new segment.
                segments[i_ind] = i_ind
        else:
            segments[i_ind] = None
    
    def obeysDominance(i_ind, j_ind):
        # if the area above the diagonal is positive for all balls in [i_ind, j_ind), the left ball can be moved left without violating dominance 
        # check if i and j-1 are in the same segment. If not, then there is a 0 area between them.
        return segments[i_ind] != None and segments[j_ind-1] == segments[i_ind]
    
    return obeysDominance
            
def getLegalContractions(state):
    '''
    state: an array of length n of the indices of the n balls
    return a list of pairs of indices of balls which can legally both contract
    '''
    n = len(state)

    spaceLeftBalls = [i_ind for i_ind in range(n) if spaceToLeft(state, i_ind)]
    spaceRightBalls = [i_ind for i_ind in range(n) if spaceToRight(state, i_ind)]

    legal_contracting_moves = []
    for r_i in range(len(spaceRightBalls)):
        i_ind = spaceRightBalls[r_i]
        for j_ind in spaceLeftBalls[max(0, r_i-1):]: # take only the balls to the right of ball i
            if i_ind >= j_ind:
                continue
            i_slot, j_slot = state[i_ind], state[j_ind]
            if j_slot - i_slot <= 2: # contraction requires at least 2 slots of distance
                continue
            legal_contracting_moves.append((i_ind, j_ind))

    return legal_contracting_moves

def getLegalExpansions(state):
    '''
    state: an array of length n of the indices of the n balls
    return a list of pairs of indices of balls which can legally both expand
    '''
    n = len(state)

    spaceLeftBalls = [i_ind for i_ind in range(n) if spaceToLeft(state, i_ind)]
    spaceRightBalls = [i_ind for i_ind in range(n) if spaceToRight(state, i_ind)]
    
    obeysDominance = preconditionDominanceCheck(state) # gets a preconditioned function which checks dominance after expanding at two indices.

    legal_expanding_moves = []
    for l_i in range(len(spaceLeftBalls)):
        i_ind = spaceLeftBalls[l_i]
        for j_ind in spaceRightBalls[max(0, l_i-1):]: # take only the balls to the right of ball i
            if i_ind >= j_ind:
                continue
            if not obeysDominance(i_ind, j_ind):
                continue
            legal_expanding_moves.append((i_ind, j_ind))

    return legal_expanding_moves

def getContractionNextState(state, i_ind, j_ind):
    i = state[i_ind]
    j = state[j_ind]
    new_state = state.copy()
    new_state[i_ind] = i+1
    new_state[j_ind] = j-1
    return new_state

def getExpansionNextState(state, i_ind, j_ind):
    i = state[i_ind]
    j = state[j_ind]
    new_state = state.copy()
    new_state[i_ind] = i-1
    new_state[j_ind] = j+1
    return new_state

def getLegalNextStates(state):
    '''
    state: an array of length n of the indices of the n balls

    return a list of all possible next states from the current state.
    NOT EFFICIENT. Enumerating possible next states is rarely necessary.
    '''
    contractions, expansions = getLegalContractions(state), getLegalExpansions(state)
    legal_next_states = []
    for (i_ind, j_ind) in contractions:
        new_state = getContractionNextState(state, i_ind, j_ind)
        legal_next_states.append(new_state)

    for (i_ind, j_ind) in expansions:
        new_state = getExpansionNextState(state, i_ind, j_ind)
        legal_next_states.append(new_state)
        
    return legal_next_states

def getUniformLegalNextStateEnumerative(state):
    '''
    Sample a random legal move from the state by checking all legal moves.
    return the new state.

    If no legal moves, return None.
    '''
    contractions, expansions = getLegalContractions(state), getLegalExpansions(state)
    num_contractions = len(contractions)
    num_expansions = len(expansions)
    if num_contractions + num_expansions == 0:
        return None
    rand_int = random.randint(0, num_contractions + num_expansions - 1)
    if rand_int < num_contractions:
        i_ind, j_ind = contractions[rand_int]
        return getContractionNextState(state, i_ind, j_ind)
    else:
        i_ind, j_ind = expansions[rand_int - num_contractions]
        return getExpansionNextState(state, i_ind, j_ind)

def getUniformLegalNextStateFast(state, tries=100): 
    '''
    Sample a random legal move from the state without checking all legal moves.
    Find the balls with right and left spaces, but instead of enumerating possible moves, sample a random possible pair.

    If all tries fail, find all possible legal moves and sample one.
    
    return the new state.

    Let C, E = set of legal contractions and expansions and E' = set of expansions which are almost legal but disobey dominance.
    Assume |C| = |E'|. Randomly sample a move from C + E'. 
    If the move is from C, return the contraction. 
    If the move is from E', check dominance and return the expansion if legal.
    Else, add the move to failed_expansions and retry.
    '''
    n = len(state)

    # raise NotImplementedError("Fix uniformity of choices.")
    # TODO: NOT ACTUALLY A UNIFORM NEXT STATE SAMPLER. JUST A FAST SAMPLER.

    if random.random() > 0.5:
        # sample a contraction
        return getUniformLegalContractedStateFast(state, tries)
    else:
        # sample an expansion
        return getUniformLegalExpandedStateFast(state, tries)
    
    # screen initially for balls with spaces on the right or left side.
    # spaceLeftBalls = [i_ind for i_ind in range(n) if spaceToLeft(state, i_ind)]
    # spaceRightBalls = [i_ind for i_ind in range(n) if spaceToRight(state, i_ind)]

    # obeysDominance = preconditionDominanceCheck(state) # gets a preconditioned function which checks dominance after expanding at two indices.
    # for _ in range(tries):
    #     if random.random() < 0.5:
    #         # sample a contraction
    #         r_i = random.randint(0, len(spaceRightBalls)-1)
    #         i_ind = spaceRightBalls[r_i]
    #         j_ind = random.choice(spaceLeftBalls[max(0, r_i-1):])
    #         if i_ind >= j_ind:
    #             continue
    #         i_slot, j_slot = state[i_ind], state[j_ind]
    #         if j_slot - i_slot <= 2:
    #             continue
    #         return getContractionNextState(state, i_ind, j_ind)
    #     else:
    #         # sample an expansion
    #         l_i = random.randint(0, len(spaceLeftBalls)-1)
    #         i_ind = spaceLeftBalls[l_i]
    #         j_ind = random.choice(spaceRightBalls[max(0, l_i-1):])
    #         if i_ind >= j_ind:
    #             continue
    #         if not obeysDominance(i_ind, j_ind):
    #             continue
    #         return getExpansionNextState(state, i_ind, j_ind)
        
    # # if all tries fail, find all possible legal moves and sample one.
    # return getUniformLegalNextStateEnumerative(state)

def getUniformLegalContractedStateEnumerative(state):
    '''
    Sample a random legal contraction from the state by checking all legal contractions.
    return the new state.

    If no legal contractions, return None.
    '''
    contractions = getLegalContractions(state)
    if not contractions:
        return None
    i_ind, j_ind = random.choice(contractions)
    return getContractionNextState(state, i_ind, j_ind)

def getUniformLegalContractedStateFast(state, tries=100):
    '''
    Sample a random legal contraction from the state without checking all legal contractions.
    Find the balls with right and left spaces, but instead of enumerating possible moves, sample a random possible pair.

    If all tries fail, find all possible legal contractions and sample one.
    
    return the new state.

    If no legal contractions, return None.
    '''
    n = len(state)

    # screen initially for balls with spaces on the right or left side.
    spaceLeftBalls = [i_ind for i_ind in range(n) if spaceToLeft(state, i_ind)]
    spaceRightBalls = [i_ind for i_ind in range(n) if spaceToRight(state, i_ind)]

    # weight the left balls by the number of possible right balls to contract with.
    # this isnt perfect - we include some impossible contractions, but we just trash those and try again.
    r_i_weights = [len(spaceLeftBalls)-r_i for r_i in range(len(spaceRightBalls))]
    s = sum(r_i_weights)
    r_i_weights = [w/s for w in r_i_weights]
    for _ in range(tries):
        r_i = np.random.choice(range(len(spaceRightBalls)), p=r_i_weights)
        i_ind = spaceRightBalls[r_i]
        l_i = random.randint(r_i, len(spaceLeftBalls)-1)
        j_ind = spaceLeftBalls[l_i]
        if i_ind >= j_ind:
            continue
        i_slot, j_slot = state[i_ind], state[j_ind]
        if j_slot - i_slot <= 2:
            continue
        return getContractionNextState(state, i_ind, j_ind)
        
    # if all tries fail, find all possible legal contractions and sample one.
    return getUniformLegalContractedStateEnumerative(state)

def getUniformLegalExpandedStateEnumerative(state):
    '''
    Sample a random legal expansion from the state by checking all legal expansions.
    return the new state.
    If no legal Expansions: return None
    '''
    expansions = getLegalExpansions(state)
    if not expansions:
        return None
    i_ind, j_ind = random.choice(expansions)
    return getExpansionNextState(state, i_ind, j_ind)

def getUniformLegalExpandedStateFast(state, tries=100):
    '''
    Sample a random legal expansion from the state without checking all legal expansions.
    Find the balls with right and left spaces, but instead of enumerating possible moves, sample a random possible pair.

    If all tries fail, find all possible legal expansions and sample one.
    
    return the new state.
    '''
    n = len(state)
    spaceLeftBalls = [i_ind for i_ind in range(n) if spaceToLeft(state, i_ind)]
    spaceRightBalls = [i_ind for i_ind in range(n) if spaceToRight(state, i_ind)]
    obeysDominance = preconditionDominanceCheck(state)

    # weight the left balls by the number of possible right balls to expand with.
    # this isn't perfect - we include some impossible expansions, but we just trash those and try again.
    l_i_weights = [len(spaceRightBalls)-l_i for l_i in range(len(spaceLeftBalls))]
    s = sum(l_i_weights)
    l_i_weights = [w/s for w in l_i_weights]

    dominance_failures = 0
    index_failures = 0

    for _ in range(tries):
        l_i = np.random.choice(range(len(spaceLeftBalls)), p=l_i_weights)
        i_ind = spaceLeftBalls[l_i]
        r_i = random.randint(l_i, len(spaceRightBalls)-1)
        j_ind = spaceRightBalls[r_i]
        if i_ind >= j_ind:
            index_failures += 1
            continue
        if not obeysDominance(i_ind, j_ind):
            dominance_failures += 1
            continue
        # print(f"Index failures: {index_failures}, Dominance failures: {dominance_failures}")
        return getExpansionNextState(state, i_ind, j_ind)
        
    # if all tries fail, find all possible legal expansions and sample one.
    return getUniformLegalExpandedStateEnumerative(state)

def getPlotPoints(state):
    '''
    return the x and y coordinates of the lattice paths of a state
    '''
    binary_repr = getBinaryRepr(state, False)
    x_vals = [0]
    y_vals = [0]
    # each 1 is a horizontal step, each 0 is a vertical step.
    x = 0
    y = 0
    for step in binary_repr:
        x += step
        y += 1-step
        x_vals.append(x)
        y_vals.append(y)
    return x_vals, y_vals

def plotState(state, showDiag=False, plotArea=True):
    '''
    draw the lattice paths of a state
    '''
    n = len(state)

    if plotArea:
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 5))
    else:
        fig, ax1 = plt.subplots(1, 1, figsize=(5, 5))

    ax1.plot(*getPlotPoints(state))

    if showDiag:
        # plot the diagonal path
        ax1.plot(*getPlotPoints(maxState(n)), color='gray')

    if plotArea:
        ax2.plot(range(n+1), getAreasAtSlots(state))    

    plt.show()

def plotRandomWalk(startingState, nextStateFunc, iters, plotInterval):
    '''
    From a starting state, do a random walk by choosing a state using the next state function.
    plot the state at regular intervals.
    '''
    state = startingState
    for i in range(iters // plotInterval):
        startTime = time.time()
        for _ in tqdm.tqdm(range(plotInterval)):
            state = nextStateFunc(state)
        print(f"{plotInterval} iterations took {time.time() - startTime}")
        plotState(state, showDiag=True, plotArea=True)

def animateRandomWalk(startingState, nextStateFunc, iters, max_updates_per_draw=None, ms_per_draw=1):
    '''
    From a starting state, do a random walk by choosing a state using the next state function.
    animate the state at each step along with the net area above the diagonal.

    max_updates_per_draw: the maximum number of updates to make in a single draw.
    set to None for no limit.
    ms_per_draw: the time in milliseconds to spend drawing the state.
    '''
    if ms_per_draw < 1:
        raise ValueError("ms_per_draw must be at least 1")
    if max_updates_per_draw is None:
        max_updates_per_draw = iters

    def nestedNextStateFunc(state):
        start_time = time.time()
        updates_this_frame = 0
        while time.time() - start_time < ms_per_draw / 1000 and updates_this_frame < max_updates_per_draw:
            updates_this_frame += 1
            nextState = nextStateFunc(state)
            if nextState is None: # there are no legal moves from this state
                return state, 0
            state = nextState
        return state, updates_this_frame
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 5))
    UL = UpdateLattice(ax1, ax2, startingState, nextStateFunc=nestedNextStateFunc, n=len(startingState))
    latticeAnim = FuncAnimation(fig, UL.update_lattice, init_func=UL.start_lattice, frames=iters, interval=ms_per_draw, blit=True)
    areasAnim = FuncAnimation(fig, UL.update_areas, init_func=UL.start_areas, frames=iters, interval=ms_per_draw, blit=True)
    plt.show()

class UpdateLattice:
    def __init__(self, lattice_ax, area_ax, startingState, nextStateFunc, n=100, diag=True):
        self.state = startingState
        self.nextStateFunc = nextStateFunc
        self.n = n

        self.lattice_ax = lattice_ax
        self.diag, = lattice_ax.plot(*getPlotPoints(maxState(n)), '-', color='gray')
        self.path, = lattice_ax.plot(*getPlotPoints(startingState), 'k-')

        self.lattice_ax.set_xlim(0, n)
        self.lattice_ax.set_ylim(0, n)
        self.lattice_ax.grid(True)
        self.updates_per_draw_recent = [0]*100

        self.area_ax = area_ax
        self.area_ax.set_xlim(0, n)
        self.area_ax.set_ylim(0, n)
        self.area_ax.grid(True)
        self.area_ax.set_title("Signed area above diagonal")
        self.area_line, = area_ax.plot([], [], 'k-')

    def start_lattice(self):
        return self.path,

    def start_areas(self):
        return self.area_line,
        
    def update_lattice(self, i):
        nextState, updates_done = self.nextStateFunc(self.state)
        # num_contractions = len(getLegalContractions(self.state))
        # percent_moves_contract = num_contractions / (len(getLegalExpansions(self.state)) + num_contractions)
        if nextState is not None:
            self.state = nextState
        self.path.set_data(*getPlotPoints(self.state))

        self.updates_per_draw_recent.pop(0)
        self.updates_per_draw_recent.append(updates_done)

        # self.lattice_ax.set_xlabel(f"{np.mean(self.updates_per_draw_recent)} updates per frame")
        if i % 100 == 0:
            # print(f"Percent moves contractions: {percent_moves_contract}")

            if np.mean(self.updates_per_draw_recent) > 0:
                print(f"{np.mean(self.updates_per_draw_recent)} updates per frame")

        return self.path,

    def update_areas(self, i):
        areas = getAreasAtSlots(self.state)
        self.area_line.set_data(range(self.n+1), areas)

        self.area_ax.set_ylim(0, self.n**2/10)

        return self.area_line,




