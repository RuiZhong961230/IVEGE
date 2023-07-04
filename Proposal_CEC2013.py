import os
import numpy as np
from opfunu.cec_based import cec2013
import warnings
warnings.filterwarnings("ignore")


PopSize = 10                                                  # the number of individuals (PopSize > 4)
DimSize = 30                                                    # the number of variables
LB = [-100] * DimSize                                      # the maximum value of the variable range
UB = [100] * DimSize                                                 # the minimum value of the variable range
TrialRuns = 30                                                   # the number of independent runs
MaxFEs = 1000 * DimSize                      # the maximum number of fitness evaluations
GC = 6                                                                # the maximum growth cycle of an individual
GR = 1                                                                # the maximum growth radius of an individual
MS = 2                                                                # the moving scale
SeedNum = 6                                                          # the number of generated seeds by each individual

Pop = np.zeros((PopSize, DimSize))
PopFit = np.zeros(PopSize)                        # the fitness value of all individuals
curLife = 0                                                  # the current growth cycle of an individual
Func_num = 1                                                           # the serial number of benchmark function
curFEs = 0

SuiteName = "CEC2013"
def CheckIndi(Indi):
    for i in range(DimSize):
        range_width = UB[i] - LB[i]
        if Indi[i] > UB[i]:
            n = int((Indi[i] - UB[i]) / range_width)
            mirrorRange = (Indi[i] - UB[i]) - (n * range_width)
            Indi[i] = UB[i] - mirrorRange
        elif Indi[i] < LB[i]:
            n = int((LB[i] - Indi[i]) / range_width)
            mirrorRange = (LB[i] - Indi[i]) - (n * range_width)
            Indi[i] = LB[i] + mirrorRange
        else:
            pass


def Initialization(func):
    global Pop, PopFit, curLife, curFEs
    for i in range(PopSize):
        for j in range(DimSize):
            Pop[i][j] = np.random.uniform(LB[j], UB[j])
        PopFit[i] = func.evaluate(Pop[i])
    curLife = 1
    curFEs = PopSize


def Growth(func):
    global Pop, PopFit, curFEs
    offspring = np.zeros((PopSize, DimSize))
    offspring_fitness = np.zeros(PopSize)

    for i in range(PopSize):
        for j in range(DimSize):
            offspring[i][j] = Pop[i][j] + GR * (np.random.random() * 2.0 - 1.0)
        CheckIndi(offspring[i])
        offspring_fitness[i] = func.evaluate(offspring[i])
        curFEs += 1
        if offspring_fitness[i] < PopFit[i]:
            PopFit[i] = offspring_fitness[i]
            Pop[i] = offspring[i].copy()


def DynamicAllocation():
    global PopSize, PopFit
    fitness = np.zeros(PopSize)  # dynamic allocation of the computational budget
    cp = np.zeros(PopSize + 1)
    allocation = np.zeros(PopSize)
    for i in range(PopSize):
        allocation[i] = 3
    seed_FEs = SeedNum * PopSize - sum(allocation)
    for i in range(PopSize):
        fitness[i] = 1 / PopFit[i]
    fitness_sum = 0
    for f in fitness:
        fitness_sum += np.exp(f)
    c = 0
    for i in range(1, PopSize + 1):
        c += np.exp(fitness[i - 1]) / fitness_sum
        cp[i] = c
    for i in range(int(seed_FEs)):
        r = np.random.uniform(0, 1)
        for j in range(1, PopSize + 1):
            if cp[j - 1] <= r < cp[j]:
                allocation[j - 1] += 1
    return allocation


def Maturity(func):
    global Pop, PopFit, DimSize, curFEs
    seedIndi = np.zeros((PopSize * SeedNum, DimSize))
    seedFit = np.zeros(PopSize * SeedNum)

    allocation = DynamicAllocation()
    pointer = 0
    for i in range(PopSize):
        for j in range(int(allocation[i])):
            index1 = index2 = 0
            while index1 == i:
                index1 = np.random.randint(0, PopSize)
            while index2 == i or index2 == index1:
                index2 = np.random.randint(0, PopSize)

            seedIndi[pointer] = Pop[i] + MS * (np.random.random() * 2.0 - 1.0) * (
                        Pop[index1] - Pop[index2])

            r = np.random.uniform(0, 1)
            if r < 1 / 3:
                for m in range(DimSize):
                    if np.random.random() < 0.1:
                        seedIndi[pointer][m] += 0.05 * np.random.normal() * (UB[m] - LB[m])
            elif 1 / 3 <= r < 2 / 3:
                for m in range(DimSize):
                    if np.random.random() < 0.5:
                        seedIndi[pointer][m] = seedIndi[pointer][m]
                    else:
                        seedIndi[pointer][m] = Pop[i][m]
            else:
                for m in range(DimSize):
                    if np.random.random() < 0.01:
                        seedIndi[pointer][m] = np.random.uniform(LB[m], UB[m])
            CheckIndi(seedIndi[pointer])
            seedFit[pointer] = func.evaluate(seedIndi[pointer])
            curFEs += 1
            pointer += 1
    tmpIndi = np.vstack((Pop, seedIndi))
    tmpFit = np.hstack((PopFit, seedFit))

    tmp = list(map(list, zip(range(len(tmpFit)), tmpFit)))
    small = sorted(tmp, key=lambda x: x[1], reverse=False)
    for i in range(PopSize):
        key, _ = small[i]
        PopFit[i] = tmpFit[key]
        Pop[i] = tmpIndi[key].copy()


# the implementation process of differential evolution
def VEGE(func):
    global curLife, GC
    if curLife < GC:
        Growth(func)
        curLife += 1
    elif curLife == GC:
        Maturity(func)
        curLife = 0
    else:
        print("Error: Maximum generation period exceeded.")


def RunVEGE(func):
    global curFEs, Func_num, PopFit
    All_Trial_Best = []
    for i in range(TrialRuns):
        Best_list = []
        curFEs = 0
        np.random.seed(2022 + 88*i)
        Initialization(func)
        Best_list.append(min(PopFit))
        while curFEs < MaxFEs:
            VEGE(func)
            Best_list.append(min(PopFit))
        All_Trial_Best.append(Best_list)
    np.savetxt('./Proposal_Data/CEC2013/F{}_{}D.csv'.format(Func_num, DimSize), All_Trial_Best, delimiter=",")


def main(Dim):
    global Func_num, DimSize, MaxFEs, Pop, LB, UB, SuiteName
    DimSize = Dim
    Pop = np.zeros((PopSize, DimSize))
    MaxFEs = DimSize * 1000
    LB = [-100] * DimSize
    UB = [100] * DimSize
    SuiteName = "CEC2013"
    CEC2013Funcs = [cec2013.F12013(Dim), cec2013.F22013(Dim), cec2013.F32013(Dim), cec2013.F42013(Dim),
                    cec2013.F52013(Dim), cec2013.F62013(Dim), cec2013.F72013(Dim), cec2013.F82013(Dim),
                    cec2013.F92013(Dim), cec2013.F102013(Dim), cec2013.F112013(Dim), cec2013.F122013(Dim),
                    cec2013.F132013(Dim), cec2013.F142013(Dim), cec2013.F152013(Dim), cec2013.F162013(Dim),
                    cec2013.F172013(Dim), cec2013.F182013(Dim), cec2013.F192013(Dim), cec2013.F202013(Dim),
                    cec2013.F212013(Dim), cec2013.F222013(Dim), cec2013.F232013(Dim), cec2013.F242013(Dim),
                    cec2013.F252013(Dim), cec2013.F262013(Dim), cec2013.F272013(Dim), cec2013.F282013(Dim)]

    Func_num = 0
    for i in range(len(CEC2013Funcs)):
        Func_num = i + 1
        RunVEGE(CEC2013Funcs[i])


if __name__ == "__main__":
    if os.path.exists('./Proposal_Data/CEC2013') == False:   # Automatically create a folder when the folder does not exist
        os.makedirs('./Proposal_Data/CEC2013')
    Dims = [30, 50]
    for Dim in Dims:
        main(Dim)
