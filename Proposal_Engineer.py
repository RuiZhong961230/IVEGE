import os
from enoppy.paper_based.pdo_2022 import *
import warnings
warnings.filterwarnings("ignore")


PopSize = 10
DimSize = 30
LB = [-100] * DimSize
UB = [100] * DimSize
TrialRuns = 30
MaxFEs = 20000
GC = 6
GR = 1
MS = 2
SeedNum = 6

Pop = np.zeros((PopSize, DimSize))
PopFit = np.zeros(PopSize)
curLife = 0
Func_num = 1
curFEs = 0


SuiteName = "Engineer"
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
    np.savetxt('./Proposal_Data/Engineer/{}.csv'.format(Func_num), All_Trial_Best, delimiter=",")


def main():
    global Func_num, DimSize, Pop, MaxFEs, SuiteName, LB, UB

    SuiteName = "Engineer"
    MaxFEs = 20000

    # Probs = [WBP(), PVP(), CSP(), SRD(), TBTD(), GTD(), CBD(), IBD(), TCD(), PLD(), CBHD(), RCB()]
    # Names = ["WBP", "PVP", "CSP", "SRD", "TBTD", "GTD", "CBD", "IBD", "TCD", "PLD", "CBHD", "RCB"]

    Probs = [GTD()]
    Names = ["GTD"]

    for i in range(len(Probs)):
        DimSize = Probs[i].n_dims
        LB = np.array(Probs[i].bounds)[:, 0]
        UB = np.array(Probs[i].bounds)[:, 1]

        Pop = np.zeros((PopSize, DimSize))
        Func_num = Names[i]
        RunVEGE(Probs[i])


if __name__ == "__main__":
    if os.path.exists('./Proposal_Data/Engineer') == False:
        os.makedirs('./Proposal_Data/Engineer')
    main()
