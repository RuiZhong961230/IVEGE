# IVEGE
Vegetation Evolution with Dynamic Maturity Strategy and Diverse Mutation Strategy for Solving Optimization Problems

## Abstract
We introduce two new search strategies to further improve the performance of vegetation evolution (VEGE) for solving continuous optimization problems. Specifically, the first strategy, named the dynamic maturity strategy, allows individuals with better fitness to have a higher probability of generating more seed individuals. Here, all individuals will first become allocated to generate a fixed number of seeds, and then the remaining number of allocatable seeds will be distributed competitively according to their fitness. Since VEGE performs poorly in getting rid of local optima, we propose the diverse mutation strategy as the second search operator with several different mutation methods to increase the diversity of seed individuals. In other words, each generated seed individual will randomly choose one of the methods to mutate with a lower probability. To evaluate the performances of the two proposed strategies, we run our proposal (VEGE + two strategies), VEGE, and another seven advanced evolutionary algorithms (EAs) on the CEC2013 benchmark functions and seven popular engineering problems. Finally, we analyze the respective contributions of these two strategies to VEGE. The experimental and statistical results confirmed that our proposal can significantly accelerate convergence and improve the convergence accuracy of the conventional VEGE in most optimization problems.

## Citation
@Article{Zhong:23,  
AUTHOR = {Rui Zhong and Fei Peng and Enzhi Zhang and Jun Yu and Masaharu Munetomo},  
TITLE = {Vegetation Evolution with Dynamic Maturity Strategy and Diverse Mutation Strategy for Solving Optimization Problems},  
JOURNAL = {Biomimetics},  
VOLUME = {8},  
YEAR = {2023},   
NUMBER = {6},  
ARTICLE-NUMBER = {454},  
ISSN = {2313-7673},  
DOI = {10.3390/biomimetics8060454 }  
}

## Datasets and Libraries
CEC benchmarks are provided by the opfunu library and engineering problems are provided by the enoppy library.
