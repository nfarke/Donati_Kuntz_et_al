# Donati_Kuntz_et_al: Multi-omics analysis of CRISPRi-knockdowns identifies mechanisms that buffer decreases of enzymes in E. coli metabolism

These codes were used to create the model figures in the paper.

"create_SS_solutions_par.m" sets up the model structure and samples 1000 different parameter sets. The resulting stable parameter sets are then saved as PAR.mat

"Figure" replicates the model figures from the paper. It uses PAR.mat as input. The code is compatible with the parallel computing toolbox which speeds up the computation.

FigureS13_carAB and FigureS13_ppc create the supplemental figure S13
