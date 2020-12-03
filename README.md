# co_culture_device
This repository contains code and data related to the BioMe device.


/modeling - This directory contains all of the information related to the computational modeling portion of the paper.

/modeling/analyses - This directory contains four different analyses that were performed in the computational modeling section of the paper. Within each analysis directory is a MATLAB script to run the analysis and several intermediate variables that can be loaded to speed up running/re-running the analysis. In the analyses folder titles, Equal_Leakage refers to inferring the leakage stoichiometries assuming that they have the same value for both amino acids and Unequal_Leakage infers these values separately. Lit_Params refers to assuming all other parameters from literature derived values and Dis_Params adds some noise to these parameters. These different analyses are described further in the paper.

/modeling/data - This directory contains the raw data for the syntrophic E. coli mutant co-culture experiment that was used for the computational modeling.

/modeling/model - This directory contains the code for the computational model used in this work, written as a MATLAB function.
