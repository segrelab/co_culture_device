# co_culture_device
This repository contains code and data related to the BioMe device.

/BioMe - Distribution Files - This directory contains relevant files for the development of the BioMe device.

/BioMe - Distribution Files/Drawings - This directory contains PDF files of the technical drawings for the BioMe Base, Body, and Gasket parts.

/BioMe - Distribution Files/Fabrication Files - This directory contains files for 3D printing of the Body and laser cutting of the Gasket parts.

/BioMe - Distribution Files/Solidworks Files - This directory contains the SolidWorks computer aided drafting (CAD) files for the BioMe Base, Body, and Gasket parts.

/modeling - This directory contains all of the information related to the computational modeling portion of the paper.

/modeling/analyses - This directory contains four different analyses that were performed in the computational modeling section of the paper. Within each analysis directory is a MATLAB script to run the analysis and several intermediate variables that can be loaded to speed up running/re-running the analysis. In the analyses folder titles, Equal_Leakage refers to inferring the leakage stoichiometries assuming that they have the same value for both amino acids and Unequal_Leakage infers these values separately. Lit_Params refers to assuming all other parameters from literature derived values and Dis_Params adds some noise to these parameters. These different analyses are described further in the paper.

/modeling/model - This directory contains the code for the computational model used in this work, written as a MATLAB function.

/raw_data_and_plots - This directory contains all of the raw data and MATLAB scripts for generating all of the plots in the figures.
