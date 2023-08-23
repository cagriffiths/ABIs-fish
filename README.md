# ABIs-fish
Data and code to produce and visualise the age-based indicators (ABImsy and ABI0) in Griffiths et al. 2023. 

Data is split into three .rdata file:
(1) FLR_stock_objects - FLR stock assessment objects for 81 Category 1 stocks in the Northeast Atlantic
(2) age_structure_at_eq - estimated age structure at equilibrium for 72 stocks under both Fmsy and F0
(3) indicator_data - ABImsy and ABI0 values for each stock (n = 72) and year

Code is split into three primary files:
(1) age_structure_at_eq - code to approximate the age structure of a stock at equilibrium under a given F. Includes SRR specification and simulation from Griffiths et al. 2023. 
(2) calc_ABI - code to calculate ABImsy and ABI0 (as well as Amsy, Pmsy, A0 and P0). Includes threshold option (> 0.9) and can be used to extract total abundance of old fish. 
(3) plot_ABI - code to recreate Figures 3 and 4 from Griffiths et al. 2023 which contain all the information neccessary for stock status classification.

Simulation code will follow in the next few weeks. 

# FLR
The calculation of ABImsy and ABI0 requires several packages from the Fisheries Library in R (FLR) framework. 

To install these packages, start R and type:
  install.packages("FLCore", repos="http://flr-project.org/R")
  install.packages("FLBRP", repos="http://flr-project.org/R")
  install.packages("FLasher", repos="http://flr-project.org/R")

or directly directly from the github repository by using:
  remotes::install_github("flr/FLCore")
  remotes::install_github("flr/FLBRP")
  remotes::install_github("flr/FLasher")
  
For documentation or troubleshooting on FLR, we refer users to [https://github.com/flr]. The version of each package used is detailed at the start of each code file. 

# Citations
Griffiths, C. A., Winker, H., Bartolino, V., Wennhage, H., Orio, A. and Cardinale, M. (2023). Including older fish in fisheries management: a new age-based indicator and reference point for exploited fish stocks. Fish and Fisheries, accepted. 

Kell, L. T., I. Mosqueira, P. Grosjean, J-M. Fromentin, D. Garcia, R. Hillary, E. Jardim, S. Mardle, M. A. Pastoors, J. J. Poos, F. Scott, R. D. Scott. 2007. FLR: an open-source framework for the evaluation and development of management strategies. ICES J Mar Sci, 64 (4): 640-646. [doi: 10.1093/icesjms/fsm012].




