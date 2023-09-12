# ABIs-fish
Data and code to calculate the age-based indicators (ABImsy and ABI0) in Griffiths et al. 2023. 

*Data* is split into three .rdata files:

(1) `FLR_stock_objects` - FLR stock objects for 81 ICES Category 1 stocks in the Northeast Atlantic. Includes all aspects of each stock's age based analytical assessment including estimates of SSB, R and F, and details of any reference points (e.g., Fmsy or MSY Btigger) and methods (e.g., SAM or SS3) used. 

(2) `age_structure_at_eq` - estimated age structure at equilibrium for 72 stocks under both Fmsy and F0.

(3) `indicator_data` - ABImsy and ABI0 values for each stock (n = 72) and year. Seperate .csv files for each stock are also provided in the folder `indicator data by stock`. 

A fourth .rdata file `sim6stks` is generated by the code `create_stock_OMs_for_simulations` and is used in the code `simulation_evaluations_and_ROCs`, however, the file is too large [470MB] to upload to GitHub and therefore we request that users run the OM creation code and generate a local copy of the .rdata file. `sim6stks` contains all the simulation output for 6 selected stocks as per Griffiths et al. 2023.  

*Code* is split into five primary files:

(1) `age_structure_at_eq` - code to approximate the age structure of a stock at equilibrium under a given F. Includes SRR specification and simulation from Griffiths et al. 2023. Also produces the example .csv file `output/age_structure_Eq_hke.27.3a46-8abd.csv` and the plot `output/age_structure_Eq_hke.27.3a46-8abd`.  

(2) `calc_ABI` - code to calculate ABImsy and ABI0 (as well as Amsy, Pmsy, A0 and P0). Includes threshold option (> 0.9) and can be used to extract total abundance of old fish. Also produces the files listed in `data/indicator data by stock`. 

(3) `plot_ABI` - code to recreate Figures 3 and 4 from Griffiths et al. 2023 which contain all the information neccessary for stock status classification. Produces the plot `output/stock_plot_her.27.5a`.

(4) `create_stock_OMs_for_simulations` - code to create operating models and run simulation tests for 6 selected stocks. Produces the plot `output/sim_plot_pra.27.3a4a`

(5) ...

# FLR
The calculation of ABImsy and ABI0 requires several packages from the Fisheries Library in R (FLR) framework. 

To install these packages, start R and type:
 ``` 
  install.packages("FLCore", repos="http://flr-project.org/R")
  
  install.packages("FLBRP", repos="http://flr-project.org/R")
  
  install.packages("FLasher", repos="http://flr-project.org/R")
 ``` 
or directly directly from the github repository by using:
 ``` 
  remotes::install_github("flr/FLCore")
  
  remotes::install_github("flr/FLBRP")
  
  remotes::install_github("flr/FLasher")
 ``` 
For documentation or troubleshooting on FLR, we refer users to https://github.com/flr. The version of each package used is detailed at the start of each `code` file. 

# Citations
Griffiths, C. A., Winker, H., Bartolino, V., Wennhage, H., Orio, A. and Cardinale, M. (2023). Including older fish in fisheries management: a new age-based indicator and reference point for exploited fish stocks. *Fish and Fisheries*, accepted. 

Kell, L. T., I. Mosqueira, P. Grosjean, J-M. Fromentin, D. Garcia, R. Hillary, E. Jardim, S. Mardle, M. A. Pastoors, J. J. Poos, F. Scott, R. D. Scott. 2007. FLR: an open-source framework for the evaluation and development of management strategies. *ICES Journal of Marine Science*, 64 (4): 640-646. doi: 10.1093/icesjms/fsm012.




