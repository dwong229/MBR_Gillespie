READ ME for denwong_BE567_final_project submission, Dec 5 2012 
Directions to run chemotaxis model for e.coli in different scenarios.  

updated: 1/13/2014: section 2.5  

- MAIN FUNCTION TO RUN -
EcoliSimMain.m : main function to run simulation.  Relies on all other functions in folder.  Set 'SimMode' variable to change the different cases of the simulation.  
1: Single cell in free solution
2: 'simIterations' times simulations of 1
3: Single cell in chemoattractant environment
4: 'simIterations' times simultion of 3.
5: MBR with 4 cells

- 2 SUBFUNCTIONS -
----- 2.1 Gillespie functions:
ecoli_gillespie_func.m : Gillespie simulation of reactions for e.coli
cell_gillespie.m : Same equations a in ecoli_gillespie_func.m but with different inputs and outputs suitable for MBR case.
MBR_gillespie_func.m : keeps track of state of MBR, calls a gillespie function within to compute which bacterium has a reaction.

----- 2.2 chemo-attractant chemical gradient
attract_lin.m  : for linear deacy of concentration from source   
attract_exp.m  : for case of exponential decaying concentration from source
concgradientvis.m : visualize the concentration gradient on trajectory plot

----- 2.3 data visualizaiton
plotStates.m   : makes plots state for 1 iteration of chemical levels
runtumble.m    : compute run and tumble ratios for single run models  

----- 2.4 for run time visualization (from CIS 520)
CTimeleft.m 
sec2timestr.m       

----- 2.5 for converting output from OpenCV
cellposn_ends2headangle.m  :  takes in headtail coordinates from OpenCV contour finder to head-angle notation for simulation
