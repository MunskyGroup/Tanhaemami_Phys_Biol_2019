### Codes
In this folder, 2 scripts can be run by the user:
1. **BODIPY\_signal\_histPlot\_ColorCoded.m**:  
Run this script to see the BODIPY signals collected by the flow cytometer.
2. **RunFCSC\_3M\_Log.m (THE MAIN SCRIPT)**:  
Run this script (which automatically uses the *FCSC\_3M\_Log.m* file as its object) to analyze our main flow cytometry data and see the results at different computation steps. This script provides a weighted model at its last step in our strategy.

### The Directories You See Here
After you run the main script as instructed above (number 2), you can go to the following directories to discover more and see the results:  
- **Codes\_for\_figures**  
Go to this directory to generate the figures (as you see in the paper and the supplementary materials).  
- **Codes\_for\_NewAnalyses**  
Go to this directory to test our optimized weighted modeling stratey on a new FCM data set.  
- **Codes\_for\_regular\_approach**  
Go to this directory to see how a model based on averaged information would look like. We called this the regular approach, in which the model over predicts for lower, and under predicts for higher signal intensities.  
- **Codes\_to\_load\_and\_process\_Main_Data**  
If you would like to see how we called the *\*.csv* files from the flow cytometer and gerenated a MATLAB structure in the format of a *.mat file, go to this directory and follow the instructions.  
- **Functions**  
This directory contains all the functions needed for data analysis and modeling, which will be automatically called and executed upon running the main script (RunFCSC\_3M\_Log.m).  
- **Thomas\_Blasi\_Cell\_Cycle\_Analysis**  
This directory contains the codes for analysis of the GBML method, provided by Blasi et al (for publications references, please see the paper and the supplementary information).

