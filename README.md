# ThesisProject
In my thesis, we've studied star systems, and the relation between the orbital period (PERIOD) and rotational velocity at the equator (VBROAD). We've attempted to find the linear relation between the two variables, but our population depicted a very larger scatter in the data, as seen in 'Figures/scatter without solution.jpg'. Common robust method techinques were ineffectual in determining the actual slope of the densest are in the graph, hence we developed a novel MCMC algorithm which fits a 2D distribution with a 2D Gaussian to derive the slope of the data. 

A comprehensive descroption of the algorithm can be found in this article: 'https://ui.adsabs.harvard.edu/abs/2025A%26A...701A.195H/abstract'. 

Following is a descrioption of each file in the project: 
1. 'EBs with vbroad on MS.csv' is the data file. It includes information about 1050 star systems, including the rotational velocity (VBROAD) & orbital periods (PERIOD), which are the main focus.
2. 'global_variable.py' includes general functions & variables to create graphs, folders, etc., used throught the project. 
3. 'Gaussin_Function_Fitting.py' is the algorithm's code. It accepts two data vectors (PERIOD and VBROAD) and fits them with a 2D Gaussian using MCMC method. The inclination of the Gaussian represents the linear relation between the variables. This file also generates the relevant graphs (distribution of the MCMC parameters chain, corner plot, etc.). These graphs are available in 'Figures/depth ratio below 0.7 temperature in 5600 8000 rad in 1.25 3'. Also, shows the original population, fitted with the 2D Gaussian. This is shown in 'Figures/scatter with solution.pdf'.
4. 'Activate_Gaussian_Function.py' filters the initial population (based on radii, temperature, etc.) and calls 'Gaussian_Function_Fitting.py'. It also geneates several plots that were used by the study, but are discarded from this repository.
5. 'Figures/MCMC figures' provides the MCMC plots: corner plot, log chain and distribution of the parameters. 'SI hist.jpg' is the distribution of the linear fit parameters (slope & inclination) derived by the algorithm.
   

