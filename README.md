# Plotting/processing codes for slowness vector and multipathing observations 

## Code descriptions

* Add_Locus_Information.py - adds the locus results to the plotting file. Takes in frequency band limits and depth. 
* Both_Vec_Multi_Tomo.py - A misleading name. It plots the slowness vector deviations + divergence, multipathing proportion and tomography map separately.
                           The user needs to give the frequency band limits, depth and region. The tomography model needs to be from the IRIS EMC database
                           [here](https://ds.iris.edu/ds/products/emc-earthmodels/) and the path to this defined in the code.
* divergence_map.py - A misleadingly named code which plots the summary of observations where subplot A shows the slowness vector observations with divergence,
                      B a tomography model from the IRIS EMC database [here](https://ds.iris.edu/ds/products/emc-earthmodels/). C shows the multipathing proportion.
                      The user needs to download the tomography model beforehand and hard-code the path to the tomography model in the code. The user gives the
                      lower and upper frequencies to plot, the region and depth. 
* extract_variance_depth_sections.py - 
* histogram_variance.py - Plot the histograms of the number of low-variance bins with depth for the given frequency band and region. 
* histogram_variance_plot_3.py - Same as above except it will plot for three specific frequency bands.
* histograms_multi_and_mag_frequencies.py - 
* locus_variance.py - 
* plot_variance_depth.py - 
* variance_neigbourhood_depth.py - 
* vec_multi_tomo_zmodel.py - Plot the summary plot specifically beneath Europe using the EU60 model. Subplot A will show the slowness vector observations
                             and calculate the divergence which is then plotted in the background. Then the EU60 model is plotted in subplot B.
                             Finally, in subplot C the multipathing proportion is plotted. The bin size can be adjusted in the code. The EU60 model is required
                             for the code to run and the model is available [here](https://academic.oup.com/gji/article/201/1/18/724841?login=false). The user gives the
                             lower and upper frequencies to plot, the region and depth. 
