# Behav_ED
Series of MATLAB Live scripts and functions for reading and processing nonhuman primate behavioral data from the Egly, Driver, and Rafal object based spatial attention paradigm.

Data are stored in a standard csv format, in a specific file folder. This directory is to be specified within the initial file read script, "MAINlive". MAINlive also generates the initial datatable which is an aggregated table containing all trials across sessions that are stored in the data directory.

Once the table "Data" is generated, a series of functions can be implemented to perform statistical analysis and data visualization. A description of each function is listed below:

"PerfMetric.mlx" - Function which generates a series of behavioral metrics such as hit rate and segments the performance based on trial condition and spatial location of the cue.

"daily.mlx" - Function which generates daily performance charts to observe whether rhythmic changes in behavioral performance are present. 

"cumulative.mlx" - Function which generates aggregated charts to observe whether rhythmic changes in behavioral performance are present.

"Object.mlx" - Function comparing behavioral performance in the Egly task depending on whether a trial's target location and viewing bars represent the cued, uncued-same-object, or uncued-different-object location.
