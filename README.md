# Behav_ED
Series of MATLAB Live scripts and functions for reading and processing nonhuman primate behavioral data from the Egly, Driver, and Rafal object based spatial attention paradigm.

Data are stored in a standard csv format, in a specific file folder. This directory is to be specified within the initial file read script, "MAINlive". MAINlive also generates the initial datatable which is an aggregated table containing all trials across sessions that are stored in the data directory.

Once the table "Data" is generated, a series of functions can be implemented to perform statistical analysis and data visualization. A description of each function is listed below:

daily.mlx
cumulative.mlx
Object.mlx
rhythm.mlx
rhythmTwo.mlx
