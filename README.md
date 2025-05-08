This repository contains code that computes and plots the composite difference of lagged SST anomalies between blocked and non-blocked days in the Northern Atlantic region (30N-80N, 70W-0E).

THE FOLLOWING ORDERED LIST DESCRIBES THE STEPS TO COMPUTE AND TEST SIGNIFICANCE OF THE COMPOSITE DIFFERENCE:

1. IB_trimmed_overpole.py - Calls in ERA5 reanalysis data for geopotential height and sea surface temperature and identifies points where instantaneous blocking (IB) occurs at each time step using the AGP threshold metric (see paper). This includes a trimming function to remove low-latitude blocks and altered bounds to identify IBs at high latitudes.

2. StitchBlobs.sh - Stitches IB points into 'blobs' with unique IDs that span multiple time steps. These blobs each identify a blocking event.
  - blocktag_file.txt - Directory for the input file for the function (calls in the IB_trimmed_overpole output file).
  - blockid_file.txt - Directory for the output file for the function.

3. region_block.py - Determines the central coordinate of blocking events at their first time step. Subsequently isolates blocks that form in the Northern Atlantic reigon.

4. blockinfo.py - Creates a .csv file with information for each block including duration, mininum and maximum size, etc. Useful for verifying the amount of unique blocking events. (NOTE: StitchBlobs gives each block a sequential ID based on start time. The region filtering does not reassign ID numbers, so non-consecutive block IDs are normal and expected.)

5. block_freq.py - Visualizes the frequency of blocking events over the Northern Atlantic region.

6. sst_climo.py - Calculates SST climatology.

7. sst_anom.py - Calcualtes SST anomalies based on climatology.

8. sst_anom_statprep.py - Detrends and standardizes SST anomalies to account for the change in average SST over time and to better assess statistical anomalies.

9. composite.py - Calculates composite difference between blocked and non-blocked days and plots the resulting maps for each lag. Additionally, it creates a list of days considered blocked an non-blocked for each lag (this is important for significance testing).

10. bootstrap.py - Tests significance of the composite differences using bootstrapping. Plots sample PDF histogram for a random point to visualize significance. Plots composite difference maps for each lag with stippling for stastistically significant points.

If you have any questions, please contact Jacob Fields at jacobfields@ou.edu
