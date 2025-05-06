This repository contains code that computes and plots the composite difference of lagged SST anomalies between blocked and non-blocked days in the Northern Atlantic region (30N-80N, 70W-0E).

THE FOLLOWING ORDERED LIST DESCRIBES THE STEPS TO COMPUTE AND TEST SIGNIFICANCE OF THE COMPOSITE DIFFERENCE:

1. IB_trimmed_overpole.py - Identifies points where instantaneous blocking (IB) occurs at each time step using the AGP threshold metric. This                                includes a trimming function to remove low-latitude blocks and altered bounds to identify IBs at high latitudes.

2. StitchBlobs.sh - Stitches IB points into 'blobs' with unique IDs that span multiple time steps. These blobs each represent a blocking event.
  - blocktag_file.txt - Directory for the input file for the function (calls in the IB_trimmed_overpole output file).
  - blockid_file.txt - Directory for the output file for the function.

4. region_block.py - Determines the central coordinate of blocking events at their first time step. Subsequently isolates blocks that form in                         the Northern Atlantic reigon.

5. blockinfo.py - 

6. sst_climo.py

7. sst_anom.py

8. sst_anom_statprep.py

9. composite.py

10. bootstrap.py
