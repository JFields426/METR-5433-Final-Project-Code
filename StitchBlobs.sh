#StitchBlobs joins instantaneous blocks into continuous 'blobs' that persist through multiple time steps. 
#These blobs represent unique blocking events, with each event having a unique ID.

StitchBlobs --in_list "blocktag_file.txt" \ #Input file (File created by IB_trimmed_overpole.py)
            --out_list "blockid_file.txt" \ #Output file
            --var "block_tag" \ #Variable in input file being called in (instantaneous blocks)
            --latname "latitude" \ 
            --lonname "longitude" \
            --mintime "5d" \ #Threshold for the minimum duration of blocking events
            --min_overlap_prev 50 \ #Threshold for percentage of blocking 'blob' that spatially overlaps between consecutive time steps
