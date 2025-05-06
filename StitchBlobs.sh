StitchBlobs --in_list "hgt.19502022_blocktag_files.txt" \
            --out_list "hgt.19502022_blockid_files.txt" \
            --var "Blocking" \
            --latname "latitude" \
            --lonname "longitude" \
            --mintime "5d" \
            --min_overlap_prev 50 \

#            --in_list "hgt.19912020_blocktag_files.txt" \
#           --out_list "hgt.19912020_blockid_files.txt" \