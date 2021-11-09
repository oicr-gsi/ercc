#!/bin/bash
cd $1

#  .csv rpkmData = rpkmTable.rpkmTable
#  .png image - we are not checking this, binary content
#  .pdf reportPdf file to be checked, need to remove creation date from the header
#  .json is tested with jq

echo ".csv files:"
# calculate md5sum, no stochastic content expected
find . -name "*.csv" | xargs md5sum

echo ".pdf Files:"
# Need to exclude some date-specific information from the file
for f in *_DoseResponse.pdf;do tail -n +11 $f | md5sum;done

echo ".json files"
# jq complains about formatting, but our json does not have date-specific metadata. We may just calculate md5
find . -xtype f -name "*.json" | xargs md5sum | sort -V
