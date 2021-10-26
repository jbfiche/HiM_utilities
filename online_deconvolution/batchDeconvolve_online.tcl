# batchDeconvolve.tcl
set tcl_files "/home/fiche/huygensDeconvolutionScripts/"

set dataList [list /mnt/grey/DATA/users/JB/Deconvolution_test/Data/scan_002_DAPI_002_ROI.tif ] 

set destDir /mnt/grey/DATA/users/JB/Deconvolution_test/ 

#source "${tcl_files}merfish.tcl"
source "${tcl_files}merfish_online.tcl"

decon_segment "$tcl_files" "$destDir" "$dataList"
