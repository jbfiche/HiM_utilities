# This Huygens script deconvolves and segments images for MERFISH analysis.
# source "/home/marcnol/Repositories/huygensDeconvolutionScripts/merfish.tcl"


proc decon_segment {tcl_files destDir fileList} {

# running examples	
# decon_segment /mnt/PALM_dataserv/DATA/Andres/MERFISH/20170623_MERFISH_2_Readouts/005_merFISH_merFISH_slide4_RT7_8/ 
# decon_segment /mnt/PALM_dataserv/DATA/Andres/MERFISH/20170623_MERFISH_2_Readouts/004_merFISH_merFISH_slide4_DAPI/ 
# decon_segment /home/marcnol/Dropbox/projects/methodological/SVI/Images/raw
# decon_segment /mnt/PALM_dataserv/DATA/JB/2017/20171025_merFISH/DAPI/
# decon_segment /mnt/PALM_dataserv/DATA/JB/2017/20171025_merFISH/merFISH/

# decon_segment /mnt/grey/DATA/rawData_2020/Experiment_5/DAPI/005_merFISH_RAMM_DAPI

# Stop questions and process the images in my folder.
set oldVerbosity [huOpt verb -mode silent]

# source Tcl functions
# set tcl_files "/home/marcnol/Repositories/huygensDeconvolutionScripts/"
source "${tcl_files}merfish_parameters.tcl"
source "${tcl_files}convert_to_Tiff.tcl"
source "${tcl_files}log_file.tcl"
source "${tcl_files}parse_analyse_output.tcl"
source "${tcl_files}do_TIFF_fileList_from_folder.tcl"
source "${tcl_files}Utilities.tcl"

# sets parameters
set TIME_start [clock clicks -milliseconds]
set log_filename "${destDir}output_${TIME_start}.log"
set writeMode "a"
set thresholds_filename "${destDir}thresholds_${TIME_start}.log"
set singleDir "false"
set HiMexperimentDir "false"

# Retrieve the list of files stored in the source folder.
if {$singleDir == "true"} {
	huOpt report "> Reading images from a single directory"
	set fileList [glob "${rootDir}*.tif"]
} elseif {$HiMexperimentDir == "true"} {
	huOpt report "> Reading images from a HiM experiment directory"
	set fileList [do_TIFF_fileList_from_folder "$rootDir" "$log_filename" \
	"$destDir" "$tcl_files"]
	set fileAttributes0 [retrieve_file_attributes $fileList]
	array set fileAttributes $fileAttributes0
} else {
	huOpt report "> Reading images from a list of predefined names"
}


# Loop over the images located in the source folder.
set Nimages 0	
set totalImages [llength $fileList]
foreach imgInDisk $fileList {

	# Starting time of iteration
	set TIME_start_iteration [clock clicks -milliseconds]
	
	# Open this image.
	set rawImg0 [img open $imgInDisk -foreignTo float -cmode scale]

	huOpt report "> File [expr $Nimages + 1] out of $totalImages loaded"	
	log_file $log_filename "---------------------------------------------- \
		\n > File [expr $Nimages + 1] out of $totalImages loaded" $writeMode

	# Sets parameters
#	$rawImg0 setp -dx $dx -dy $dy -dz $dz -na $NA -ri $ri -ril $ril
#	for {set ch 0} {$ch < $Nchannels} {incr ch} {
#		$rawImg0 setp -chan $ch -ex $wavelength(ch$ch,ex) -chan $ch -em $wavelength(ch$ch,em)
#	}
	
	# converts image from 2-channel TIFF
		
	set rawImg [$rawImg0 repl "${rawImg0}_converted"] 
	convert2Tiff "$rawImg0" "$rawImg" "$Nchannels"
	huOpt report "> TIFF stack converted..."	
	log_file $log_filename "> TIFF stack converted..." $writeMode    

#############################################
	# orders channels so that DAPI/RT is ch0 and beads is ch1
#	set rawImg [img create rawImg -dim [$rawImg1 getdims] -type [$rawImg1 getconfig -mode type] ]
#	set rawImg [$rawImg1 repl "${rawImg0}_conv"]

#	set image [$rawImg1 split -mode all]
#	set im0 [lindex $image 0]
#	set im1 [lindex $image 1]
	
#	if {$channel_order($experiment_type) == "0"} {
#		$im0 join $im1 -> $rawImg	
#	} elseif {$channel_order($experiment_type) == "1"}	{
#		$im1 join $im0 -> $rawImg	
#	}

#############################################

	# Sets parameters
	$rawImg setp -dx $dx -dy $dy -dz $dz -na $NA -ri $ri -ril $ril
	for {set ch 0} {$ch < $Nchannels} {incr ch} {
		$rawImg setp -chan $ch -ex $wavelength(ch$ch,ex) -chan $ch -em $wavelength(ch$ch,em)
	}
	
	# Get the microscopic parameters of this raw image.
	array set rawImgParams [$rawImg setp -tclReturn]

	# Create an empty image for the deconvolved version of the image.
	set destImg [$rawImg repl "deconvolved"]

	# Deconvolve the image with default settings. 
		# Add arguments to change the SNR, nr. of iterations, etc. 
	$rawImg cmle psf -> $destImg -bgMode wf -bgRadius 0.7 -blMode off 

	# Save the deconvolved image.
#	$rawImg save "${destDir}${rawImg}_converted" -type $destType -tiffMultiDir
	$destImg save "${destDir}${rawImg}_decon" -type $destType -tiffMultiDir
	huOpt report "> Deconvolved Image saved to: ${destDir}${rawImg}_decon"	
	log_file $log_filename "> Deconvolved Image saved \
		to: ${destDir}${rawImg}_decon" $writeMode
		

	if {$thresholding_method != "none"} {	
		
	# Defines parameters to be used for this specific file
	switch $thresholding_method {
	"otsu" {
		# sets global parameters for everyone
		set garbage_value0(ch0) $garbage_value(ch0)
		set garbage_value0(ch1) $garbage_value(ch1)
		set garbage_value0(ch2) $garbage_value(ch2)	
	}
	"search" {
		# sets global parameters for everyone
		set threshold_value0(ch0) $threshold_value(ch0)
		set garbage_value0(ch0) $garbage_value(ch0)
		set threshold_value0(ch1) $threshold_value(ch1)
		set garbage_value0(ch1) $garbage_value(ch1)
		set threshold_value0(ch2) $threshold_value(ch2)
		set garbage_value0(ch2) $garbage_value(ch2)
	}
	"fixed" {
		if {$singleDir == "false"} { 
				if {$globalSegmPars == "true"} { 
					# sets global parameters for everyone
					set threshold_value0(ch0) $threshold_value(ch0)
					set garbage_value0(ch0) $garbage_value(ch0)
					set threshold_value0(ch1) $threshold_value(ch1)
					set garbage_value0(ch1) $garbage_value(ch1)
					set threshold_value0(ch2) $threshold_value(ch2)
					set garbage_value0(ch2) $garbage_value(ch2)
				} else {
					# sets specific parameters for RTs/DAPI/beads
					set threshold_value0(ch0) $threshold_value(RT$fileAttributes($Nimages.RT))
					set garbage_value0(ch0) $garbage_value(RT$fileAttributes($Nimages.RT))
					set threshold_value0(ch1) $threshold_value(beads)
					set garbage_value0(ch1) $garbage_value(beads)
					huOpt report "RT=$fileAttributes($Nimages.RT); \
						roi = $fileAttributes($Nimages.roi) \
						Thresh=$threshold_value0(ch0); garb=$garbage_value0(ch0)"
					log_file $log_filename "> RT=$fileAttributes($Nimages.RT); \
						roi = $fileAttributes($Nimages.roi); \
						Thresh=$threshold_value0(ch0); garb=$garbage_value0(ch0)" $writeMode
				}
			} else {
				# sets global parameters for everyone in singleDir mode
				set threshold_value0(ch0) $threshold_value(ch0)
				set garbage_value0(ch0) $garbage_value(ch0)
				set threshold_value0(ch1) $threshold_value(ch1)
				set garbage_value0(ch1) $garbage_value(ch1)
				set threshold_value0(ch2) $threshold_value(ch2)
				set garbage_value0(ch2) $garbage_value(ch2)
			}
	}
	}

	# allocates memory for processed images
	set thresholded_image [$destImg repl "thresholded"] 
	#set Nchannels [llength $thresholds]

	# defines thresholds to use
	switch $thresholding_method {
		"fixed" {
			# conducts segmentation with thresholds fixed in merfish_parameters.
			log_file $log_filename "? Entering pre-defined thresholding routine" $writeMode
	
			for {set ch 0} {$ch < $Nchannels} {incr ch} {
				append best_thresholds "$threshold_value0(ch$ch) "
			}			
			log_file $log_filename "> Fixed thresholds used: $best_thresholds" $writeMode

		}
		"search" {
			# searches for automatic, optimal threshold
			set threshold_split [split $threshold_value0(ch0) ","]
			set N_iterations [lindex $threshold_split 2]
			set threshold_spacing [expr [ expr [lindex $threshold_split 1] - [lindex $threshold_split 0] ] / $N_iterations ]
			log_file $log_filename "> Finding optimal threshold, $N_iterations iterations" $writeMode
			set dims [$destImg getdims]
			set type [$destImg getconfig -mode type]  
	
			if {$Nimages > 0} {
	#			log_file $log_filename "> Erased old list: $Nobj" $writeMode
				unset threshold_list
				unset N_obj
				unset Nobj
				unset best_thresholds
			}
	
			#		huPrint info "? Segments with different thresholds"
			for {set i_thresh [lindex $threshold_split 0]} {$i_thresh < [lindex $threshold_split 1]} \
				{incr i_thresh $threshold_spacing} {
	
				set thresholds ""
				set garbages ""
	
				for {set ch 0} {$ch < $Nchannels} {incr ch} {
					append thresholds "$i_thresh "
					append garbages "$garbage_value0(ch$ch) "
				}
	
				append Nobj [thresholds_labels "$destImg" "$thresholds" "$garbages" "$dims" "$border_fill" "$type"]
	
#				huPrint info "Threshold: $i_thresh"	
	
			}
	#		log_file $log_filename "> Threshold list: $Nobj" $writeMode	
	#		huPrint info "? Decodes Nobj to define N_obj(ch)"
			for {set ch 0} {$ch < $Nchannels} {incr ch} {
			
				for {set i $ch} {$i < [llength $Nobj]} {incr i [expr $Nchannels] } {
					set iNobj [lindex $Nobj $i]
					array set Nobj_array $iNobj
					append threshold_list($ch) "$Nobj_array(thresh) "
					append N_obj($ch) "$Nobj_array(Nobj) "
				}
			}
			
			for {set ch 0} {$ch < $Nchannels} {incr ch} {
				append best_thresholds "[lindex $threshold_list($ch) [find_best_threshold $N_obj($ch) $minDifference]] "
				log_file $thresholds_filename "$Nimages $ch $N_obj($ch)" $writeMode
	
			}
			log_file $log_filename "> Optimal thresholds found: $best_thresholds" $writeMode
		} 
		"otsu" {
			log_file $log_filename "? Entering Otsu routing to set thresholds" $writeMode
			set best_thresholds ""
			for {set ch 0} {$ch < $Nchannels} {incr ch} {
				set otsu_threshold [otsu "$destImg" "$ch"]
				append best_thresholds "$otsu_threshold "
			}		
			log_file $log_filename "> Optimal thresholds (Otsu) used: $best_thresholds" $writeMode

		}
	}
	
	$destImg thres -> $thresholded_image "$best_thresholds" -mode all
		
	for {set ch 0} {$ch < $Nchannels} {incr ch} {
		# gets the channel from the thresholded image
		set thresholded_image_ch [$thresholded_image split -mode one -chan $ch] 
		$thresholded_image_ch borderFill $border_fill 
		set labeled_image_ch [$thresholded_image_ch repl "labeled_ch"] 
		
		# labels thresholded channel image
		$thresholded_image_ch label -> $labeled_image_ch
	
		# gets output of analysis of clusters detected
		set segm_output_ch [ $thresholded_image_ch analyze $labeled_image_ch \
			 -garb $garbage_value0(ch$ch) -style tcl ] 
#		log_file $log_filename "> N_obj for ch${ch}: [lindex [split [$labeled_image_ch range] " "] 1]" $writeMode
		log_file $log_filename "> N_obj for ch${ch}: [llength $segm_output_ch]" $writeMode
			 
		# Save the labeled image for further processing in matlab
		$thresholded_image_ch save "${destDir}${rawImg}_thresholded_ch$ch" \
			-type tiff16 -tiffMultiDir -cmode clip		
		$labeled_image_ch save "${destDir}${rawImg}_labeled_ch$ch" \
			-type tiff16 -tiffMultiDir -cmode clip		

		# saves output data from label
		set filename "${destDir}${rawImg}_label_output_ch${ch}.dat"
		parse_analyse_output "$segm_output_ch" "$filename"	

		$labeled_image_ch del
		$thresholded_image_ch del
	}
	}
	
	#clears up memory
	$destImg del
	$rawImg0 del
	$rawImg del
#	$im0 del
#	$im1 del
#	$rawImg1 del
	
	huOpt report "> Segmented Images saved"	
	log_file $log_filename "> Segmented Images and objects saved " $writeMode

	set TIME_taken [expr [clock clicks -milliseconds] - $TIME_start_iteration]
	set TIME_human [hh:mm:ss [expr $TIME_taken / 1000 ] ]
	huOpt report "> time taken for this iteration: $TIME_human" 
	log_file $log_filename "> Cumulative time: $TIME_human" $writeMode

	# Continue to the next iteration of the loop: another raw image.
	incr Nimages
}   

# set total processing time
set TIME_taken [expr [clock clicks -milliseconds] - $TIME_start]
set TIME_human [hh:mm:ss [expr $TIME_taken / 1000 ] ]

huOpt report ">number of images: $Nimages" 
huOpt report ">time taken $TIME_human" 
	
# Reset verbosity level to old value
huOpt verb -mode $oldVerbosity

huOpt report "done."

}



# Run the batch deconvolution.
#set merfishDir "/mnt/PALM_dataserv/DATA/Andres/MERFISH/"
#set rootDir "${merfishDir}20170623_MERFISH_2_Readouts/005_merFISH_merFISH_slide4_RT7_8/"
# set destDir "/mnt/grey/DATA/rawData_2020/Experiment_5/deconvolved_DAPI_test/"
#decon_segment $rootDir

#decon_segment $sourceDir $destDir $tcl_files $srcType $destType $rootDir


