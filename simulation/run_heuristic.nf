#!/usr/bin/env nextflow

/*
Example to run (from root of repository):
nextflow run simulation/run_heuristic.nf --simDirPath sim_data/50_scenarios --rscript reporting/heuristic_classifier.R -profile cluster -resume 
*/


// parameters 
params.simDirPath = ""
params.rscript = ""

// initial logging
log.info "\n" 
log.info "Perform inference over a collection of simulations -- version 0.1"
log.info "Simulation directory           :  ${params.simDirPath}"
log.info "Path to R script               :  ${params.rscript}"
log.info "\n"

Channel
	.fromPath("${params.simDirPath}/*", type: 'dir')
	.set {simList}

rscript = file(params.rscript)


process runH {

	input:
		file simDir from simList
	
	"""
	$rscript -d $simDir
	"""

}



