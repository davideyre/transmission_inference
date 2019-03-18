#!/usr/bin/env nextflow

/*
Example to run (from root of repository):
nextflow run simulation/compare_simulations.nf --simDirPath sim_data/5_scenarios --checkMCMC reporting/checkMCMC.R -profile cluster
*/


// parameters 
params.simDirPath = ""
params.checkMCMC = ""

// initial logging
log.info "\n" 
log.info "Perform analyse inferences from a collection of simulations -- version 0.1"
log.info "Simulation directory           :  ${params.simDirPath}"
log.info "Path to comparison R script    :  ${params.checkMCMC}"
log.info "\n"

Channel
	.fromPath("${params.simDirPath}/*", type: 'dir')
	.set {simList}

checkMCMC = file(params.checkMCMC)


process compare {

	input:
		file simDir from simList
	
	"""
	Rscript $checkMCMC -d $simDir
	"""

}

