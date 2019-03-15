#!/usr/bin/env nextflow

/*
Example to run (from root of repository):
nextflow run simulation/analyse_simulations.nf --simDirPath sim_data/5_scenarios --infBin bin/transmission_test -resume -profile standard
*/


// parameters 
params.simDirPath = ""
params.infBin = ""

// initial logging
log.info "\n" 
log.info "Perform inference over a collection of simulations -- version 0.1"
log.info "Simulation directory        :  ${params.simDirPath}"
log.info "Path to binary for MCMC     :  ${params.infBin}"
log.info "\n"

Channel
	.fromPath("${params.simDirPath}/*", type: 'dir')
	.set {simList}

infBin = file(params.infBin)


process runInference {

	input:
		file simDir from simList
	
	"""
	$infBin -i 5000 -p $simDir
	"""

}

