#!/usr/bin/env nextflow

/*
Example to run (from root of repository):
nextflow run simulation/analyse_simulations.nf --simDirPath sim_data/5_scenarios --infBin bin/transmission_test --iter 5000 -profile standard
*/


// parameters 
params.simDirPath = ""
params.infBin = ""
params.iter = "10000"

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
iter = params.iter


process runInference {

	input:
		file simDir from simList
	
	"""
	$infBin -i $iter -p $simDir
	"""

}

