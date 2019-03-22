#!/usr/bin/env nextflow

/*
Example to run (from root of repository):
nextflow run simulation/run_simulations.nf --simDirPath sim_data/50_scenarios --infBin bin/transmission_test --iter 50000 --checkMCMC reporting/checkMCMC.R -profile cluster -resume
*/


// parameters 
params.simDirPath = ""
params.infBin = ""
params.checkMCMC = ""
params.iter = "10000"

// initial logging
log.info "\n" 
log.info "Perform inference over a collection of simulations -- version 0.1"
log.info "Simulation directory           :  ${params.simDirPath}"
log.info "Path to binary for MCMC        :  ${params.infBin}"
log.info "MCMC iterations                :  ${params.iter}"
log.info "Path to comparison R script    :  ${params.checkMCMC}"
log.info "\n"

Channel
	.fromPath("${params.simDirPath}/*", type: 'dir')
	.set {simList}

infBin = file(params.infBin)
iter = params.iter
checkMCMC = file(params.checkMCMC)

process runInference {

	input:
		file simDir from simList
	output:
		file simDir into simRun
	
	"""
	$infBin -i $iter -p $simDir
	"""

}

process compare {

	input:
		file simDir from simRun
	
	"""
	Rscript $checkMCMC -d $simDir
	"""

}


