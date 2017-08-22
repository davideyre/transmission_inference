# Transmission Inference

This repository contains a series of healthcare-associated tranmission tools currently under development.

1. Inference code - C++ source code can found in src
2. Simulation code - R and Python code for simulating transmission events and associated epidemiological and genetic data
3. Reporting - R scripts for checking inferrence outputs against simulations, and visualisation scripts

For further details please contact david.eyre@ndm.ox.ac.uk

---

## Running the code
 
### Input files
All times are numbered for input files starting from 1. Assumes only single infection per patient. Assumes only single sample per patient.

1. input/patientLog.csv : file with headers -   
    patient\_id [string] any patient identfier  
    t\_inf [int] infection time if know (only used for testing)  
    source [string] infection source patient identifier (only used for testing), must be "-1" for backrgound source types and starting the simulation infected  
    source_type [int] infection source type (only used for testing: BGROUND\_HOSP = 0, WARD = 1, HOSP = 2, BGROUND\_COMM = 3, START\_POS = 4, SPORE = 5)  
    t\_sample [int] sample time
    t\_recover [int] recovery time if known (only used for testing)

    where testing values are not known enter "NA" in the cell

2. input/wardLog.csv : file with headers - 
    patient\_id [string] any patient identfier  
    ward [string] any ward identfier  
    hospital [string] any hospital identfier  
    t\_admit [int] ward admission time  
    t\_discharge [int] ward discharge time  

3. input/simDistancesSNPs.txt - matrix of SNP distances, seperated by spaces, with column and row labels using patient ids, column header has no leading space

Example input files are provided in the example_input folder
