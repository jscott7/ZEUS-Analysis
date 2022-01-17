# ZEUS-Analysis
Code from Year 2000 Physics analysis for the ZEUS Particle Physics Detector

## Build
 * make clean
 * make ana

## Run
* analysis.exe : Executable

_Other files_

Input
* fort.41 : List of input ntuple files for '96 analysis
* fort.45 : List of input ntuple files for '97 analysis
* fort.46 : List of Background ntuple files
* fort.47 : List of Monte Carlo (MC) ntuple files (97?)
* st96_n : Steering cards with list of parameters. There can be many of these.
* st97_n : Steering cards for '97 analysis

Output
* fort.999 : Log file
* israll96.hbook : Output ntuple
* fort.81 : Background bin
* fort.82 : MC bin
* fort.83 : Data bin
* fort.84 : F2 ISR 
* fort.85 : Errors
