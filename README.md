# ZEUS-Analysis
Code from Year 2000 Physics analysis for the ZEUS Particle Physics Detector

## Prerequisties
CERNLIB needs to be installed

## Build
 * make clean
 * make ana

## Run
* analysis.exe : Executable

_Other files_

Input
* fort.41 : List of Data ntuple files for '96 analysis
* fort.42 : List of Background ntuple files for '96 analysis
* fort.43 : List of Monte Carlo (MC) ntuple files for '96 analysis
* fort.45 : List of Data ntuple files for '97 analysis
* fort.46 : List of Background ntuple files for '97 analysis
* fort.47 : List of Monte Carlo (MC) ntuple files for '97 analysis
* st96_n : Steering cards with list of parameters. There can be many of these.
* st97_n : Steering cards for '97 analysis

Output
* fort.998 : F2 and Error
* fort.999 : Log file
* israll96/97.hbook : Output ntuple
* fort.81 : Background bin
* fort.82 : MC bin
* fort.83 : Data bin
* fort.84 : F2 ISR 
* fort.85 : Errors
* fort.500 : xtilts (96)
* fort.501 : xspreads (96)
* fort.502 : ytilts (96)
* fort.503 : yspreads (96)
* fort.504 : xtilts (97)
* fort.505 : yspreads (97)
* fort.506 : ytilts (97) 
* fort.507 : yspreads (97)

## Interpret Results

This is done using PAW. Open PAW and Enter:

\>hi/file 1 israll97.hbook

\>hi/li 

\>hi/pl 105
