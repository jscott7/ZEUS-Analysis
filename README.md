# Introduction
Code from Year 2000 Physics analysis for the ZEUS Particle Physics Detector
[See the wiki for more detailed information](https://github.com/jscott7/ZEUS-Analysis/wiki)

# Prerequisites
CERNLIB needs to be installed. On a 32bit Unbuntu based Linux Distro it is sufficient to install from the Software Manager. 
This has been tested on Linux Mint 19.3 Cinnamon

# Components
## data_preparation
This code loads reconstructed events from the ZEUS datastore and selects a coarse population. 
Output is a set of ntuples. 
If can also be run over generated MonteCarlo datasets. 

## fl_code
This code runs over the ntuples generates by the data_preparation step. 

### Build
 * make clean
 * make analysis

### Run
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

This is done using  Physics Analysis Workstation (PAW). 

To open an ntuple, run paw from the command line in the exe folder and enter:

```
\>hi/file 1 israll97.hbook
\>hi/li 
\>hi/pl 105
```

### Kumacs
Kumac files are macros that can be loaded to run multiple complex commands in PAW

As an example, open PAW in the kumacs/statis_plots folder and run 

```
exec isr
```

This will create a plot of an DIS-ISR Feynman diagram
