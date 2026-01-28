# Quantum Pulse Recovery from OPA + SFG-XFROG (MATLAB)

## Overview
- Phase retrieval algorithm for separable state generalized projections algorithm (decomposing spectrogram that is a weighted sum of independent spectrograms built from orthonormal modes)
- Reproduce recovered quantum pulse figures used in Quantum FROG Main Manuscript

## Repo Contents
This repository contains MATLAB code and example data for recovering ultrafast quantum optical pulses from measured XFROG data. 

### Folders: 
- MeasuredRawData: contains raw measured spectrogram data used to generate plots in Quantum FROG paper. 
- RecoveredModesData: contains recovered amplified vacuum and amplified squeezed modes information

### Codes: 
- svdFROGMM: Runs multimode iterative phase retrieval algorithm given spectrogram input
- svdexFROGMM: Called by svdFROGMM - applies nonlinearity constraint at each iteration
- makeXFROGMM: Called by svdFROGMM, ExpQuantumPulseRecovery, TheoryExpMacroscopicPulseRecovery

### Demos: 
- ExpQuantumPulseRecovery: Loads recovered ampified vacuum and amplified squeezed vacuum information (see Quantum FROG main manuscript) and recovers quantum pulse
- TheoryExpMacroscopicPulseRecovery: Loads theoretical or experimental spectrogram (user set) and runs phase retrieval algorithm

## System requirements
### Software requirements: 
Tested on MATLAB R2022B

### Operating system requirements: 
Tested on Windows 11

## Installation guide
- Download repository. 
- Open MATLAB and add repository to path. 
- Typical setup time: few minutes 

## Running instructions: 
### Demo
#### Run ExpQuantumPulseRecovery to reproduce quantum pulse measurement data displayed main Quantum FROG manuscript.
- Expected output: MATLAB figures displaying results. 
- Expected run time: Less than one minute

#### Run TheoryExpMacroscopicPulseRecovery to see examples of phase retrieval. 
- Expected output: MATLAB figures displaying results.
- Expected run time: Minutes to hours depending on number of iterations, size of spectrogram (can be defined by user)

#### Run on user defined spectrogram 
- Call svdFROGMM with your spectrogram to run phase retrieval on your data. Use TheoryExpMacroscopicPulseRecovery as an example on how to call svdFROGMM.

## Software licensing notice
The software is provided at no cost for academic and non-commercial use.

Copyright © 2026, California Institute of Technology.
All rights reserved.
Redistribution and use in source and binary forms for academic and other non-commercial purposes, with or without modification, are permitted provided that the following conditions are met:
•	Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
•	Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
•	Neither the name of the California Institute of Technology (Caltech) nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 
Copyright © 2026, California Institute of Technology, based on research funded under U.S. Government Grants FA9550-23-1-0755, W911NF-23-1-0048, and D23AP00158-00.
All rights reserved.

