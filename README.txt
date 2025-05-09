What is GeneConv?
GeneConv is the source code used to estimate the probabilities of single/multi-SNV gene conversion events.
The probabilities were used to estimate the distribution of gene conversion tract length from HiFi PacBio datasets. 
For more details, please see Charmouh et al. 2025, Estimating gene conversion tract length and rate from PacBio HiFi data, Molecular Biology & Evolution, Volume 42, Issue 2, February 2025, https://doi.org/10.1093/molbev/msaf019.

Parameters
All parameter of the model are defined in the parameter.h file. Parameter values can be changed in parameters.cpp before producing an executable. 
The default parameters are those used in the main text.

How to compile GeneConv:
Simulations were conducted on a high performance computer (Linux). The executable was compiled in Release mode (slurm 23.02.5).
To compile the code on a Linux system, the macro LINUX in Distributions.h must be set to LINUX 0. To compile on a Windows system, LINUX in Distributions.h must be set to LINUX 0.

One way of creating an executable on a Linux system is:
Collect all cpp and h files from this repository into a folder named src. 
Add the Makefile to the directory containing src.
While in this directory, run the command "make all CFG=Release". 
This will create a folder called bin.Release which will contain your executable.

If you find the model useful, please be so kind as to cite our paper in MBE: 

Anders Poulsen Charmouh, Peter Sørud Porsborg, Lasse Thorup Hansen, Søren Besenbacher, Sofia Boeg Winge, Kristian Almstrup, Asger Hobolth, Thomas Bataillon, Mikkel Heide Schierup, Estimating Gene Conversion Tract Length and Rate From PacBio HiFi Data, Molecular Biology and Evolution, Volume 42, Issue 2, February 2025, msaf019, https://doi.org/10.1093/molbev/msaf019.
