--RADIAL CHARGE DENSITY DESCRIPTOR
--

EMA4915 Senior Design

2024-2025


Language: Python


Required Dependancies / Libraries:
  numpy
  matplotlib
  scipy
  sklearn


--Directory Summary--
--
VASP Input Files - The input files for VASP that generated CHG files where RCD was calculated from
--
Descriptor - Code dedicated to generating a RCD descriptor from VASP CHG files.

  calc_rrho_edited.py - Calculates charge density, plots & writes to file.
  
  rho_gen.py - Runs calc_rrho_edited for a set of elements in directory. 
  
  descriptor.py - Generates RCD descriptor file.
  
  descriptor_win.py - Windows compat. version of descriptor.py.
--
ML - Machine Learning Code & Neccessary Input Files

  SR4.py - KRR and k-fold validation machine learning code. Plots & writes to file.
  
  DescriptorV1.txt - Pre-generated descriptor.
  
  (AtomicNumber.txt
  AtomicRadius.txt
  DescriptorEditedElectroAff.txt
  DescriptorEditedElectroNeg.txt
  EditedElectroAff.txt
  EditedElectroNeg.txt
  FirstIonEn.txt
  ZEditedElectroAff.txt
  ZEditedElectroNeg.txt) - Input files for atomic properties neccessary for code to run. 
--
