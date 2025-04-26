--RADIAL CHARGE DENSITY DESCRIPTOR
--

EMA4915 Senior Design

2024-2025


__Language:__ Python


__Required Dependancies / Libraries:__
  numpy
  matplotlib
  scipy
  sklearn

--
Directory Summary
--

--
VASP Input Files
--
*The input files for VASP that generated CHG files where RCD was calculated from*

--
Descriptor
--
*Code dedicated to generating a RCD descriptor from VASP CHG files.*

  __calc_rrho_edited.py__ - Calculates charge density, plots & writes to file.
  
  __rho_gen.py__ - Runs calc_rrho_edited for a set of elements in directory. 
  
  __descriptor.py__ - Generates RCD descriptor file.
  
  __descriptor_win.py__ - Windows compat. version of descriptor.py.

--
ML
--
*Machine Learning Code & Neccessary Input Files*

  __SR4.py__ - KRR and k-fold validation machine learning code. Plots & writes to file.
  
  __DescriptorV1.txt__ - Pre-generated descriptor.
  
  (AtomicNumber.txt
  AtomicRadius.txt
  DescriptorEditedElectroAff.txt
  DescriptorEditedElectroNeg.txt
  EditedElectroAff.txt
  EditedElectroNeg.txt
  FirstIonEn.txt
  ZEditedElectroAff.txt
  ZEditedElectroNeg.txt) - Input files for atomic properties neccessary for code to run. 
