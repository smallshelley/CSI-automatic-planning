# CSI-automatic-planning
The purpose of this script is to create robust planning automatically for CranialSpinal irradiation (CSI).
The Script is based on RayStation TPS. We have tested it on Version 9A.
Targets named PTV_brain and PTV_spinal were defined before running the script.
A database was created to predict OARs sparing previously. However, this step can be ignored by using templated OARs planning goals.
A self-adjusting method was implied in optimization in this script.
The robust optimization tool was implied. The PTV dose constraints were set as robustness with independent isocenters position uncertainty of 5 mm in superior/inferior directions.
It's recommend to standardize the OARs' naming. In our institute, the OARs naming follow the rules below: a. The first letter is uppercase and the rest are lowercase. If there were multiples words in an organ's name, underscore character ( _ ) was used to joins words. For example, Hear, Spinal_Cord. b. Laterality OAR is indicated by appending an underscore character ( _ ), followed by L or R, respectively. For example, Lung_R, Lung_L. c.  PRVs are specified as _PRV, For example, Brain_Stem_PRV, Spinal_Cord_PRV.
