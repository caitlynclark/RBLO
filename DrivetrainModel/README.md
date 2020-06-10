This directory contains all the files you need to calculate planet bearing forces and L10 life from an output ASCII FAST.Farm file. This directory does not contain any optimization tools, and is instead meant to isolate the L10 calculation used in the reliability-based design optimization code. For optimization, refer to the optimization folder. The planet bearing force code has been adapted from Eq. 9 in [1] for the 5MW NREL reference wind turbine, with parameters provided by [3]. Any discrepancies from Eq. 9 have been confirmed and validated with Yi Guo (NREL) most recently in June 2020.

You will need the drivetrain.py file, which contains all classes and functions to calculate planet bearing forces and speeds, and the corresponding L10 life of the bearings. The Calc_L10.py file calls those classes and functions and plots the forces on the bearing. The Case_A_8_dist630.0_off0.0_T2.csv file contains sample data file from FAST.Farm which will provide the basis for the life calculation. This data file is compressed for download, so it will need to be unzipped before use.

Equations and parameters referenced in the code come from these sources:

1] Guo, Yi, Jonathan Keller, and William LaCava. "Planetary gear load sharing of wind turbine drivetrains subjected to non‐torque loads." Wind Energy 18.4 (2015): 757-768.

2] Zaretsky, Erwin V. "Rolling bearing life prediction, theory, and application." (2013).

3] Nejad, Amir Rasekhi, et al. "Development of a 5 MW reference gearbox for offshore wind turbines." Wind Energy 19.6 (2016): 1089-1106.
