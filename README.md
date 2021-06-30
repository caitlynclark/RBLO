# RBLO
A reliability-based layout optimization model for offshore wind turbines. Currently under development.

For technical questions regarding FLORIS, please contact Caitlyn Clark (caitlyn.clark@nrel.gov).

A couple of publications with practical information on using floris as a modeling and simulation tool for controls research are

Clark, C.E., Barter, G., DuPont, B.: Reliability-Based Layout Optimization of Offshore Wind Energy Systems, Wind Energy, 2021.

# Citation

If RBLO played a role in your research, please cite it. This software can be cited as:

RBLO. Version 1.0.0 (2019). Available at https://github.com/caitlynclark/RBLO.
For LaTeX users:

@misc{RBLO_2019,
author = {Caitlyn E. Clark},
title = {{RBLO. Version 1.0.0}},
year = {2019},
publisher = {GitHub},
journal = {GitHub repository},
url = {https://github.com/caitlynclark/RBLO}
}

# Using RBLO

To use the RBLO code, you may want to clone this repository for development, or simply copy and past to use pieces in your own code. This repository contains three directories:

1] The DrivetrainModel directory enables the calculation of planet bearing forces and L10 life of planet bearings based on FAST.Farm or FAST output ascii files. 

2] The RBLO_BearingExample directory contains an example to complete reliability-based layout optimization based on planet bearing reliablity.

3] The Surrogate directory contains the code to process FAST.Farm data generated over 36 inflow wind conditions to complete RBLO with the RBLO_BearingExample directory's code.  

There are a few equations referenced in the code, which come from these sources:

1] Guo, Yi, Jonathan Keller, and William LaCava. "Planetary gear load sharing of wind turbine drivetrains subjected to non‐torque loads." Wind Energy 18.4 (2015): 757-768.

2] Zaretsky, Erwin V. "Rolling bearing life prediction, theory, and application." (2013).

The source code repository must be cloned directly from GitHub:

git clone https://github.com/caitlynclark/RBLO

# Dependencies
RBLO has dependencies on various math, statistics, and plotting libraries in addition to other general purpose packages. For the simulation and tool modules, the dependencies are: numpy, scipy, pandas, math, and matplotlib. 

# License
Copyright 2019 NREL

Licensed under the Apache License, Version 2.0 (the "License"); you may not use this file except in compliance with the License. You may obtain a copy of the License at

http://www.apache.org/licenses/LICENSE-2.0
Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the License for the specific language governing permissions and limitations under the License.
