# RBLO
A reliability-based layout optimization model for offshore wind turbines. Currently under development.

For technical questions regarding FLORIS, please contact Caitlyn Clark.

A couple of publications with practical information on using floris as a modeling and simulation tool for controls research are

Clark, C.E., Barter, G., DuPont, B.: Reliability-Based Layout Optimization of Offshore Wind Turbines, In Preparation, 2020.

Clark, C.E., Barter, G., Shaler, K., Guo, Y., DuPont, B.: Surrogate Modeling for Wind Turbine Bearing Lifetime Estimation Based on FAST.Farm Simulations, In Preparation, 2020.

# Citation

If FLORIS played a role in your research, please cite it. This software can be cited as:

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

To use the RBLO code, you may want to clone this repository for development, or simply copy and past to use pieces in your own code.

The source code repository must be cloned directly from GitHub:

git clone https://github.com/caitlynclark/RBLO

# Dependencies
RBLO has dependencies on various math, statistics, and plotting libraries in addition to other general purpose packages. For the simulation and tool modules, the dependencies are: numpy, scipy, pandas, math, and matplotlib. 

# License
Copyright 2019 NREL

Licensed under the Apache License, Version 2.0 (the "License"); you may not use this file except in compliance with the License. You may obtain a copy of the License at

http://www.apache.org/licenses/LICENSE-2.0
Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the License for the specific language governing permissions and limitations under the License.
