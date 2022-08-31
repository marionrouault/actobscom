# Changes-of-mind in controllable and uncontrollable environments

This repository contains analysis code and data associated with the following paper:

Rouault M., Weiss A., Lee J. K., Drugowitsch J., Chambon V.* and Wyart V.* Controllability boosts neural and cognitive correlates of changes-of-mind in uncertain environments. BioRxiv (2022) https://doi.org/10.7554/eLife.75038

Script and data files are included in the repository to enable replication of the main data analyses and statistics reported in the paper.

The folder SCRIPTS contains:


• actobscom_psychometry.m, the main Matlab code that provides performance and confidence analyses, and psychometric analyses.


• fit_model_prev and fit_model_crev, two helper scripts for fitting psychometric parameters.


• a folder DATA contains anonymised behavioral data files for Experiment 1 and Experiment 2A of the paper (labelled "ACTOBS_C" and "ACTOBS_D_rule1" respectively).


• run_fit_actobs_sim.m, the main Matlab code that allows simulations and fitting of the computational model, along with fit_actobs_sim.m and getc.m. We rely on Bayesian Adaptive Direct Search: https://github.com/lacerbi/bads




License.


This code is being released with a permissive open-source license. You should feel free to use or adapt the code as long as you follow the terms of the license, which are enumerated below. If you make use of or build on the analyses, we would appreciate that you cite the paper.


Copyright (c) 2022, Marion Rouault and Valentin Wyart.


Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
