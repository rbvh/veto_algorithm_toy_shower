[![arXiv](https://img.shields.io/badge/arXiv-2205.01697-<COLOR>.svg)](https://arxiv.org/pdf/1605.09246.pdf)

`veto_algorithm_toy_shower` is a reimplementation of the code associated with [Competing Sudakov Veto Algorithms](https://arxiv.org/pdf/1605.09246.pdf), for a talk at the [Taming the accuracy of event generators](https://indico.cern.ch/event/876082/) workshop. It consists of a simple toy parton shower, and an implementation of a variety of veto algorithms.

## Usage
The `Makefile` in `ToyShower` generates an executable `ShowerRunner`, which can be ran with option `-shower [SHOWER]`, where `[SHOWER]` is one of `vetomaxdumb`, `vetomax`, `maxvetodumb`, `maxveto`, `maxvetosmart`, `generateselectdumb`, `generateselect`. Explanations of the algorithms can be found in the slides in `Talk/CompetingVetoAlgorithms.pdf`. 

The code is dependent on [ROOT](https://root.cern/) for histogramming. 