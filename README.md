# CONTUR Cross-Section Validation
This project was carried out in the summer of 2023. It is an extension to be used in conjunction with [CONTUR](https://hepcedar.gitlab.io/contur-webpage/). Currently, it is only able to evaluate points in the [Heavy Dark Mesons](https://hepcedar.gitlab.io/contur-webpage/results/HeavyDarkMesons/index.html) model; however, it can be extended to other Beyond the Standard Model theories as detailed below.

## Operation
In order for the relevant files to be created ensure that CONTUR has been run on the .yoda file inside the point's directory and `contur-mkhtml --all --all-errbars` has also been run.

This procedure works by using the branching ratios from the .log file and the cross sections from the .out file to find all the decay paths that lead to the products contributing to the dominant subpool and finding its net cross-section and its uncertainty. Then the area under all of the histograms for a particular analysis is retrieved and the most sensitive subpool is compared against the calculated cross-section. If all is well the calculated cross-section should be larger than the histogram cross-section. Note: other decay products not implemented could also contribute to the cross-section and these should be investigated if the histogram cross-section seems higher than the calculated cross-section.

Additionally, the histogram areas are compared against the .yoda file cross-sections to check for inconsistencies.

A crude approximation is also available that uses the parent particles and fixed branching ratios. If this is to be used then it should be manually adjusted for each point.

### Adding additional detection modes
When adding new detection modes the possible combinations of products should first be added in the `DataReader` class initialisation in the possibleProductComb dictionary in the format: `'detection_mode' : ((comb1), (comb2), ...)`. If there are further conditions, e.g. on the mass of the particle that affects what is being detected, also implement them here.

In `dataReader.getBranchingRatios` determine whether standard model decays of W and Z are to be included. If **not** add the detection mode name to the tuple just below `SM_decays`. If the resultant particles strictly must decay from the same parent, implement a branch off the if-elif chain which merges the string of the resultant particles.

Note: in some models, particles may not be found by `dataReader.doInvert`, these should be implemented in the same manner as the other branches to get the antiparticle version if Dirac or `pass` if Majorana. 

In the `CrossSectionComputer` class initialiser if particles must strictly decay from the same parent (as implemented above in `dataReader.getBranchingRatios`) then the desired products must be changed to this. The final branching ratios for the last level(s) of the decay chain should be implemented below this, e.g. $Z^0 \rightarrow \text{leptons}$.

Note: the branching ratio uncertainty is an arbitrary magic number (currently set to $10^{-4}$) that should be adjusted depending on the minimum precision of the available branching ratios tabulated, generated or otherwise present.

### Background type

The histogram areas are currently using SM as background data (read from the `ANALYSIS.Summary.txt` file) and are set at the end of `DataReader.getPool`. To select DATABG, EXP, and SMBG use 0,1 and 2 respectively in the `INDEX` of `exclusions = np.asarray(exclusions).T[INDEX]`. This will affect the detection method used when computing the cross sections.



