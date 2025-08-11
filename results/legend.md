# Legend for results.csv

The csv file `results.csv` contains a summary of the analysis of the networks from ODEbase. 

The columns appearin in the file have the following meaning:

* `ID` is the name of the network in the BioModels database.
* `NumberReactions` is the number of reactions after irrelevant species (that do not participate in any reactions) have been removed.
* `Consistent` is true if the network is consistent (admits positive steady states for some choice of rate constants).
* `Linear` is true if all rectant complexes are monomolecular.
* `Nondegenerate` is true if the the steady state system admits a nondegenerate zero.
* `FullDimensionalInvarianceSpace` is true if the dimension of the toric invariance group equals $n-\operatorname{rank}(C)$.
* `Toric` is true if we were able to verify that the network is toric with our methods.
* `LocallyToric` is true if we have able to verify that the networ is locally toric with our methods.
* `Multistationarity` is true if we were able to verify that the network has the capacity for multistationarity with our methods.
* `ACR` is true if we were able to vertify that the network has ACR with our methods.
* `LocalACR` is true if we were able to vertify that the network has local ACR with our methods.
