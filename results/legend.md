# Legend for results.csv

The csv file `results.csv` contains a summary of the analysis of the networks from ODEbase. The columns appearing in the file have the following meaning.

## Basic information about the network
* `ID` is the name of the network in the BioModels database.
* `OriginallyMassAction` is true if the network is registered with mass-action kinetics in BioModels (regardless of this, we use mass-action kinetics for our analysis).
* `NumberReactions` is the number of reactions after irrelevant species (that do not participate in any reactions) have been removed.
* `IrrelevantSpecies` is a list of those species that does not participate in any reaction (and therefore is ignored in our analysis).
* `NumberOfSpecies` is the number of species of the network (after removing the irrelevant ones).
* `NumberOfReactions` is the number of reactions. 
* `Consistent` is true if the network is consistent (i.e., admits positive steady states for some choice of rate constants).
* `Nondegenerate` is true if the (vertically parametrized) steady state system admits a nondegenerate zero.

## Toricity
* `ToricRank` is the toric rank of the network (i.e., the rank of the toric invariance group)
* `MaximalToricRank` is true if the toric rank equals the number of species minus the rank of the network.
* `GenericallyLocallyToric` is true if the network is generically locally toric in the sense of Definition 3.4 and Theorem 5.3 of the paper. 
* `LocallyToric` is true if we were able to verify that the network is locally toric with Algorithm 6.7 in the paper.
* `Toric` is true if we were able to verify that the network is toric with Algorithm 6.7 in the paper.
* `TrivialToricinvariancePartition` is true if the matroid partition consists of a single block.
* `Quasihomogeneous` is true if the toric rank is positive and the steady state system is quasihomogeneous with respect to the elements of the toric invariance group.

## Consequences of toricity
* `Multistationarity` is true if we were able to verify that the network has the capacity for multistationarity with the methods in Section 7.
* `ACR` is true if we were able to vertify that the network has ACR with the methods in Section 7.
* `LocalACR` is true if we were able to vertify that the network has local ACR with the methods in Section 7.

## Other properties of the network
* `GenericallyBinomial` is true if we were able to prove that the steady state ideal is generically binomial [^1]. 
* `Binomial` is true if we were able to prove that the steady state ideal is binomial for all positive choices of rate constants [^1].
* `CoveredByDZT` is true if the network is covered by the deficiency zero theorem. 
* `CoveredByDOT` is true if the network is covered by the deficiency one theorem.

[^1]: We tested for binomiality of the steady state ideal in the following way:

    1. We checked if putting the coefficient matrix C of the steady state system in row reduced form gives a binomial system. If this is the case, the ideal is both generically binomial and binomial for all positive rate constants. 
    2. We plugged in a random choice of rate constants, and checking if the resulting ideal is binomial (by computing a reduced Gröbner basis with respect to the reverse lexicographic ordering). If this is the case, the ideal is generically binomial. 
    3. We checked for binomiality using the method described in Remark 3.7 in the paper.

    Missing values indicate that we were not able to obtain a Gröbner basis after 30 minutes in step 3.


