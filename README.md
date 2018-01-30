The project contains files for statistical analysis of Parth models of microbial communities. 

Problem setup:

Each bacteria reqire one carbon (C) and one nitrogen (N) source for growth, one of them could be limiting. One resource may be utilized by several bacterial species, but could not be limiting for more than one. At the initial setup we have M carbon and M nitrogen sources and set of MxM possible species. Bacteria is added one by one in a random order. After thousands of addition rounds each specie was introduced multiple times, but only subset of them (2M at max) was able to survive.

analyze_C_N_models.R:
contains script for analysis of a bacterial networks constructed in a following way: vertices represent nutrients,
                       edges - survivng bacterial species, arrow direction shows limitation of a particular specie - from limiting
                       nutrient to the enriched one

functions.R:          
contains all additional functions used by analyze_C_N_models.R
