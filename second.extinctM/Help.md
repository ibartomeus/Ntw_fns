second.extinctM
===============

Fork of second.extinct function on bipartite package, R

Description: 
------------
This function is an update to the ```second.extinct``` function in the bipartite package to add a method for removing species according to its topological role in the network based on Modularity analysis. The method removes first network hubs, next module hubs, next connectors and finnally peripherial species. Species are removed randomly withinn a role.

Note: If you want to invert the order of the simulated exctinctions (i.e. peripherial first), just atribute the periferial species to the Network Hub label, and so forth. I know this is not elegant, but it works.

Usage:
------
```
second.extinctM(web, participant = "higher", method = "abun", nrep
= 10, details = FALSE, pol_ntw_hubs = NA, pol_mod_hubs = NA,
pol_conn = NA, pol_periph = NA, plant_ntw_hubs = NA, plant_mod_hubs
= NA, plant_conn = NA, plant_periph = NA)
```

Arguments:
----------

**web:** Web is a matrix representing the interactions observed between higher trophic level species (columns) and lower trophic level species (rows). Usually this will be number of pollinators on each species of plants or number of parasitoids on each species of prey.

**participant:** high (default) or low or both, depending if you want to exterminate higher trophic level species, lower trophic level species, or a randomly chosen species of both levels; partial matching.

**method:** Random deletion of a species (```random```); according to its abundance, with least abundant going extinct first (```abundance```; default) or "```degree```" to build a sequence from the best-to-least connected species. This is the most extreme case, where the most generalist species goes extinct first (see Memmott et al. 1998). New method "```modularity```" is added, and removes first network hubs, next module hubs, next connectors and finnally peripherial species. Species are removed randomly withinn a role (see Albretch et al in prep.)

**nrep:** Number of replicates of extermination sequence. Will not be considered for method abundance.

**details:** Logical; returns details, i.e. for each replicate the sequence of secondary extinctions. If set to ```FALSE``` (default), replicated runs will be averaged.

**pol_ntw_hubs, pol_mod_hubs, ...:** Vector of species names comprising each topological role (matching web colnames and rownames). If method is different from "modularity" is depreciated.


Details:
--------

Internally, each extermination is achieved by a call to ```extinctionM```, followed by a call to ```empty```, which counts the number of all-zero columns and rows.

**Note:** The length of an extinction sequence is obviously given by the number of species in the selected trophic level. When setting participant="```both```", lengths will differ for each replicate run, since it is unpredictable in which sequence species go extinct, and hence how many secondary extinctions will pre-empt further primary extinctions.

Author(s)
---------
Carsten F. Dormann wrote the initial function and Ignasi Bartomeus addedd the modularity method.

Example
-------
```
#load data example
data(Safariland)
#need to define network roles (here I did that randomly)
pol_ntw_hubs <- colnames(web)[1:3]
pol_mod_hubs <- colnames(web)[4:10]
pol_conn <- colnames(web)[11:18]
pol_periph <- colnames(web)[19:27]
plant_ntw_hubs <- rownames(web)[1:3]
plant_mod_hubs <- rownames(web)[4:6]
plant_conn <- rownames(web)[7:8]
plant_periph <- rownames(web)[9]

x <- second.extinctM(web, participant = "lower", method = "modularity", nrep = 10, details = FALSE, pol_ntw_hubs = pol_ntw_hubs, pol_mod_hubs = pol_mod_hubs, pol_conn = pol_conn, pol_periph = pol_periph, plant_ntw_hubs = plant_ntw_hubs, plant_mod_hubs = plant_mod_hubs, plant_conn = plant_conn, plant_periph = plant_periph)
robustness(x)
```