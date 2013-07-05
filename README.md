Ntw_fn
===============

Some code and functions to work with networks (R).

Includes:
---------
### Fork of second.extinct function on bipartite package, R

Description: 
This function is an update to the ```second.extinct``` function in the bipartite package to add a method for removing species according to its topological role in the network based on Modularity analysis. The method removes first network hubs, next module hubs, next connectors and finnally peripherial species. Species are removed randomly withinn a role.

### Hierarchihcal models

Description: Example with real data on how to run hierarchichal models that incorporate pollinator detectability and plant traits following my [PLoS One article](http://www.plosone.org/article/info%3Adoi%2F10.1371%2Fjournal.pone.0069200). It also has a simulation to apply the same idea to do food webs from gut content analysis. The basics of the simulation are aplicable tothe p-p ntw as well. This file is easier to follow than the real data example if you are interested to see how it works.

### Simulated Networks

Description: I tried to create biologically meaninful random networks, but I end up not using them. Sim_ ntw.R does it in a quite elaborated way. I also run some hierarchichal models on them. Sim_ntw2.R has simplier alternatives and try to reproduce Bluthegen et al. 2008 approach, but I didn't quite succeed (or invest much time).

