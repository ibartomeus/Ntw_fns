Ntw_fn
===============

Some code and functions to work with networks (R).

Includes:
---------
### Fork of second.extinct function on bipartite package, R

Description: 
This function is an update to the ```second.extinct``` function in the bipartite package to add a method for removing species according to its topological role in the network based on Modularity analysis. The method removes first network hubs, next module hubs, next connectors and finnally peripherial species. Species are removed randomly withinn a role.