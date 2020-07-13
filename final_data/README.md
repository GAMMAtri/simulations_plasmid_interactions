The script all_simulations_steady_state_AND_time_to_extinction_REVIEW_commented.R was used to perform the simulations for the manuscript.

Simulations were run using R version 3.4.4 (2018-03-15) -- "Someone to Lean On", on a Platform: x86_64-pc-linux-gnu (64-bit) with processor: 32-core processor × 64; memory: 125,8 GiB; os: Ubuntu 18.04.4 LTS.

WARNING: the notation for the equation parameters below differs from that in the manuscript.

All the simulation output files can be found in the same directory as this script.
There are 3 types of output files:  
gama_parameters_....txt --> a dataframe of the combinations (conjugation, fitness, loss, intracellular conjugation interactions, epistasis, etc) used as input for the simulations  
gama_steady_state_....txt --> a dataframe of the steady states of each of the combinations  
gama_time_to_extinction_....txt --> a dataframe of the time each entity in combination takes to go extinct (conditions where the plasmid was maintained are obviously not included)

There are 4 files per type:  
..singles_expanded_all.. --> simulations where there's only 1 type of interaction per combination (or no interaction at all)  
..singles_fixed_competitor.. --> simulations where there's all possible interactions per combination, but the competitor plasmid is always the same.  
..solo.. --> simulations where there's only 1 plasmid (includes controls of conjugation 10-10 and fitness >= 1)  
..solo_density --> control simulations where there's only 1 plasmid and the initial density is 1e6 instead of 5e5
