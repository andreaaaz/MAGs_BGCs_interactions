# MAGs_BGCs_interactions

Code used to identify ecological interactions in the ocean microbiome mediated by biosynthetic gene clusters.

## Scripts
### Workflow_interactions.r
Main pipeline to identify ecological patterns of co-exclusion and co-ocurrence

To run the script you need the tables of the metadata of BGCs and MAGs previously prepared in ExploratoryAnalisis_omd.rmd. Also you can define the way yo want to group de BGCs (eg. families: gcf, clans: gcc, products) and Genomes by microbial lineages (eg. mOTUs, genus, species, family, etc). To decrease the running time and ignore groups of BGCS and MAGs that are just in few sites you can put a treshold of minimum number of sites

#### Manual:
```bash
-m --microbial_lineage 
-b --bgc_groups
-s --minimum_sites
-i --workdir (where the metadata is)
-o --outdir (where the output is saved)
```

Example of a command:

```bash
Rscript workflow_interactions.r -m mOTUs_Species_Cluster -b gcf -s 5 -i /home/user/MAGs_BGCs_interactions/ -o /home/user/exp/2025-interacions/
```
 

