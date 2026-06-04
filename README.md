# MAGs_BGCs_interactions

Code used to identify ecological interactions mediated by biosynthetic gene clusters (BGCs) in the ocean microbiome.

## Scripts
### Main pipeline
The script ```pipeline.nf``` has the main pipeline to identify cases of co-occurrence and bulid networks with these cases.

#### Identify cases
There are 2 ways to identify interactions:
- ```interactions_MAG-MAG.r```: Co-occurrence of a microbial lineage with other microbial lineage 
- ```interactions_MAG-BGC.r```: Co-occurrence of a microbial lineage with a cluster of BGCs

For this pipeline we paralelize with nextflow both scripts. Also the pipeline run the scripts for all the different temperatures (***global***, ***high***, ***mid*** and ***low***). Each result produce a network.

#### Networks
We build networks with the results of co-occurrence cases. Depending if the cases are MAG-MAG or MAG-BGC, they go to diferent prosesing scirpts:
- ```networks_mm.r```: Create tables of edges and nodes of the co-occurrence network MAG-MAG.
- ```networks_mb.r```: Create tables of edges and nodes of the co-occurrence network MAG-BGC. Also create networks for MAG-BGC-MAG and MAG-MAG.


#### Workflow so far...

Temperatures:     
global ----|  |--> interactions_MAG-MAG.r --> networks_mm.r --|
high   ----|->|                                               |---> info_networks.r
mid    ----|  |--> interactions_MAG-BGC.r --> networks_mb.r --| 
log    ----|


#### Manual

```bash
--indir (input directory) 
--outdir (output directory)
--microbial_lineage (mOTUs, genus, familiy)
--bgc_groups (GCFs, GCCs)
--method (binomial, mutual)
```

Example of a command:

```bash
$ nextflow run MAGs_BGCs_interactions/pipeline.nf \
    --indir /home/user/data/metadata/ \
    --outdir /home/user/exp/interactions/ \
    --microbial_lineage mOTUs_Species_Cluster \
    --bgc_groups gcc \
    --method binomial
```

#### Output saving

As we have 2 way of ifentify interaction and 4 temperatures, we have 8 reults. These results are saved in diferent directories: ```exp/interactions/${microbial_lineage}_${bgc_groups}/${temperature}/``` or ```exp/interactions/${microbial_lineage}_${microbial_lineage}/${temperature}/```. The output of the command example will be saved in this order:


interactions/ -
              |-mOTUs_Species_Cluster_gcc/ -
              |                            |- global/ -
              |                            |          |- all_cases.csv
              |                            |          |- ex_cases.csv
              |                            |          |- oc_cases.csv
              |                            |
              |                            |- high/
              |                            |- mid/
              |                            |- low/ 
              |
              |-mOTUs_Species_Cluster_mOTUs_Species_Cluster/ 
                                                              
                                     
## Acknowledgments

The work in this repository was supported by DGAPA-PAPIIT grant IA200824
