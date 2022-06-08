# BlastMovie from blast all vs all network to evolution of threshold

This script allows to take a blast file of all vs all comparison and to create an image by image file containing the evolution of the network according to the evolution of a specific threshold. For a specific threshold you will have a view od the network, but also the link between groups. Proteins of different groups that are link in the network will be shown as red.

You wil have to provide: 
- a blast file (output format 6)
- an annotation file with the protein id, the length of the protein and the column containing the name of the group the protein belong(, and the columns for the coloring of the node that could differ from the grouping column [e.g. you want the color of the node based on the subfamily but the red edges based on the name of the protein]).

## Dependencies

### Python dependencies 
- matplotlib
- networkx
- numba
- numpy
- pandas
- pygraphviz
- tqdm
  
### Other dependencies
- graphviz

### Install all dependencies using mamba/conda

```
mamba create -n blast_movie matplotlib=3.5.2 networkx=2.8.3 numba=0.55.1 numpy=1.21.6 pandas=1.4.2 pygraphviz=1.9 tqdm=4.64.0 graphviz=3.0.0
```

## Usage 

```
usage: blast_movie.py [-h] -b <file> -s {score,pident,coverage,evalue} -a <annotation> -c <column_name> [-lcc {mean,subject,query,shortest,longest}]
                      [-id {mean,subject,query,shortest,longest,HSP}] [-o <OUTPUT>] -t <num_threads>

See the effect of alignement threshold based on annotation color

options:
  -h, --help            show this help message and exit

General input dataset options:
  -b <file>, --blastfile <file>
                        Blast all vs all file
  -s {score,pident,coverage,evalue}, --selected_threshold {score,pident,coverage,evalue}
                        Choose the threshold you want to see evolve between ['score', 'pident', 'coverage', 'evalue']
  -a <annotation>, --annotation <annotation>
                        Tabulated file with the information about the sequence need to have at least, 'protein_id', 'length' and the columns to group the protein
  -g {dot,neato,fdp,sfdp,twopi,circo}, --graph_layout {dot,neato,fdp,sfdp,twopi,circo}
                        Name of the layout to plot the graph between: dot, neato, fdp, sfdp, twopi, circo
  -ce <column_name>, --cluster_egdes <column_name>
                        Name of the column in annotation file that contains the information about the clustering that will be used to know if two proteins are in the same cluster
  -c <column_name>, --column_name <column_name>
                        Name of the column in annotation file that contain the group you want to highligh in the figure. The node will be colored based on this group
  -keep, --keep_singleton
                        Option to add or not the singleton in the final image
  -lcc {mean,subject,query,shortest,longest}, --length_choice_cov {mean,subject,query,shortest,longest}
                        Length used for percentage overlap calculation between 2 sequences: 'mean'=mean of the 2 lengths (default), 'subject'=subject length, 'query'=query length,
                        'shortest'=shortest length, 'longest'=longest length
  -id {mean,subject,query,shortest,longest,HSP}, --length_choice_id {mean,subject,query,shortest,longest,HSP}
                        Length used for percentage identity calculation between 2 sequences: 'mean'=mean of the 2 lengths (default), 'subject'=subject length, 'query'=query length,
                        'shortest'=shortest length, 'longest'=longest length 'HSP'=HSP length
  -o <OUTPUT>, --output <OUTPUT>
                        Name of the output file (default: In the same folder as the blast output)
  -t <num_threads>, --threads <num_threads>
                        Number of threads to use (default:1)
```