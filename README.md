[URD](https://github.com/farrellja/URD) is to infer the cell fate trajectory especially for time series single cell experiment . However, for some circumstance, samples from far apart time points (far away in terms of biology) could be very near computationally, like our case of human primordial germline cell like cell (PGCLC) development, which will seriously bias our trajector results.

Here I modified URD as time constraint URD (tcURD) to make sure that nearest neighbors could only come from samples at the same or nearby time points for the time series experiment, so that we could avoid artificial connections of two cells from samples far away from each other in temporal order.

In brief, in URD the diffusion map is calculated by destiny package from a sparse distance matrix, with 0 for most of element except those from nearest neighbors. in tcURD Instead of calculating distance from top K nearest neighbors for each cell among all other cells, I put a time constraint : the distances are also set to 0 if two cells are far away (defined by user) in terms of the time points they came from.


There are two input files in the script are needed:
-file.in.txt: <br/>
    |_gene expression data for all samples <br/>
    |_tab seprated file with gene names as header<br/>
    |_each column for each gene, each row for each sample.<br/>

-file.in.meta:<br/>
    |_meta data for each sample<br/>
    |_tab seprated file with each row for each sample<br/>
    |_should contain information to identify which timepoint sample is from<br/>


Output files:
-seriers of RDS files record URD object, which could be used for visualize/analyses in original URD package

-pdf files to visualize lineages


check our paper "Human Primodial Germ Cells are Specified from Lineage Primed Progenitors" reently acepted on Cell Reports, with previewd version here https://papers.ssrn.com/sol3/papers.cfm?abstract_id=3443016.
