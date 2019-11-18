URD is to infer the cell fate trajectory especially for time series single cell experiment.
Here I modified URD as tcURD where the diffusion map is calculated by destiny package from a sparse distance matrix, with 0 for most of element except those from nearest neighbors. Instead of calculating distance from top K nearest neighbors for each cell among all other cells, I put a time constraint that nearest neighbors could only come from samples at the same or nearby time point for the time series experiment so that we could avoid artificial connections of two cells from samples far away from each other in temporal order.


There are two input files in the script are needed:
-file.in.txt: 
    gene expression data for all samples 
    tab seprated file with gene names as header
    each column for each gene, each row for each sample.

-file.in.meta:
    meta data for each sample
    tab seprated file with each row for each sample
    should contain information to identify which timepoint sample is from