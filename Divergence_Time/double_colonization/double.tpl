//Number of population samples (demes)
3 samples to simulate : 
//Population effective sizes (number of genes) 
48501
NPOP1
NPOP2
//Samples sizes and samples age 
32
16
16 
//Growth rates: negative growth implies population expansion 
0 
0 
0
//Number of migration matrices : 0 implies no migration between demes 
0 
//historical event: time, source, sink, migrants, new deme size, growth rate, migr mat index 
2 historical event 
TDIV1 0 1 1 RESIZE1 0 0
TDIV2 2 1 1 RESIZE2 0 0 
//Number of independent loci [chromosome] 
1 0 
//Per chromosome: Number of contiguous linkage Block: a block is a set of contiguous loci 
1 
//per Block:data type, number of loci, per gen recomb and mut rates 
FREQ 1 0 3.67e-8 OUTEXP