# SNPtoSNAPP

This script will generate an input nexus file for *BEAUTi* (Bouckaert et al. 2019) which can then be used to create an XML file for *SNAPP* (Bryant et al. 2012).  Input data is a SNP file in phylip format. SNPs are assumed to be unlinked. Script will phase all ambiguity codes and identify all bi-allelic SNPs. Then all SNPs with at least one individual represented per OTU will be retained. SNPs will be converted to [0,1,2] and written to nexus file. Major (0) and minor (2) alleles are coded to estimate u and v in BEAUti. This script allows subsampling. If you want to subsample the number of individuals in the analysis, you need to specify the maximum number of individuals per OTU. You will also want to set how you subsample - either by selecting individuals within OTUs at random (**random**) or by selecting individuals in OTUs in decreasing order of sequence length (**length**). The latter is useful if you want to minimize missing data.

usage:  
```python
    python SNPtoSNAPP.py infile.snps traits.file seqs_per_otu random|length
```
You need to specify a phylip data file, traits file, maximum number of individuals per OTU, and how you want the individuals to be subsampled. If you do not want to subsample your data, place a number higher than your maximum number of individuals in an OTU.

***
The traits file is tab-delimited and assigns individuals to populations:

ind1&nbsp;&nbsp;&nbsp;&nbsp;otu1  
ind2&nbsp;&nbsp;&nbsp;&nbsp;otu1  
ind3&nbsp;&nbsp;&nbsp;&nbsp;otu2  
ind4&nbsp;&nbsp;&nbsp;&nbsp;otu2  
ind5&nbsp;&nbsp;&nbsp;&nbsp;otu3  
ind6&nbsp;&nbsp;&nbsp;&nbsp;otu3  
ind7&nbsp;&nbsp;&nbsp;&nbsp;otu4  
ind8&nbsp;&nbsp;&nbsp;&nbsp;otu4  
...
***
References:

Bouckaert R, Vaughan TG, Barido-Sottani J, Duchêne S, Fourment M, Gavryushkina A, et al. (2019) BEAST 2.5: an advanced software platform for Bayesian evolutionary analysis. PLoS Computational Biology, 15, p.e1006650.  

Bryant D, Bouckaert R, Felsenstein J, Rosenberg NA, RoyChoudhury A (2012) Inferring species trees directly from biallelic genetic markers: bypassing gene trees in a full coalescent analysis. Molecular Biology and Evolution, 29, 1917-1932.
