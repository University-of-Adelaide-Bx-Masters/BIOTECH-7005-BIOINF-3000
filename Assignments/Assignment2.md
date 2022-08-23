# Assignment 2 [*29 marks*]

**Due before 12pm, Friday 26th August, 2022**

Your answers to all questions should be submitted to myUni as a `.zip` file containing three files:
1) a bash script for Q1 2) a bash script for Q2, and 3) the answers to the statistics questions (Q3 and Q4) in a single Rmarkdown document.
Note that the file `my_species_gff_features.txt` is not required as part of your submission for Q1, **only the script which will generate this file**!
Similarly, for Q2, only the script is required.

## Required scripts [*14 marks*]

1. Write a script to:
    + Download the gff3 file for your assigned species ([see bottom of page](#species-for-question-1)) to your current directory from Ensembl [*1 mark*]
    + Count how many of each feature type there is, sorted in numerical order [*3 marks*]
    + Export the results to a file with a name of the form `my_species_gff_features.txt` **where you use your assigned species name instead of** `my_species` [*1 mark*].
    NB: If your actual species is not included in the name, no marks will be given.
    + The script must also include code to generate one or more comment lines in the output file/table before the table with the genome-build used, (hint: grep your gff to find the genome build info as the header is very large in most cases)
    + The script must also write the code used used to generate the summary (counts) data to the output file as part of the file header. [*2 marks*]
    
2. For the file we used in the practicals (Drosophila_melanogaster.BDGP6.ncrna.fa), add to the final practical script provided so that:
    + the output contains a meaningful header [*1 mark*]
    + the output contains column names [*2 marks*]
    + the output includes: a) gene id; b) chromosome; c) start; d) stop; e) strand and f) gene_biotype [*3 marks*]
    + Appropriate comments which make the script easier to understand [*1 mark*]
    
NB: If identical comments are identified in any submissions, a mark of zero will be given for this question for all suspicious submissions.

## Statistics questions [*15 marks*]

In a single rmarkdown file answer the following questions:

3 Two groups of people have volunteered to take part in a genetic study. Group 1 (n = 126) are volunteers with no history of Type I Diabetes in their immediate family, whilst Group 2 (n = 183) have all been diagnosed with Type I Diabetes. A genotyping study was undertaken on these volunteers using 25,786 SNPs selected due to their proximity to key immune genes.
Researchers are looking to identify any SNP genotypes which may increase the risk of Type I Diabetes. In your answer, consider the reference SNP allele as `A` and the alternate SNP allele as `B`, using the genotypes `AA`, `AB` and `BB`.

a. For an individual SNP, what test would be appropriate for this comparison? [*1 mark*]  
b. Define H₀ and Hₐ for the genotype at each individual SNP. [*2 marks*]  
c. If there was no true difference in any genotypes between the two groups, how many p-values would you expect to see < 0.05? [*1 mark*]  
d. Using Bonferroni's method, what would a suitable cutoff value be to consider a SNP as being associated with an increased risk of Type I diabetes, i.e., to reject H₀ [*1 mark*]  
e. Given the following genotype table, would you reject or fail to reject H₀? Provide your working and a full explanation. [*3 marks*]

| Group | AA   | AB  | BB |
| ----- | ---- | --- | --- |
| Control | 25 | 60  | 41 |
| T1D     | 21 | 55 | 103 |


4 An experiment was repeated multiple times, in which GFP fluorescence was measured in a cell culture as a measurement of gene expression, both *before* and *after* viral transfection.
GFP was present on a plasmid as a reporter for activity at a specific promoter.
The change in fluorescence values obtained for each repeat are given [below as the vector `x`](#values-for-question-4), presented on the log2 scale for your individual subset of experiments.  

a. Define H₀ and Hₐ [*2 marks*]  
b. Calculate the sample mean and sample variance in `R` [*2 marks*]  
c. Calculate the *T*-statistic using `R`. [*1 mark*]  
d. What would the degrees of freedom be for your *t*-test? [*1 mark*]  
e. Calculate the *p*-value using `R` [*1 mark*]

Show all working & code.

## Species For Question 1

*If your student number is not listed, please contact Dave to ensure you are added to the list*

You can download your assigned species here: 'http://ftp.ensembl.org/pub/release-100/gff3/' of course you will have to add the relevant additional information to specify your species and the '.100.gff3.gz' file. 

| ID            | Species                                  | Taxonomy ID         | Common Name                               | 
|:-------------:|:----------------------------------------:|:-------------------:|:-----------------------------------------:|  
|        1717954|                          Pogona vitticeps|               103695|                     Central Bearded Dragon|  
|        1762492|                           Geospiza fortis|                48883|                        Medium Ground-Finch|          
|        1800454|                      Lepidothrix coronata|               321398|                       Blue-Crowned Manakin|          
|        1670158|                              mus musculus|                10090|                                           |          
|        1796775|                              Pan paniscus|                 9597|                           Pygmy Chimpanzee|          
|        1807015|                           Clupea harengus|                 7950|                           Atlantic Herring|          
|        1795963|                                 Mola mola|                94237|                              Ocean Sunfish|          
|        1770922|                         sus scrofa usmarc|                 9823|                                           |          
|        1770630|                         Oryzias javanicus|               123683|                          Javanese Ricefish|          
|        1645188|                   panthera tigris altalca|                 9694|                                 Amur Tiger|          
|        1802655|                        Fukomys damarensis|               885580|                            Damara Mole-Rat|          
|        1774582|                     chrysemys picta belli|                 8479|                     Western Painted Turtle|          
|        1774707|                     Fundulus heteroclitus|                 8078|                                  Mummichog|          
|        1773581|                       Macaca fascicularis|                 9541|                        Crab-Eating Macaque|          
|        1737759|                         Poecilia mexicana|                48701|                                           |          
|        1772501|              terrapene carolina triunguis|               158814|                      Three-Toed Box Turtle|          
|        1810020|                           Serinus canaria|                 9135|                              Common Canary|          
|        1749559|                       Larimichthys crocea|               215358|                       Large Yellow Croaker|          
|        1794893|                         Otolemur garnetti|                30611|                         Small-Eared Galago|          
|        1798390|                       Poecilia reticulata|                 8081|                                      Guppy|          
|        1772000|                     sus scrofa berrkshire|                 9823|                                           |          
|        1794885|                     mus musculus c57bl6nj|                10090|                                           |          
|        1707609|                           Aotus nancymaae|                37293|                          Ma's Night Monkey|          
|        1781459|                                Ovis aries|                 9940|                                      Sheep|


## Values For Question 4

*If your student number is not listed, please contact Dave to ensure you are added to the list*

The results you are analysing for Q4 are as follows.
You can simply paste these values into your RMarkdown document as the object `x` and perform all of your analysis on these values.

| ID       | Values                                                                                     |
|:---------|:-------------------------------------------------------------------------------------------|
| a1645188 | x <- c(0.4135, 2.0186, 2.3253, -0.9883, 2.0445)                                            |
| a1670158 | x <- c(3.4078, -0.362, -0.1241, 0.17, -0.2975, -1.1878, 2.5201, -0.7974, 0.6078, 1.8075)   |
| a1707609 | x <- c(0.4166, 2.1107, -0.3169, 0.5991, 0.3768, -1.5983, -0.8858)                          |
| a1717954 | x <- c(0.9675, -0.07, -0.7075, -0.316, 0.5211, 0.897, 3.0854, 0.6515, -3.9977)             |
| a1737759 | x <- c(-3.1958, 1.8503, -0.2122, -0.567, 1.1772, 1.0677)                                   |
| a1749559 | x <- c(0.2091, 1.0473, 0.1945, 1.5549, -1.2558, 0.5538, 0.9481, 0.6317, 1.0508, 0.5358)    |
| a1762492 | x <- c(1.378, 1.8755, 1.9961, -0.6889, 1.3229, 0.2232, 0.8853, -0.9496, 3.2412, 3.1686)    |
| a1770630 | x <- c(1.5092, -1.7895, -2.064, 0.8686, -0.7827, 3.7699, 2.2285, 1.0469, 3.081, 0.7944)    |
| a1770922 | x <- c(-1.3034, 0.3264, 0.105, 3.1695, 0.7723, 1.5884, -2.5724, 0.8057)                    |
| a1772000 | x <- c(0.5582, 1.348, -0.6235, 0.4487, 2.8151, 0.5023, 3.2686)                             |
| a1772501 | x <- c(-0.4398, 0.8425, 3.6124, 1.9122)                                                    |
| a1773581 | x <- c(-1.6214, -1.4619, 0.6893, -0.5449, 0.6709)                                          |
| a1774582 | x <- c(-0.4948, -0.3286, -3.8963, 0.1328, -0.113)                                          |
| a1774707 | x <- c(2.5423, 1.5786, -0.2554, 2.6932, 1.3397, 0.3383)                                    |
| a1781459 | x <- c(0.6505, -0.2976, 2.0017, 2.8812, 0.4608, 0.1071, -1.9532)                           |
| a1794885 | x <- c(-2.4143, 0.3796, 3.9181, -0.1258)                                                   |
| a1794893 | x <- c(2.1058, 1.5335, 3.7292, 1.1219, 1.9459, -1.8184, 3.6827, -0.5124, -0.4709, 2.1627)  |
| a1795963 | x <- c(0.8398, 2.1176, 0.4784, -0.7237, 0.0217)                                            |
| a1796775 | x <- c(0.1375, 1.7714, -0.1547, -2.5228, 0.7608, -2.3798, 0.1486, -0.3732)                 |
| a1798390 | x <- c(-0.7753, 0.2618, 0.1188, -0.8174, 3.6561, 2.0872, 0.5155, -2.0312, -1.3291, 0.0941) |
| a1800454 | x <- c(1.4611, 1.9262, 1.6159, 1.9525, 1.659)                                              |
| a1802655 | x <- c(-1.0547, 0.927, -1.7525, 1.0935, 1.7637, -1.2442)                                   |
| a1807015 | x <- c(0.6538, 2.5586, -0.5123, -2.505, 0.1218, 0.1247, 0.3958, -2.2144, 1.5001)           |
| a1810020 | x <- c(-0.4628, 1.0185, -0.4619, 3.7815, 1.4636, -0.148, 0.6344, -1.0034, 1.5235)          |
