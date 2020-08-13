# Assignment 2 [*29 marks*]

**Due before 12pm, Tuesday 25th August**

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
    + Include one or more comment lines before the table detailing which build of the genome was used, and the code executed to generate the summary [*2 marks*]
    
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
d. Using Bonferroni's method, what would a suitable cutoff value be to consider a SNP as being associated with an increased risk of Type I diabetes, i.e. to reject H₀ [*1 mark*]  
e. Given the following genotype table, would you accept or reject H₀? Provide your working and a full explanation. [*3 marks*]

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

*If your student number is not listed, please contact Dan to ensure you are added to the list*

| ID       | Species                       | Taxonomy ID | Common Name                    |
|:---------|:------------------------------|------------:|:-------------------------------|
| a1627307 | Apteryx rowi                  |      308060 | Okarito Brown Kiwi             |
| a1675025 | Ursus americanus              |        9643 | American Black Bear            |
| a1680671 | Calidris pygmaea              |      425635 | Spoon-Billed Sandpiper         |
| a1681781 | mus musculus c3hhej           |       10090 |                                |
| a1686550 | Felis catus                   |        9685 | Domestic Cat                   |
| a1686683 | Rhinopithecus bieti           |       61621 | Black Snub-Nosed Monkey        |
| a1704339 | mus musculus pwkphj           |       10090 |                                |
| a1717593 | Paramormyrops kingsleyae      |     1676925 |                                |
| a1720273 | Bos mutus                     |       72004 | Wild Yak                       |
| a1723812 | Ovis aries                    |        9940 | Sheep                          |
| a1724529 | Periophthalmus magnuspinnatus |      409849 |                                |
| a1725889 | Melopsittacus undulatus       |       13146 | Budgerigar                     |
| a1731521 | Meriones unguiculatus         |       10047 | Mongolian Gerbil               |
| a1740160 | Monodelphis domestica         |       13616 | Gray Short-Tailed Opossum      |
| a1747903 | Oryzias melastigma            |       30732 | Indian Medaka                  |
| a1748603 | mus musculus lpj              |       10090 |                                |
| a1749756 | Sus scrofa                    |        9823 | Pig                            |
| a1754274 | equus asinus asinus           |        9793 |                                |
| a1754994 | Salvator merianae             |       96440 | Argentine Black And White Tegu |
| a1756655 | marmota marmota marmota       |        9993 | Alpine Marmot                  |
| a1758223 | Microtus ochrogaster          |       79684 | Prairie Vole                   |
| a1765286 | Anser brachyrhynchus          |      132585 | Pink-Footed Goose              |
| a1767858 | Mus spretus                   |       10096 | Western Wild Mouse             |
| a1776797 | Mus musculus                  |       10090 | House Mouse                    |
| a1778718 | Cavia aperea                  |       37548 | Brazilian Guinea Pig           |
| a1783527 | Cercocebus atys               |        9531 | Sooty Mangabey                 |
| a1786514 | mus musculus dba2j            |       10090 |                                |
| a1786706 | canis lupus dingo             |        9612 | Dingo                          |
| a1786921 | Xiphophorus maculatus         |        8083 | Southern Platyfish             |
| a1787867 | Takifugu rubripes             |       31033 | Torafugu                       |
| a1789951 | Xiphophorus couchianus        |       32473 | Monterrey Platyfish            |
| a1791929 | Ciona intestinalis            |        7719 | Vase Tunicate                  |
| a1792295 | Manacus vitellinus            |      328815 | Golden-Collared Manakin        |
| a1792442 | Crocodylus porosus            |        8502 | Australian Saltwater Crocodile |
| a1794835 | Neovison vison                |      452646 | American Mink                  |
| a1795646 | Propithecus coquereli         |      379532 | Coquerel's Sifaka              |
| a1805535 | Capra hircus                  |        9925 | Goat                           |
| a1805954 | Procavia capensis             |        9813 | Cape Rock Hyrax                |
| a1811380 | Anolis carolinensis           |       28377 | Green Anole                    |
| a1811637 | Astatotilapia calliptera      |        8154 | Eastern Happy                  |
| a1812900 | Oreochromis niloticus         |        8128 | Nile Tilapia                   |
| a1817954 | mus musculus cbaj             |       10090 |                                |
| a1819113 | Pelodiscus sinensis           |       13735 | Chinese Soft-Shelled Turtle    |
| a1820124 | Apteryx haastii               |        8823 | Great Spotted Kiwi             |

## Values For Question 4

*If your student number is not listed, please contact Dan to ensure you are added to the list*

The results you are analysing for Q4 are as follows.
You can simply paste these values into your RMarkdown document as the object `x` and perform all of your analysis on these values.


| ID       | Values                                                                                      |
|:---------|:--------------------------------------------------------------------------------------------|
| a1627307 | x <- c(-0.1153, 2.0603, 0.9036, 0.6667, -1.5101, -1.5915, -0.6284, 1.9533)                  |
| a1675025 | x <- c(0.7065, 3.7166, -0.1085, -1.1471, 1.9624, -0.2922, 0.1988)                           |
| a1680671 | x <- c(1.7152, -1.0975, 2.2799, 2.1761)                                                     |
| a1681781 | x <- c(-0.1976, -0.7497, 1.1868, 3.0562, -0.6506, 1.3049)                                   |
| a1686550 | x <- c(3.5552, 1.2924, -0.7957, 1.0645, 1.5927, -1.0291, 1.7104, 0.4379, -2.7134)           |
| a1686683 | x <- c(3.6991, 0.1258, 1.0976, 0.1942, 2.2363, 0.0034, -1.1038, 2.7508)                     |
| a1704339 | x <- c(1.1883, 1.9448, 0.967, -0.3096)                                                      |
| a1717593 | x <- c(1.3697, 2.969, 0.0338, 3.0862, 0.8136, -0.1801, 1.1177)                              |
| a1720273 | x <- c(4.1117, 0.1749, -2.1728, 0.5281, 1.6203, 3.9644, 2.3855, 2.8199, 2.123)              |
| a1723812 | x <- c(1.3829, -1.3075, -2.081, -2.5665, 1.4863, -1.2643, 0.8084, 1.7293)                   |
| a1724529 | x <- c(-1.4031, -0.702, 3.1607, -0.1583, 0.4303, 2.8029, -1.6023)                           |
| a1725889 | x <- c(2.7939, -0.2713, 1.1828, -0.1258, 1.2275, 1.3663)                                    |
| a1731521 | x <- c(0.5213, -0.2224, 1.3412, 2.5935, -0.0679, -1.2923)                                   |
| a1740160 | x <- c(-2.9478, -0.2247, -1.7512, 2.061)                                                    |
| a1747903 | x <- c(-0.95, 1.0276, -2.2685, 2.1789, 2.2589)                                              |
| a1748603 | x <- c(0.1066, -1.3714, 4.3887, 1.429, 2.4142, 2.5031, -1.4259, -0.212, 1.356, -1.3977)     |
| a1749756 | x <- c(1.8281, 0.6641, 0.4096, 1.536, 2.8583, 1.6934, 3.0151, -0.4723, -0.1255, 1.6146)     |
| a1754274 | x <- c(-2.6509, 3.6036, 2.2295, 1.7117, 1.4966, 0.3151, 0.4794)                             |
| a1754994 | x <- c(-0.7017, 2.2433, -1.0716, 2.2364, -0.1599, 2.0544, 4.6681, -1.8806, 1.2852)          |
| a1756655 | x <- c(0.831, 2.9776, -0.1775, 1.8246, 3.0893, -0.688, 2.004, 1.6759)                       |
| a1758223 | x <- c(2.4081, -0.1234, 1.0343, 0.1866, 1.2536, 3.0342, 1.3065, 1.9026, 1.9277, 0.8871)     |
| a1765286 | x <- c(-1.4421, 1.4301, 1.9076, 1.9282, -1.4257)                                            |
| a1767858 | x <- c(-1.3356, -1.4163, 3.043, 0.7998, -0.0624, -0.4481, 2.9551, 1.7119)                   |
| a1776797 | x <- c(-0.0907, 0.9193, 1.5532, -1.0124, 0.0248, 2.5968, 2.8069, -0.2281)                   |
| a1778718 | x <- c(0.7971, -0.822, -3.8884, -1.4897, 0.3142, 3.4211, 0.3878, 2.1248)                    |
| a1783527 | x <- c(0.9781, -2.229, -0.227, 0.8415, -1.6074, 1.0363, 1.6418, 0.3249, 0.4077)             |
| a1786514 | x <- c(0.1062, -2.2914, 2.1582, -0.2737, 1.5965, 0.3698)                                    |
| a1786706 | x <- c(0.7194, -1.4352, -1.9274, -1.1021, -0.7713)                                          |
| a1786921 | x <- c(0.1195, 1.2584, -0.7798, 0.8793)                                                     |
| a1787867 | x <- c(1.7349, -1.0616, -0.1644, 1.2586, 1.2131, 0.9119, -0.0653, 1.6789, 0.846, 2.2)       |
| a1789951 | x <- c(1.334, 2.6306, 2.6329, 2.4017, -0.4444)                                              |
| a1791929 | x <- c(0.4146, -1.3487, -1.8892, -0.5764, 1.0161, 0.1756, 1.5125, -0.6207, -0.7289, 1.9736) |
| a1792295 | x <- c(0.8889, 1.1257, 1.1197, 3.7577, -0.298, 1.8663, 2.2503, 1.1981)                      |
| a1792442 | x <- c(-1.4627, -1.9357, 3.3449, 4.4048, 0.0861)                                            |
| a1794835 | x <- c(1.5051, 0.0913, 0.0845, 0.4578, 3.1312, 0.5306, -0.7336, 1.0874, 0.1603, -1.1161)    |
| a1795646 | x <- c(1.1928, 2.8593, 1.0546, -1.9531)                                                     |
| a1805535 | x <- c(0.021, 3.9584, 1.9891, 2.0794, 1.058)                                                |
| a1805954 | x <- c(0.3024, 1.437, 0.079, 0.7626, 0.2164, 2.0368, 1.0834, 2.914, 2.0474, 0.4952)         |
| a1811380 | x <- c(0.9845, 1.1896, 1.2549, -2.4693, -0.2006, -1.0952)                                   |
| a1811637 | x <- c(-0.711, -1.0063, 1.5289, 1.447, 1.7196, 0.6219, -3.3202, 0.6254, 2.0943)             |
| a1812900 | x <- c(-0.3665, -1.0665, -2.9329, 0.0777, 1.6746, 0.0382, -0.668, 1.6509, -1.0028, -2.7258) |
| a1817954 | x <- c(-3.0869, 2.2351, 1.1912, 0.9019, 2.5681)                                             |
| a1819113 | x <- c(-0.6187, 1.2712, 2.0475, 0.9763, 1.2547, -2.4799, 0.4566)                            |
| a1820124 | x <- c(1.6177, -1.0221, -0.4349, 0.2803, 0.8798, 2.7866)                                    |
