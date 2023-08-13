# Assignment 2 [*29 marks*]

**Due before 12pm, Friday 25th August, 2023**

Your answers to all questions should be submitted to myUni as a single `.zip` file containing three files:
1) a bash script for Q1 2) a bash script for Q2, and 3) a single Rmarkdown document containing the answers to the statistics questions (Q3 and Q4).
Note that the file `my_species_gff_features.txt` is not required as part of your submission for Q1, **only the script which will generate this file**!
Similarly, for Q2, only the script is required.

## Required scripts [*14 marks*]

**Q1.** Write a script to:
    + Download the gff3 file for your assigned species ([see bottom of page](#species-for-question-1)) to your current directory from Ensembl [*1 mark*]
    + Count how many of each feature type there is, sorted in numerical order [*3 marks*]
    + Export the results to a file with a name of the form `my_species_gff_features.txt` **where you use your assigned species name instead of** `my_species` [*1 mark*].
    NB: If your actual species is not included in the name, no marks will be given.
    + The script must also include code to generate one or more comment lines in the output file/table before the table with the genome-build used, (hint: grep your gff to find the genome build info as the header is very large in most cases)
    + The script must also write the code used used to generate the summary (counts) data to the output file as part of the file header. [*2 marks*]
    
**Q2.** For the file we used in the practicals (Drosophila_melanogaster.BDGP6.ncrna.fa), add to the final practical script provided so that:
    + the output contains a meaningful header [*1 mark*]
    + the output contains column names [*2 marks*]
    + the output includes: a) gene id; b) chromosome; c) start; d) stop; e) strand and f) gene_biotype [*3 marks*]
    + Appropriate comments which make the script easier to understand [*1 mark*]
    
NB: If identical comments are identified in any submissions, a mark of zero will be given for this question for all suspicious submissions.

## Statistics questions [*15 marks*]

In a single rmarkdown file answer the following questions:

**Q3.** Two groups of people have volunteered to take part in a genetic study. Group 1 (n = 126) are volunteers with no history of Type I Diabetes in their immediate family, whilst Group 2 (n = 183) have all been diagnosed with Type I Diabetes. A genotyping study was undertaken on these volunteers using 25,786 SNPs selected due to their proximity to key immune genes.
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


**Q4.** An experiment was repeated multiple times, in which GFP fluorescence was measured in a cell culture as a measurement of gene expression, both *before* and *after* viral transfection.
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

| ID       | Species                      | Taxonomy ID | Common Name                    |
|:---------|:-----------------------------|------------:|:-------------------------------|
| a1670158 | Hippocampus comes            |      109280 | Tiger Tail Seahorse            |
| a1685584 | sus scrofa bamei             |        9823 |                                |
| a1711004 | Ictidomys tridecemlineatus   |       43179 | Thirteen-Lined Ground Squirrel |
| a1731685 | Ictalurus punctatus          |        7998 | Channel Catfish                |
| a1745654 | Gambusia affinis             |       33528 | Western Mosquitofish           |
| a1749310 | Anabas testudineus           |       64144 | Climbing Perch                 |
| a1771951 | Zonotrichia albicollis       |       44394 | White-Throated Sparrow         |
| a1779334 | panthera tigris altaica      |        9694 | Amur Tiger                     |
| a1781749 | Fundulus heteroclitus        |        8078 | Mummichog                      |
| a1791787 | Ovis aries                   |        9940 | Sheep                          |
| a1793463 | Pelodiscus sinensis          |       13735 | Chinese Soft-Shelled Turtle    |
| a1809569 | Geospiza fortis              |       48883 | Medium Ground-Finch            |
| a1813640 | Poecilia formosa             |       48698 | Amazon Molly                   |
| a1818796 | bison bison bison            |        9901 |                                |
| a1820572 | sus scrofa pietrain          |        9823 |                                |
| a1826048 | Castor canadensis            |       51338 | American Beaver                |
| a1826245 | astyanax mexicanus pachon    |        7994 |                                |
| a1827429 | Strigops habroptila          |     2489341 | Kakapo                         |
| a1827586 | Xiphophorus maculatus        |        8083 | Southern Platyfish             |
| a1828691 | Lepidothrix coronata         |      321398 | Blue-Crowned Manakin           |
| a1830214 | canis lupus dingo            |        9612 | Dingo                          |
| a1838426 | mustela putorius furo        |        9668 | Domestic Ferret                |
| a1866470 | Aotus nancymaae              |       37293 | Ma's Night Monkey              |
| a1871383 | Fukomys damarensis           |      885580 | Damara Mole-Rat                |
| a1872578 | Poecilia reticulata          |        8081 | Guppy                          |
| a1873748 | mus musculus akrj            |       10090 |                                |
| a1878299 | terrapene carolina triunguis |      158814 | Three-Toed Box Turtle          |
| a1879743 | Poecilia mexicana            |       48701 |                                |
| a1880903 | Pogona vitticeps             |      103695 | Central Bearded Dragon         |
| a1881885 | Serinus canaria              |        9135 | Common Canary                  |
| a1882601 | chrysemys picta bellii       |        8479 | Western Painted Turtle         |
| a1893135 | mus musculus c57bl6nj        |       10090 |                                |
| a1893784 | Larimichthys crocea          |      215358 | Large Yellow Croaker           |
| a1894438 | Macaca fascicularis          |        9541 | Crab-Eating Macaque            |
| a1894742 | Mola mola                    |       94237 | Ocean Sunfish                  |
| a1895152 | sus scrofa usmarc            |        9823 |                                |
| a1895218 | Pan paniscus                 |        9597 | Pygmy Chimpanzee               |
| a1896450 | Otolemur garnettii           |       30611 | Small-Eared Galago             |
| a1896652 | sus scrofa berkshire         |        9823 |                                |
| a1899420 | Clupea harengus              |        7950 | Atlantic Herring               |
| a1901571 | Oryzias javanicus            |      123683 | Javanese Ricefish              |

## Values For Question 4

*If your student number is not listed, please contact Dave to ensure you are added to the list*

The results you are analysing for Q4 are as follows.
You can simply paste these values into your RMarkdown document as the object `x` and perform all of your analysis on these values.

| ID       | Values                                                                                   |
|:---------|:-----------------------------------------------------------------------------------------|
| a1670158 | x <- c(3.4078, -0.362, -0.1241, 0.17, -0.2975, -1.1878, 2.5201, -0.7974, 0.6078, 1.8075) |
| a1685584 | x <- c(1.783, 0.7986, 0.6015, 1.042, 1.3967, 3.4582, 3.87, 2.4096, 0.8027)               |
| a1711004 | x <- c(0.3384, -1.239, -0.4653, -1.2589, 1.9388, -1.0695, -1.0622, 2.3031)               |
| a1731685 | x <- c(2.2962, 2.3225, -0.6423, -0.2588)                                                 |
| a1745654 | x <- c(2.5999, 2.4134, 0.7226, -0.5628)                                                  |
| a1749310 | x <- c(0.9669, 0.9458, -1.0238, -0.537, 0.0973, 0.8976, -1.2348, 0.3315)                 |
| a1771951 | x <- c(-3.6227, -0.5281, 0.8609, 1.5519, -1.9202, 0.111)                                 |
| a1779334 | x <- c(3.0414, 0.3347, 0.2598, 0.0669, -2.273, 2.3805, -0.2023, 0.0795)                  |
| a1781749 | x <- c(-4.8139, -0.693, -0.4226, 2.3468, -0.6309, 0.1241)                                |
| a1791787 | x <- c(0.4395, 1.3864, -0.0849, 2.2884, 0.6172, 0.9241, 0.9488)                          |
| a1793463 | x <- c(-1.5012, -0.3761, 1.6316, 1.5747, 0.5813, -2.4942, 3.8002, 3.0263, 1.0505)        |
| a1809569 | x <- c(-0.4614, 0.8025, 1.1472, 0.1815, -1.51, 1.2096)                                   |
| a1813640 | x <- c(3.5413, 1.1081, -1.2902, -0.492, 1.2689, 0.2239, 2.7458, 0.5778, 4.1346, -0.1029) |
| a1818796 | x <- c(1.4362, 0.374, 0.8476, 0.5466, 0.7706, -2.423)                                    |
| a1820572 | x <- c(0.379, 2.9524, 3.2554, -0.277)                                                    |
| a1826048 | x <- c(0.3902, 1.0562, 1.8898, -1.5771)                                                  |
| a1826245 | x <- c(0.6447, 0.884, 2.1958, 0.285)                                                     |
| a1827429 | x <- c(1.3417, 0.8312, 2.3665, 2.2907)                                                   |
| a1827586 | x <- c(-0.0582, 2.0452, 1.2181, -1.3813, 2.3816, 0.6085, 0.3904, 0.1846, 1.9378)         |
| a1828691 | x <- c(2.0259, -0.57, 0.4786, -0.146, 1.4909, 1.3815, -0.9265, 0.9515)                   |
| a1830214 | x <- c(1.2928, 0.0834, 2.385, 2.4948, 1.792)                                             |
| a1838426 | x <- c(0.3964, 1.8519, 0.7053, -0.4396, 0.7592, -3.41, 2.0735)                           |
| a1866470 | x <- c(-0.4671, 1.2932, -0.1467, 1.3883, 0.1519, -1.8978, 0.1778)                        |
| a1871383 | x <- c(0.165, -0.9126, 2.5109, -0.1687, 1.1595, -3.136, 0.7482)                          |
| a1872578 | x <- c(1.916, 0.3201, 0.7413, -0.521)                                                    |
| a1873748 | x <- c(0.1494, 1.4849, -0.7284, 2.2843, 2.5131, -3.8016, 0.6135)                         |
| a1878299 | x <- c(-0.255, 3.4432, -0.9247, 0.8383, 1.8246, -0.2636, 0.5293)                         |
| a1879743 | x <- c(-1.3717, 2.1238, 2.785, 0.0893, -0.4422, 1.8277)                                  |
| a1880903 | x <- c(0.5291, 3.2479, -0.0146, -0.129, -1.3898, 0.1013, 0.8595, 0.8511, 1.0816, 0.5843) |
| a1881885 | x <- c(0.7392, -0.4325, 0.9298, -0.1884, -1.6, 2.7959, 2.9149)                           |
| a1882601 | x <- c(1.9765, 1.8482, 0.5765, -0.7912, 0.0571, 2.2725, 0.1549)                          |
| a1893135 | x <- c(0.6247, -1.22, 1.3306, 1.391, -0.5297)                                            |
| a1893784 | x <- c(0.1346, 2.5644, 3.8117, -0.8703, 2.5529)                                          |
| a1894438 | x <- c(1.6571, 0.6298, 0.987, 3.6434)                                                    |
| a1894742 | x <- c(-4.2777, -0.5161, 1.8382, -1.8204, 0.4913, -0.4342, 2.0994, -1.3816)              |
| a1895152 | x <- c(2.1686, 1.9082, 1.5076, -0.2429, -0.7984, 2.1778, 2.802, 0.3417, 2.3734)          |
| a1895218 | x <- c(0.0324, -1.5352, -0.9229, 1.1434)                                                 |
| a1896450 | x <- c(-1.647, 1.2667, 0.5816, 2.5689, -0.0199, -0.0102, -0.0266, 0.341)                 |
| a1896652 | x <- c(-0.2389, 0.2398, 0.0642, 1.5957, 2.5635, 2.7395)                                  |
| a1899420 | x <- c(-1.558, -1.943, 1.7205, 0.7086, 0.3317, 1.361, -0.1501, 3.0054, 0.5291)           |
| a1901571 | x <- c(0.1868, 2.7585, 0.8184, 1.3255, 0.8488)                                           |

