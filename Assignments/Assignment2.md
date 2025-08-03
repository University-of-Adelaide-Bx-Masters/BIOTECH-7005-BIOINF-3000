# Assignment 2 [*30 marks*] 

## Note that failure to include the coversheet in the R markdown document will incur a 5 mark penalty.

**Due before 12pm, Friday 29th August, 2025**

Your answers to all questions should be submitted to myUni as a single `.zip` file containing 3 files:
1) a bash script for Q1 2) a bash script for Q2, and 3) a single Rmarkdown document that includes the required coversheet and contains the answers to the statistics questions (Q3 and Q4).
Note that the file `my_species_gff_features.txt` is not required as part of your submission for Q1, **only the script which will generate this file**!
Similarly, for Q2, only the script is required.

## Required scripts [*15 marks*]

**Q1.** Write a script to:
    + Download the gff3 file for your assigned species ([see bottom of page](#species-for-question-1)) to your current directory from Ensembl [*1 mark*]
    + Count how many of each feature type there is, sorted in numerical order [*4 marks*]
    + Export the results to a file with a name of the form `my_species_gff_features.txt` **where you use your assigned species name instead of** `my_species` [*1 mark*].
    NB: If your actual species is not included in the name, no marks will be given.
    + The script must also include code to generate one or more comment lines in the output file/table before the table with the genome-build used, (hint: grep your gff to find the genome build info as the header is very large in most cases)
    + The script must also write the code used to generate the summary (counts) data to the output file as part of the file header. [*2 marks*]
    
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

| ID       | Species                         | Taxonomy ID | Common Name                    |
|:---------|:--------------------------------|------------:|:-------------------------------|
| a1708771 | Geospiza fortis                 |       48883 | Medium Ground-Finch            |
| a1720668 | Pan paniscus                    |        9597 | Pygmy Chimpanzee               |
| a1755331 | Hippocampus comes               |      109280 | Tiger Tail Seahorse            |
| a1763968 | equus asinus asinus             |        9793 |                                |
| a1793955 | Xiphophorus maculatus           |        8083 | Southern Platyfish             |
| a1819179 | Stachyris ruficeps              |      181631 | Rufous-Capped Babbler          |
| a1822446 | Equus caballus                  |        9796 | Horse                          |
| a1824288 | Clupea harengus                 |        7950 | Atlantic Herring               |
| a1836916 | Ovis aries                      |        9940 | Sheep                          |
| a1843731 | Gambusia affinis                |       33528 | Western Mosquitofish           |
| a1845468 | Gasterosteus aculeatus          |       69293 | Three-Spined Stickleback       |
| a1849952 | Zonotrichia albicollis          |       44394 | White-Throated Sparrow         |
| a1852174 | Sarcophilus harrisii            |        9305 | Tasmanian Devil                |
| a1853492 | Macaca fascicularis             |        9541 | Crab-Eating Macaque            |
| a1864358 | Gorilla gorilla                 |        9593 | Western Gorilla                |
| a1869687 | Poecilia formosa                |       48698 | Amazon Molly                   |
| a1872943 | mustela putorius furo           |        9668 | Domestic Ferret                |
| a1878445 | astyanax mexicanus pachon       |        7994 |                                |
| a1879522 | bison bison bison               |        9901 |                                |
| a1879968 | Ictidomys tridecemlineatus      |       43179 | Thirteen-Lined Ground Squirrel |
| a1883243 | Chinchilla lanigera             |       34839 | Long-Tailed Chinchilla         |
| a1885087 | Otolemur garnettii              |       30611 | Small-Eared Galago             |
| a1885559 | Cebus capucinus                 |        9516 | White-Faced Sapajou            |
| a1885715 | terrapene carolina triunguis    |      158814 | Three-Toed Box Turtle          |
| a1886568 | Mus pahari                      |       10093 | Shrew Mouse                    |
| a1886602 | Neovison vison                  |      452646 | American Mink                  |
| a1887315 | Serinus canaria                 |        9135 | Common Canary                  |
| a1889389 | canis lupus familiarisgreatdane |        9612 |                                |
| a1889603 | Seriola dumerili                |       41447 | Greater Amberjack              |
| a1892205 | mus musculus c57bl6nj           |       10090 |                                |
| a1894721 | Electrophorus electricus        |        8005 | Electric Eel                   |
| a1895247 | Monopterus albus                |       43700 | Swamp Eel                      |
| a1896083 | Pogona vitticeps                |      103695 | Central Bearded Dragon         |
| a1896180 | Anas platyrhynchos              |        8839 | Mallard                        |
| a1897497 | sus scrofa usmarc               |        9823 |                                |
| a1897522 | Mola mola                       |       94237 | Ocean Sunfish                  |
| a1901527 | mus musculus lpj                |       10090 |                                |
| a1902999 | Rattus norvegicus               |       10116 | Norway Rat                     |
| a1903474 | Poecilia reticulata             |        8081 | Guppy                          |
| a1904400 | Oryzias javanicus               |      123683 | Javanese Ricefish              |
| a1907116 | Salmo salar                     |        8030 | Atlantic Salmon                |
| a1909495 | Anabas testudineus              |       64144 | Climbing Perch                 |
| a1911607 | Vulpes vulpes                   |        9627 | Red Fox                        |
| a1914291 | Notechis scutatus               |        8663 | Mainland Tiger Snake           |
| a1917758 | Rhinolophus ferrumequinum       |       59479 | Greater Horseshoe Bat          |
| a1919050 | chrysemys picta bellii          |        8479 | Western Painted Turtle         |
| a1929378 | Neolamprologus brichardi        |       32507 |                                |
| a1930651 | Castor canadensis               |       51338 | American Beaver                |
| a1930685 | Lepidothrix coronata            |      321398 | Blue-Crowned Manakin           |
| a1933020 | panthera tigris altaica         |        9694 | Amur Tiger                     |
| a1933380 | Pelodiscus sinensis             |       13735 | Chinese Soft-Shelled Turtle    |
| a1937503 | Sphenodon punctatus             |        8508 | Tuatara                        |
| a1940517 | Strigops habroptila             |     2489341 | Kakapo                         |
| a1944175 | Vombatus ursinus                |       29139 | Common Wombat                  |
| a1945399 | Aotus nancymaae                 |       37293 | Mas Night Monkey              |
| a1945739 | Parambassis ranga               |      210632 | Indian Glassy Fish             |
| a1948226 | sus scrofa berkshire            |        9823 |                                |
| a1948692 | sus scrofa bamei                |        9823 |                                |
| a1949596 | heterocephalus glaber male      |       10181 |                                |
| a1951109 | Callithrix jacchus              |        9483 | White-Tufted-Ear Marmoset      |
| a1951803 | Panthera pardus                 |        9691 | Leopard                        |
| a1955540 | sus scrofa pietrain             |        9823 | Pig                            |
| a1957056 | Ictalurus punctatus             |        7998 | Channel Catfish                |
| a1959486 | Pygocentrus nattereri           |       42514 | Red-Bellied Piranha            |
| a1963706 | Myripristis murdjan             |      586833 | Pinecone Soldierfish           |
| a1974970 | Erinaceus europaeus             |        9365 | Western European Hedgehog      |
| a1975156 | mus musculus akrj               |       10090 |                                |
| a1975386 | Spermophilus dauricus           |       99837 | Daurian Ground Squirrel        |
| a1976637 | Monodelphis domestica           |       13616 | Gray Short-Tailed Opossum      |
| a1979858 | canis lupus dingo               |        9612 | Dingo                          |
| a1983144 | Fukomys damarensis              |      885580 | Damara Mole-Rat                |
| a1983151 | Fundulus heteroclitus           |        8078 | Mummichog                      |
| a1983611 | Larimichthys crocea             |      215358 | Large Yellow Croaker           |
| a1986777 | Poecilia mexicana               |       48701 |                                |
| a1990488 | Gadus morhua                    |        8049 | Atlantic Cod                   |
| a1991182 | Physeter catodon                |        9755 | Sperm Whale                    |

## Values For Question 4

*If your student number is not listed, please contact Dave to ensure you are added to the list*

The results you are analysing for Q4 are as follows.
You can simply paste these values into your RMarkdown document as the object `x` and perform all of your analysis on these values.

| ID       | Values                                                                                    |
|:---------|:------------------------------------------------------------------------------------------|
| a1708771 | x <- c(1.1676, 1.5923, -1.4217, 0.5851, -0.2704)                                          |
| a1720668 | x <- c(-0.4687, -0.8451, 0.6906, 2.2823, 1.4765, 0.9923, -3.455, 1.6946, 1.2667, 1.0236)  |
| a1755331 | x <- c(0.9452, -1.7696, -0.3727, -0.0979, 0.5222)                                         |
| a1763968 | x <- c(-1.4874, 2.817, -0.0229, 2.1672, 3.0935, 2.8635)                                   |
| a1793955 | x <- c(-0.6952, 1.9843, 0.3385, 0.7023, -1.4734, 0.6561, -0.1272, 3.7022, 0.0056)         |
| a1819179 | x <- c(2.4937, 2.5255, 0.6468, 0.6671, -0.7512, 0.34, -1.2484, -2.2171, -0.5011, -0.5522) |
| a1822446 | x <- c(1.565, -0.6666, 3.0865, -1.1691, -3.1269, 0.2802, 0.6707, -1.102)                  |
| a1824288 | x <- c(-1.368, -0.4883, -1.3923, 0.5309, 1.4886, -0.2512)                                 |
| a1836916 | x <- c(-1.0667, 2.0526, -0.8862, 1.0869, 2.1607, 4.5185, 0.3979)                          |
| a1843731 | x <- c(0.1426, -1.142, -0.9352, 2.1866, 3.4501, 0.3307, 1.4802, 1.0241, 3.0194, 0.2908)   |
| a1845468 | x <- c(2.0259, 1.0872, -0.6666, 2.0321, 1.2129, 0.345, 0.777, -0.0602, 2.9802, 1.2011)    |
| a1849952 | x <- c(-0.562, 0.507, 0.4156, 1.9321, -0.6901, 2.5057, 2.9341)                            |
| a1852174 | x <- c(-0.0818, 1.2866, -2.7919, 0.9343, 1.8301, 0.0505)                                  |
| a1853492 | x <- c(1.3008, 1.6037, -1.9613, 1.1027, 0.0643, 1.2814, 2.5214, -0.4455, 0.943, 0.0665)   |
| a1864358 | x <- c(2.8079, -0.7271, -0.0631, 0.939, -0.5395, 1.8911, -1.1499)                         |
| a1869687 | x <- c(1.975, 0.9456, -0.0759, 1.9273, 0.7749, 1.1975)                                    |
| a1872943 | x <- c(0.0049, -2.0333, 0.2881, 1.2394)                                                   |
| a1878445 | x <- c(-0.7838, -0.872, -1.5637, -0.445, -1.3435, 0.4101, -0.3361, 2.6491)                |
| a1879522 | x <- c(-2.8408, 2.9516, -0.9095, 0.6706, 1.155, 0.4556, 1.0667, 0.086)                    |
| a1879968 | x <- c(-3.2873, 2.2309, 1.2875, 2.5734)                                                   |
| a1883243 | x <- c(-0.7617, -0.4972, -0.3728, 1.9015, 0.5713, 2.4009, 0.6329)                         |
| a1885087 | x <- c(0.9291, 2.8866, 2.6037, -0.1279, 0.6521, 1.5384, 0.5273, 0.7336)                   |
| a1885559 | x <- c(0.0304, 1.365, 0.0957, 2.0407, 0.981, -0.5347, -0.0397, 0.4169, 0.3015, -0.5497)   |
| a1885715 | x <- c(0.2798, -0.0197, -0.4017, 0.5075, 2.013, 1.7263)                                   |
| a1886568 | x <- c(1.1486, 0.4362, -0.0242, -0.9989, 1.0615, 0.0905)                                  |
| a1886602 | x <- c(-0.3897, 1.3611, 0.7487, 1.1626, 0.5217, 0.8723)                                   |
| a1887315 | x <- c(0.3498, 0.0528, 2.1724, 0.8877)                                                    |
| a1889389 | x <- c(0.3159, 2.4646, 1.0241, -0.6489, 0.592, 1.7307)                                    |
| a1889603 | x <- c(0.5761, 1.0248, 1.4891, 0.946)                                                     |
| a1892205 | x <- c(1.1897, -0.2483, 2.3366, 3.6815, 0.3369, 1.0923, -1.5599, 1.4247, 0.5023)          |
| a1894721 | x <- c(-2.91, 2.5906, -1.8167, -0.9061, -0.8382, 0.2453, 0.4134, 3.5842)                  |
| a1895247 | x <- c(-0.8024, 0.4674, 1.3285, -3.3951, 0.4212, 0.5936, 1.4575)                          |
| a1896083 | x <- c(-0.2086, 2.0357, 0.5334, -1.5515, 0.1, 1.6114)                                     |
| a1896180 | x <- c(1.2173, -0.4234, 0.2251, 1.7039)                                                   |
| a1897497 | x <- c(1.9912, -0.2666, 1.6265, 3.2852, 3.4556)                                           |
| a1897522 | x <- c(-0.2623, -0.207, -2.2133, 0.7385, -0.1388, 0.5881, -2.4479, 1.5216)                |
| a1901527 | x <- c(2.8694, 3.0356, 3.1017, 0.5034, 1.3967, -0.9433, 0.3333, 4.6196, -0.4135)          |
| a1902999 | x <- c(0.6411, -1.8589, 0.8354, 0.0432, -0.6618, -0.4525, -0.2123, -2.3073, 0.2297)       |
| a1903474 | x <- c(-1.3106, 1.6855, -0.8576, 1.3556)                                                  |
| a1904400 | x <- c(0.617, 0.2247, 2.1808, 1.5618, 1.7485, 0.1866, 0.7795)                             |
| a1907116 | x <- c(-0.1745, -1.5766, -3.4593, -1.0407, 1.2241, 0.6243, 0.5142, -1.2251)               |
| a1909495 | x <- c(-0.8178, 4.6328, 0.5662, 0.572, -0.1946, 2.2051)                                   |
| a1911607 | x <- c(-1.4294, 0.5704, -0.5816, -0.0726, -1.9107, -0.6052, 0.618, 0.8001, 0.469)         |
| a1914291 | x <- c(1.7361, 2.5607, 2.1167, -0.367, 1.9165, 0.4292, 2.7977)                            |
| a1917758 | x <- c(1.0694, -2.4646, 4.2174, -1.3214, -0.169, 0.3416, -0.6334, 4.2043, -0.182, 0.5595) |
| a1919050 | x <- c(2.8964, 3.435, 0.7892, 0.3828)                                                     |
| a1929378 | x <- c(-0.4051, 1.4586, 1.3855, -0.4218, 0.6259, -0.9671, 1.1158)                         |
| a1930651 | x <- c(-0.9146, 1.4609, 4.0118, 0.8326, 2.2514, 0.8556, 1.0982, -2.9615, 3.3992)          |
| a1930685 | x <- c(0.1906, 1.6408, 0.3808, 0.7056, -1.5742, 1.357, 0.3294, 1.7013)                    |
| a1933020 | x <- c(-0.629, 1.7442, -1.2879, 0.356)                                                    |
| a1933380 | x <- c(2.4778, -1.8242, -1.3432, 2.0368, -0.1312, 1.137, 1.0001, -0.98, 1.0707, 1.1617)   |
| a1937503 | x <- c(2.6393, 2.2273, 1.4009, 0.2748, -1.1967, 0.4698, 3.5332, 0.2342)                   |
| a1940517 | x <- c(1.6712, 1.5163, -1.7917, -2.063, 3.1766, 1.0713, 0.5973, 3.4403, 0.8431, 1.878)    |
| a1944175 | x <- c(1.5373, 1.1169, 2.7838, -1.2507, 2.9092, 0.7539)                                   |
| a1945399 | x <- c(0.0553, 1.634, 1.1557, -0.2745, 0.9093, 2.0433, 2.7048)                            |
| a1945739 | x <- c(-1.3658, -0.5366, -0.0791, 2.4061, 0.2482, -1.3245, 0.5526, 0.2848)                |
| a1948226 | x <- c(1.2493, -0.4266, 0.7263, 3.5066)                                                   |
| a1948692 | x <- c(0.1302, 2.7828, -0.8858, 1.9182)                                                   |
| a1949596 | x <- c(-2.6293, -1.271, 0.0546, -0.0771, 0.453, -0.3709, 1.8781, 1.0025)                  |
| a1951109 | x <- c(-0.4134, 1.1082, 1.9014, 1.1843, -0.2008, 0.6032)                                  |
| a1951803 | x <- c(1.1136, 2.1454, -2.166, -1.9016, -0.465, -1.2018)                                  |
| a1955540 | x <- c(3.332, 2.8928, 0.4132, 0.0113, -1.9271, -0.245, -2.724, -1.0531)                   |
| a1957056 | x <- c(1.8228, 0.5039, 0.4556, -1.7818)                                                   |
| a1959486 | x <- c(1.664, 0.4377, 3.0273, 0.2671, 1.9161, 0.3188, -1.1436, 0.6747)                    |
| a1963706 | x <- c(2.4668, -0.5662, -0.3482, 1.3086)                                                  |
| a1974970 | x <- c(0.9952, 1.5271, -0.7872, -0.4163, 1.7044, 1.5675)                                  |
| a1975156 | x <- c(0.647, 2.6374, 1.4957, -1.065, 0.0135, 0.9686, -3.1608, 1.1678, 0.7309)            |
| a1975386 | x <- c(-0.0103, -0.4479, -3.0579, 2.1068, 1.2826, 1.2599, 2.5973, -0.1294, -0.6084)       |
| a1976637 | x <- c(-0.0315, -1.4609, 3.9244, 0.9584, 2.0404, 0.9408, 2.0109, 0.9732, 1.1352, 0.111)   |
| a1979858 | x <- c(0.0789, -0.338, -1.3878, 0.157, -0.5403, 0.4008)                                   |
| a1983144 | x <- c(2.846, 0.726, 0.9368, 0.6942, 0.9697, -1.0579, 0.6473, 1.1155, -1.6974)            |
| a1983151 | x <- c(-0.9731, 2.5764, 0.8868, -0.6342, 0.0407, -1.0549, 0.3435, 3.1698)                 |
| a1983611 | x <- c(-0.3212, 0.6662, 2.3223, 3.7371, 0.7526, 2.2616)                                   |
| a1986777 | x <- c(-1.22, -1.0838, -0.3691, 0.9494, -1.5468, 1.3327, 1.2183)                          |
| a1990488 | x <- c(-1.1542, 0.6236, 1.3384, -1.7706, -0.2321, 1.8568, -1.6587, 1.5129, 2.1126)        |
| a1991182 | x <- c(0.315, 0.5431, -0.5579, 0.2079)                                                    |

