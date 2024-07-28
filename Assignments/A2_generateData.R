# This is to assign each student to a different species for Assignment 2

library(AnnotationHub)
library(pander)
library(taxize)
library(magrittr)
library(tidyverse)

# Q1

# Get the species with a gtf in release 100
ah <- AnnotationHub()
sp <- ah %>%
    subset(dataprovider == "Ensembl") %>%
    subset(rdataclass == "GRanges") %>%
    query(".gtf") %>%
    mcols() %>%
    as.data.frame() %>%
    filter(grepl("release-100", rdatapath)) %>%
    distinct(species, taxonomyid) %>%
    as_tibble()

# Load the students
ids <- read_tsv("../VMs.tsv", comment = "#")

# Define a function to map taxa to the common name
sp2common <- function(sp) {
    out <- sci2comm(get_uid(sp))[[1]]
    out <- c(out, "") # Add a blank field in case it comes back empty
    out[[1]]
}


# Map each student to a species
set.seed(7005)
id2sp <- tibble(ID = paste0(ids$Emplid)) %>%
    distinct(ID) %>%
    mutate(species = sample(sp$species, size = nrow(.), replace = FALSE)) %>%
    left_join(sp) %>%
    mutate(
        Name = vapply(species, sp2common, character(1)),
        Name = str_to_title(Name),
        Name = case_when(
            species == "Nothoprocta perdicaria" ~ "Chilean tinamou",
            species != "Nothoprocta perdicaria" ~ Name
        )
    )

id2sp %>%
    rename(
        Species = species,
        `Taxonomy ID` = taxonomyid,
        `Common Name` = Name
    ) %>%
    arrange(ID) %>%
    # mutate(ID = paste0("a", ID)) %>%
    pander( justify = "llrl", split.tables = Inf, style = "rmarkdown")


## | ID       | Species                      | Taxonomy ID | Common Name                    |
## |:---------|:-----------------------------|------------:|:-------------------------------|
## | a1137364 | Monodelphis domestica        |       13616 | Gray Short-Tailed Opossum      |
## | a1608942 | Neolamprologus brichardi     |       32507 | Princess cichlid               |
## | a1645191 | Mola mola                    |       94237 | Ocean Sunfish                  |
## | a1703423 | Pelodiscus sinensis          |       13735 | Chinese Soft-Shelled Turtle    |
## | a1749842 | terrapene carolina triunguis |      158814 | Three-Toed Box Turtle          |
## | a1766804 | Hippocampus comes            |      109280 | Tiger Tail Seahorse            |
## | a1773594 | bison bison bison            |        9901 | Plains Bison                   |
## | a1792812 | canis lupus dingo            |        9612 | Dingo                          |
## | a1795973 | sus scrofa pietrain          |        9823 | Pig                            |
## | a1797428 | chrysemys picta bellii       |        8479 | Western Painted Turtle         |
## | a1823995 | Aotus nancymaae              |       37293 | Ma's Night Monkey              |
## | a1839745 | Poecilia mexicana            |       48701 | Shortfin molly                 |
## | a1843355 | mustela putorius furo        |        9668 | Domestic Ferret                |
## | a1850508 | Serinus canaria              |        9135 | Common Canary                  |
## | a1851176 | Salmo salar                  |        8030 | Atlantic Salmon                |
## | a1851451 | Poecilia reticulata          |        8081 | Guppy                          |
## | a1851815 | Macaca fascicularis          |        9541 | Crab-Eating Macaque            |
## | a1853428 | Poecilia formosa             |       48698 | Amazon Molly                   |
## | a1862759 | Pan paniscus                 |        9597 | Pygmy Chimpanzee               |
## | a1863615 | Neovison vison               |      452646 | American Mink                  |
## | a1864525 | Fukomys damarensis           |      885580 | Damara Mole-Rat                |
## | a1865749 | Gambusia affinis             |       33528 | Western Mosquitofish           |
## | a1868878 | sus scrofa usmarc            |        9823 | Pig                            |
## | a1869981 | Chinchilla lanigera          |       34839 | Long-Tailed Chinchilla         |
## | a1871922 | Rhinolophus ferrumequinum    |       59479 | Greater Horseshoe Bat          |
## | a1872208 | Monopterus albus             |       43700 | Swamp Eel                      |
## | a1873151 | mus musculus c57bl6nj        |       10090 | Mouse                          |
## | a1875420 | Fundulus heteroclitus        |        8078 | Mummichog                      |
## | a1876866 | Parambassis ranga            |      210632 | Indian Glassy Fish             |
## | a1886661 | Castor canadensis            |       51338 | American Beaver                |
## | a1893982 | astyanax mexicanus pachon    |        7994 | Pachon Cavefish                |
## | a1894721 | Pygocentrus nattereri        |       42514 | Red-Bellied Piranha            |
## | a1894991 | Oryzias javanicus            |      123683 | Javanese Ricefish              |
## | a1897552 | Clupea harengus              |        7950 | Atlantic Herring               |
## | a1899345 | Ictidomys tridecemlineatus   |       43179 | Thirteen-Lined Ground Squirrel |
## | a1900426 | Sarcophilus harrisii         |        9305 | Tasmanian Devil                |
## | a1901379 | Otolemur garnettii           |       30611 | Small-Eared Galago             |
## | a1903005 | Xiphophorus maculatus        |        8083 | Southern Platyfish             |
## | a1904204 | Anabas testudineus           |       64144 | Climbing Perch                 |
## | a1904509 | Lepidothrix coronata         |      321398 | Blue-Crowned Manakin           |
## | a1906661 | mus musculus akrj            |       10090 | Mouse                          |
## | a1908426 | Strigops habroptila          |     2489341 | Kakapo                         |
## | a1909154 | Cebus capucinus              |        9516 | White-Faced Sapajou            |
## | a1910059 | sus scrofa bamei             |        9823 | Pig                            |
## | a1930255 | sus scrofa berkshire         |        9823 | Pig                            |
## | a1932615 | Zonotrichia albicollis       |       44394 | White-Throated Sparrow         |
## | a1935423 | Geospiza fortis              |       48883 | Medium Ground-Finch            |
## | a1940870 | Pogona vitticeps             |      103695 | Central Bearded Dragon         |
## | a1947736 | Equus caballus               |        9796 | Horse                          |
## | a1947841 | Larimichthys crocea          |      215358 | Large Yellow Croaker           |
## | a1954027 | Ictalurus punctatus          |        7998 | Channel Catfish                |
## | a1954456 | panthera tigris altaica      |        9694 | Amur Tiger                     |
## | a1955686 | Ovis aries                   |        9940 | Sheep                          |


# Q4:

ids$Emplid %>%
    sort() %>%
    unique() %>%
    str_remove("a") %>%
    na.omit() %>%
    lapply(function(x) {
        set.seed(as.integer(x))
        n <- sample(4:10, 1)
        tibble(
            ID = paste0("a", x),
            Values = paste0(
                "x <- c(",
                rnorm(n, 0.7, 1.5) %>% round(4) %>% paste0(collapse = ", "),
                ")")
        )
    }) %>%
    bind_rows() %>%
    pander(
        split.tables = Inf,
        style = "rmarkdown",
        justify = "ll"
    )

## | ID       | Values                                                                                       |
## |:---------|:---------------------------------------------------------------------------------------------|
## | a1137364 | x <- c(-0.0433, 2.2891, 1.0951, 1.8736)                                                      |
## | a1608942 | x <- c(-0.1669, -0.1257, 0.5036, 0.0166, -2.4042, -0.1197, 2.1225, 1.212, 5.6835)            |
## | a1645191 | x <- c(1.4484, 1.3787, -0.4143, 3.6535, 1.8083, 0.4871)                                      |
## | a1703423 | x <- c(0.8609, 0.1178, 0.9691, -1.4772, 2.5688, 0.1644, -0.8474, -2.1161)                    |
## | a1749842 | x <- c(1.6882, 0.6029, 1.9638, -0.8317, -0.988, -2.3742, -0.4314)                            |
## | a1766804 | x <- c(4.8103, -0.8988, -0.7313, 3.5593, -0.4826, 2.4308, -0.4153, 0.0668, 2.867)            |
## | a1773594 | x <- c(0.2608, -1.4151, 2.8979, 1.4756, 2.2516)                                              |
## | a1792812 | x <- c(-0.6245, -0.2226, 1.5684, -1.2285)                                                    |
## | a1795973 | x <- c(1.8684, 0.2476, 0.5143, 0.3547, -0.1645, 1.251, -1.9733, 2.3365, -0.1596)             |
## | a1797428 | x <- c(-1.2975, 1.7608, -0.0982, 1.8631, -2.3649, 0.1515, 1.709, -1.1104)                    |
## | a1823995 | x <- c(0.0863, 0.5653, 0.4356, 2.2952, 0.8758, 0.0343)                                       |
## | a1839745 | x <- c(0.5774, 1.3214, 2.6437, -0.2024)                                                      |
## | a1843355 | x <- c(1.0581, -1.0416, -0.8031, -0.1871)                                                    |
## | a1850508 | x <- c(-0.9714, 1.6317, -1.3517, 1.1421, -0.2032, 1.0515, 1.9313, 1.8357, -1.4458)           |
## | a1851176 | x <- c(-2.08, -0.8782, 0.285, 1.3781)                                                        |
## | a1851451 | x <- c(1.2449, -0.7304, -0.1646, 2.0851, -0.7812, 0.5486, 0.9953, 0.1769, 0.3606)            |
## | a1851815 | x <- c(0.0134, 1.9545, 1.0156, 2.4974)                                                       |
## | a1853428 | x <- c(-0.0298, -1.3617, 2.1709, 0.8543, 0.2726, -0.7025)                                    |
## | a1862759 | x <- c(0.0065, 1.0768, -0.7598, 3.0874)                                                      |
## | a1863615 | x <- c(0.7639, 0.0844, 1.6849, -0.5075, -0.8914, 1.0112, 2.7931)                             |
## | a1864525 | x <- c(1.1527, 1.6518, 1.6833, -1.3038, 0.6058, 0.4401, -1.1379, 1.6957)                     |
## | a1865749 | x <- c(1.5056, 1.647, -0.0661, -2.1343, 1.423, 0.7569, 0.0609)                               |
## | a1868878 | x <- c(-2.9709, -0.45, 0.0074, 0.1217, 2.5973)                                               |
## | a1869981 | x <- c(0.8152, 1.4271, 0.4694, -0.2918, -0.286, -0.6054, 1.2336)                             |
## | a1871922 | x <- c(0.5119, -1.8475, 1.8198, -0.1338, 0.8777, 0.0604, -0.4013, 0.3594, 2.3255)            |
## | a1872208 | x <- c(0.7079, 0.3685, 0.2181, 0.5549, -2.3652, -1.2004, 3.4151)                             |
## | a1873151 | x <- c(2.4517, 1.4722, 2.6803, -1.9407, -0.6271, -0.2569, -0.4746, 2.6373, -2.0029)          |
## | a1875420 | x <- c(-1.5192, 2.833, 3.3837, -1.0578, 0.8613, 2.5236, -2.2051, -0.4782, 3.1518)            |
## | a1876866 | x <- c(-1.123, -0.3904, -0.4987, 2.0421, -0.8621, -0.7427, 0.7549, -0.5748, 0.4906, -1.0987) |
## | a1886661 | x <- c(1.129, -0.5184, 0.0759, 0.4659)                                                       |
## | a1893982 | x <- c(0.4595, 1.7962, 2.1128, 2.646, -2.3673, -0.1961)                                      |
## | a1894721 | x <- c(-2.91, 2.5906, -1.8167, -0.9061, -0.8382, 0.2453, 0.4134, 3.5842)                     |
## | a1894991 | x <- c(1.9851, -0.5493, 1.7286, 3.7736, 0.4708, 0.7642, 0.8895, 2.0142, -0.3351)             |
## | a1897552 | x <- c(3.2391, -1.7528, 0.4305, -0.6628, -1.6469, 0.782, 1.5166, 0.237)                      |
## | a1899345 | x <- c(0.1474, 1.9407, 0.5564, 1.7953, 1.2655, 0.1691, 0.4547, 0.0187, -2.1591, 0.6655)      |
## | a1900426 | x <- c(-0.3584, 3.5074, 2.9401, 2.1211, -0.1264, 0.9589, -1.1208, -1.9353, 1.7136)           |
## | a1901379 | x <- c(0.3038, -1.1985, -1.7394, 0.2135, 2.5004, 1.6915, 2.8043, 0.5944)                     |
## | a1903005 | x <- c(0.6722, -0.1956, 0.5278, 1.8678, 3.1052, 1.2281, 3.0446, -0.5916)                     |
## | a1904204 | x <- c(-0.992, 4.0864, 1.3728, 1.2436, -0.6282)                                              |
## | a1904509 | x <- c(-0.1648, 2.1814, -1.0161, 0.9513, 1.0205, -0.0862, 0.5003, -1.3463, -0.1135, 2.1983)  |
## | a1906661 | x <- c(0.706, -0.9058, 2.1544, 2.9012, 1.1209, -1.3753, 1.9729, 2.875, 3.6917, -1.077)       |
## | a1908426 | x <- c(1.2641, 1.44, 0.1889, 0.1464, -0.4497, 1.3661, 2.1486, 2.4283, 1.108, 1.2788)         |
## | a1909154 | x <- c(2.4852, 1.2269, -0.0679, -0.9611, 2.1799)                                             |
## | a1910059 | x <- c(2.8994, -3.0528, 0.1018, 1.3461, 0.9224, 1.1955, 2.3327)                              |
## | a1930255 | x <- c(2.0478, -1.8305, 1.4769, -0.1325, 1.2525, 1.7385)                                     |
## | a1932615 | x <- c(0.9234, 1.4015, -0.3947, -3.1329, -0.1542, 1.9771, -0.4997, 0.0226)                   |
## | a1935423 | x <- c(0.6379, 2.0681, 2.887, 0.844, 0.398, 0.265, 0.0612, 0.8686)                           |
## | a1940870 | x <- c(-1.0123, 2.6733, 0.0865, 2.7096, 2.1062, -0.157, 0.6258, -0.5091, 0.114, 0.0341)      |
## | a1947736 | x <- c(1.2537, 2.7525, 0.5015, 0.3354, -1.129, 1.1688, 3.6342, 1.2302, 0.8297, 1.3461)       |
## | a1947841 | x <- c(-1.5333, 2.0824, 2.7886, -0.0431, 0.0146, -1.404, 0.7452, 3.0861, -1.9524, 0.3867)    |
## | a1954027 | x <- c(0.6892, -0.1318, -0.1043, 4.0503, 0.0989, 0.8987, 1.7833)                             |
## | a1954456 | x <- c(-0.0944, -1.7223, 0.9764, 0.657, 0.6068, -1.0101, 0.5493)                             |
## | a1955686 | x <- c(0.1344, 0.0336, -1.2412, 0.3731)                                                      |
