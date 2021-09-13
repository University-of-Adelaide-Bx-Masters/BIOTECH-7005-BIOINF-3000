##############
## Packages ##
##############
install.packages("vcfR")
library(magrittr)
library(vcfR)
library(tidyverse)

###############
## Functions ##
###############
'%nin%' <- Negate('%in%')
ulength <- function(x) {x %>% unique() %>% length()}
fsummary <- function(x) {x %>% as.factor() %>% summary()}

############
## Set Up ##
############
DirIn <- "/home/student/Practical_8"

DirOut <- "/home/student/Practical_8/Out"

if(! dir.exists(DirOut)) {dir.create(DirOut)}

DirPlot <- file.path(DirOut, "Plots")
if(! dir.exists(DirPlot)) {dir.create(DirPlot)}

#########
## IGV ##
#########
# Load ref and then Bam and VCF files in to IGV

# Load for now:
# 1:15,448,200-15,448,400 # 2 SNPs

# Look at later:
# 1:15,437,129-15,437,259 # SNPs, Del, MNP, tri-allelic SNP
# 1:15,488,098-15,488,162 # Weak Indel
# 1:15,522,627-15,522,757 # Uncalled variants
# 1:16,509,693-16,526,713 # Highly variable region
# 1:17,074,700-17,210,875 # Highly variable region
# 1:20,054,605-20,190,780 # repeated region

###################
## Load GFF File ##
###################
# DO NOT EXECUTE THIS SECTION
# FileIn <- "Arabidopsis_thaliana.TAIR10.48.gff3.gz"
#
# GFF <- read.table(file.path(DirIn, FileIn), sep="\t", quote="") #%>% as_tibble()
#
# GFF %>% dim() # 791564      9

#########################
## Load Reference File ##
#########################
# DO NOT EXECUTE THIS SECTION
# FileIn <- "Arabidopsis_thaliana.TAIR10.dna.toplevel.fa"
#
# Ref <- ape::read.dna(file.path(DirIn, FileIn), format = "fasta")

###############
## Load VCF ##
###############
# Load file
FileIn <- "SRR5882792_Athaliana_TAIR10_picard.vcf"
VCF <- read.vcfR(file.path(DirIn, FileIn), verbose = TRUE)

# A quick look at components
VCF@meta %>% as_tibble() # metadata
VCF@fix %>% as_tibble() # co-ordinates etc and INFO ie data about this position in the genome
VCF@gt %>% as_tibble() # FORMAT ie genotype(s)

INFO2df(VCF[1:20,]) %>% as_tibble() # converts INFO key value pairs to a data frame

###########################
## Extract Required Data ##
###########################
# Keys for genotypes are stored in FORMAT column
# Values for genotypes are stored in a separate column for each sample - only one sample reported here
VCF@gt %>% as_tibble() #

# Convert FORMAT to DF, add variant co-ordinates etc
VCF@gt[,1] %>% unique() # all keys are present in all records, usually not the case but here it makes analysis easier

ColNames <- VCF@gt[1,1] %>% strsplit(":") %>% unlist() %>% unname() # get column names
ColNames

GenomeData <- VCF@gt %>% as_tibble %>% dplyr::select(-FORMAT) %>% separate(unknown, ColNames, sep = ":") # convert strings to cols
GenomeData$DP %<>% as.integer()
GenomeData$RO %<>% as.integer()
GenomeData$QR %<>% as.integer()
GenomeData$AO %<>% as.integer()
GenomeData$QA %<>% as.integer()

GenomeData %<>% select(GT, DP, AD, RO, AO, QR, QA, GL) # get the column order we want

# View col definitions and data to see what each column means
VCF@meta %>% as_tibble() %>% filter(grepl("FORMAT", value))
GenomeData

# Add the genotype quality (GQ), the difference between best and next best genotype likelihood - best is always 0, just need to get 2nd best
GenomeData %<>% rowwise() %>%
    mutate(GQ = GL %>%
               strsplit(",") %>% # split string into bits using comma to separate
               unlist() %>% # convert to a vector
               as.numeric() %>% # make it numeric
               sort() %>% # sort
               .[2] %>%  # take the second element
               abs()) %>% # get absolute value
    ungroup()

GenomeData %<>% select(GT, GQ, everything()) # get the column order we want

# Get INFO info we want to keep
X <- VCF@fix %>% as_tibble() %>% dplyr::select(CHROM, POS, REF, ALT) # get required INFO & other cols
colnames(X) <- c("Chr", "Pos", "Ref", "Alt") # rename

GenomeData <- bind_cols(X, GenomeData) # merge INFO and FORMAT data

GenomeData

rm(X)
##################
## Analyse Data ##
##################
## Genotypes ##
# Check values for genotypes
GenomeData$GT %>% fsummary()
#   0/0   0/1   0/2   0/3   0/4   0/5   0/6   1/1   1/2   2/2
# 15796 17399   289    57     8     4     1  2991    84     1

GenomeData %>% filter(GT == "0/2") %>% select(-GL)

GenomeData %>% filter(GT == "1/2") %>% select(-GL)

# Count entries with multiple alternative alleles
grepl(",", GenomeData$Alt) %>% summary() # test for a comma in alt allele
#    Mode   FALSE    TRUE
# logical   35256    1374

# Inspect a highly variant region including a tri-allelic position
# 1:15,437,129-15,437,259 # SNPs, Del, MNP, tri-allelic SNP
GenomeData %>% filter(Chr == "1" & Pos >= "15437135" & Pos <= "15437191")

# Remove entries with multiple alternative alleles
GenomeData %>% dim() # 36630    13
GenomeData %<>% filter(!grepl(",", Alt))
GenomeData %>% dim() # 35256    13

GenomeData$GT %>% fsummary()
#   0/0   0/1   1/1
# 15063 17204  2989

## Read Depth ##
GenomeData$DP %>% summary()
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#  2.0    10.0    14.0   111.4    33.0 15980.0

GenomeData$DP %>% hist()

# Check enrties with far too many reads
GenomeData %>% filter(DP > 1000) # positions called hom ref when there is strong evidence for an alt allele

GenomeData %>% filter(DP > 1000) %>% arrange(desc(DP))

# Look at 1:15082284 in IGV
GenomeData %>% filter(Pos == 15082284)

## Positions With 20 Reads ##
X <- GenomeData %>% filter(DP == 20) %>% arrange(AO)
X %>% slice(c(23:25))

X$GQ %>% summary()
X$GQ %>% hist()
X$GQ[X$GQ < 1.7] %>% hist()

# Relationship between GT (genotype) and GQ (genotype quality)
Title <- sprintf("Arabidopsis TAI10 Subset with 20 Reads")
SubTitle <- sprintf("Genotype Quality by Genotype Call")
X %>%
    ggplot(aes(x = GT, y = GQ, colour = GT)) +
    geom_boxplot() +
    labs(title = Title, subtitle = SubTitle, x = "GT", y = "GQ") +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
          plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5),
          legend.position = "none", legend.title=element_blank())
FileOut <- sprintf("%s%s", Title, ".jpeg")
ggsave(filename = FileOut, path = DirPlot, device = "jpeg")

# Histogram of GQ for each GT
Title <- sprintf("TAI10 variants with 20 Reads")
SubTitle <- sprintf("Genotype Call Quality Score by Genotype Histogram")
X %>%
    ggplot(aes(x = GQ, colour = GT)) +
    geom_histogram(colour = "lightblue", fill = "lightblue") +
    labs(title = Title, subtitle = SubTitle, x = "GQ", y = "Count") +
    facet_wrap(~GT) +
    theme(axis.text.x = element_text(angle = -45, hjust = 0.5),
          plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5),
          legend.position = "none", legend.title=element_blank())
FileOut <- sprintf("%s%s", Title, ".jpeg")
ggsave(filename = FileOut, path = DirPlot, device = "jpeg")


## Non-reference Calls ##
# Keep only entries withh variant genotypes ie  0/1` and 1/1
X <- GenomeData %>% filter(GT %in% c("0/1", "1/1"))  # positions with alt

## How many reads support alt genotype calls ##
X$AO %>% hist()
X$AO[X$AO < 100] %>% hist()

X$AO %>% summary()
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 1.00    2.00    4.00   19.91    8.00 5309.00

# Check those with very few reads
X %>% filter(AO <= 1) %>% View()
# Should these variants be called?
# Look at 3:12247004 - why is the variant shifted?

# Check top and bottom distribution of variants with 2 alt reads
X %>% filter(AO == 2) %>% arrange(desc(GQ)) # sort variants with 2 alt reads by GQ
X %>% filter(AO == 2) %>% arrange(GQ) # ditto

# Check top and bottom distribution of variants with 20-30 total reads
X %>% filter(DP >= 20 & DP <= 30) %>% .$GQ %>% hist() # depth 20-30 hist
X %>% filter(DP >= 20 & DP <= 30) %>% arrange(GQ) # depth 20-30, sort by GQ
X %>% filter(DP >= 20 & DP <= 30) %>% arrange(desc(GQ)) # depth 20-30, sort by GQ

Title <- sprintf("Arabidopsis TAI10 Subset")
SubTitle <- sprintf("Read Depth vs Genotype Call Quality Score (Trunc)")
X %>% filter(DP < 30) %>%
    ggplot(aes(x = DP %>% as.factor(), y = GQ, group = GT, colour = GT)) +
    geom_point(size = 0.3) +
    labs(title = Title, x = "DP", y = "GQ") +
    theme(axis.text.x = element_text(angle = -90, hjust = 0.5),
          plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5),
          legend.position = "right", legend.title=element_blank())
FileOut <- sprintf("%s%s", Title, ".jpeg")
ggsave(filename = FileOut, path = DirPlot, device = "jpeg")

# Entries with an alt allele, GQ < 200 (so plot is not squashed), calculate and plot alt read fraction
Title <- sprintf("Arabidopsis TAI10 Subset")
SubTitle <- sprintf("Alt Read Fraction vs Likelihood of Alt Allele")
X %>% mutate(AF = AO/(RO+AO)) %>% filter(GQ < 200) %>%
    ggplot(aes(x = AF, y = GQ, group = GT, colour = GT)) +
    geom_point(size = 0.3) +
    labs(title = Title, x = "Alt Allele Fraction", y = "Likelihood of Alt Allele") +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
          plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5),
          legend.position = "right", legend.title=element_blank())
FileOut <- sprintf("%s%s", Title, ".jpeg")
ggsave(filename = FileOut, path = DirPlot, device = "jpeg")

######################
## Set GQ Threshold ##
######################
GenomeData$GQ %>% summary()
#  Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
# 0.000    1.079    2.835   21.266    8.982 4217.400

# Entries with GQ less than ~ first quartile
X <- GenomeData %>% filter(GQ < 2.835)

# Keep only 0/0, 0/1 & 1/1 genotypes, calculate Alt allele frequency
X %<>% filter(GT %in% c("0/0", "0/1", "1/1")) %>% mutate(AF = AO/(RO+AO))

# Plot GQ vs DP, AO in separate panels
Title <- sprintf("Arabidopsis TAI10 Subset")
SubTitle <- sprintf("Genotype Quality vs Nr Reads grouped by Nr Alt Reads")
X %>% filter(AO <= 8) %>%
    ggplot(aes(x = DP, y = GQ, group = AO, colour = GT)) +
    geom_point(size = 0.3) +
    labs(title = Title, x = "Nr Reads", y = "Genotype Quality") +
    facet_wrap(~AO, ncol = 4) +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
          plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5),
          legend.position = "right", legend.title=element_blank())
FileOut <- sprintf("%s%s", Title, ".jpeg")
ggsave(filename = FileOut, path = DirPlot, device = "jpeg")
# That did not help!

# Plot GQ vs DP, AO in separate panels, for 0/1
Title <- sprintf("Arabidopsis TAI10 Subset")
SubTitle <- sprintf("Genotype Quality vs Nr Reads grouped by Nr Alt Reads for 0/1")
X %>% filter(AO <= 9, GT == "0/1") %>%
    ggplot(aes(x = DP, y = GQ, group = AO, colour = GT)) +
    geom_point(size = 0.3, colour = "blue") +
    labs(title = Title, x = "Nr Reads", y = "Genotype Quality") +
    facet_wrap(~AO, ncol = 4) +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
          plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5),
          legend.position = "right", legend.title=element_blank())
FileOut <- sprintf("%s%s", Title, ".jpeg")
ggsave(filename = FileOut, path = DirPlot, device = "jpeg")
# That did not help either!

# Test with set thresholds, look for doubtful calls above this
# If we still have a lot of calls with nr alt reads, may need to push higher

# Create two datasets, one below and one above threshold
X1 <- GenomeData %>% filter(GQ < 0.1)
X2 <- GenomeData %>% filter(GQ > 0.1 & GQ < 2.835)

# Sizes of datasets
X1 %>% dim() # 953  13
X2 %>% dim() # 16675    13

# Genotypes
X1$GT %>% fsummary()
# 0/0 0/1 1/1
# 494 451   8

X2$GT %>% fsummary()
# 0/0  0/1  1/1
# 6995 7678 2002

# Genotypes as percentage
((X1$GT %>% fsummary()) / length(X1$GT) * 100) %>% as.integer()
# 51 47  0

((X2$GT %>% fsummary()) / length(X2$GT) * 100) %>% as.integer()
# 41 46 12

# Nr alt reads
X1$AO %>% summary()
#  Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 2.000   2.000   2.000   2.598   2.000  74.000

X2$AO %>% summary()
#  Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 1.000   2.000   2.000   3.331   3.000 556.000

# The same median value - may not be an appropriate threshold





