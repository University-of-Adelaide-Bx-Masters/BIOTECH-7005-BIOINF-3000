##############
## Packages ##
##############
library(magrittr)
library(vcfR)
# library(ape)
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
```{r results='hide', eval=FALSE}
DirIn <- "/home/student/Practical_8"
DirOut <- "/home/student/Practical_8/Out"
if(! dir.exists(DirOut)) {dir.create(DirOut)}

DirPlot <- file.path(DirOut, "Plots")
if(! dir.exists(DirPlot)) {dir.create(DirPlot)}
```
#########
## IGV ##
#########
# Load ref and then Bam and VCF files in to IGV

# Look at the following locations
# 1:15,437,129-15,437,259 # MNPs
# 1:15,488,098-15,488,162 # Weak Indel
# 1:15,522,627-15,522,757 # Uncalled variants
# 1:16,509,693-16,526,713 # Highly variable region
# 1:17,074,700-17,210,875 # Highly variable region
# 1:20,054,605-20,190,780 # Non-unique read region


###################
## Load GFF File ##
###################
# DO NOT EXECUTE THIS SECTION
FileIn <- "Arabidopsis_thaliana.TAIR10.48.gff3.gz"

GFF <- read.table(file.path(DirIn, FileIn), sep="\t", quote="") #%>% as_tibble()

GFF %>% dim() # 791564      9

#########################
## Load Reference File ##
#########################
# DO NOT EXECUTE THIS SECTION
FileIn <- "Arabidopsis_thaliana.TAIR10.dna.toplevel.fa"

Ref <- ape::read.dna(file.path(DirIn, FileIn), format = "fasta")

###############
## Load VCF ##
###############
# Load file
FileIn <- "SRR5882792_Athaliana_TAIR10_picard.vcf"
VCF <- read.vcfR(file.path(DirIn, FileIn), verbose = TRUE)

VCF@meta %>% as_tibble() # metadata
VCF@fix %>% as_tibble() # co-ordinates etc and INFO ie data about this position in the genome
VCF@gt %>% as_tibble() # FORMAT ie genotypes

INFO2df(VCF[1:20,]) %>% as_tibble() # converts INFO key value pairs to df

#######################
## Analyse Genotypes ##
#######################
# Genotype info stored in separate columns for each sample - only one sample reported here
# Keys are stored in FORMAT col

# Convert FORMAT to DF, add co-ords etc
VCF@gt[,1] %>% unique() # all keys are present in all records, usually not the case

ColNames <- VCF@gt[1,1] %>% strsplit(":") %>% unlist() # get col names
GenomeData <- VCF@gt %>% as_tibble %>% dplyr::select(-FORMAT) %>% separate(unknown, ColNames, sep = ":") # convert strings to cols
GenomeData$Pos %<>% as.integer()
GenomeData$Qual %<>% as.numeric()
GenomeData$DP %<>% as.integer()
GenomeData$RO %<>% as.integer()
GenomeData$QR %<>% as.integer()
GenomeData$AO %<>% as.integer()
GenomeData$QA %<>% as.integer()

X <- VCF@fix %>% as_tibble() %>% dplyr::select(CHROM, POS, REF, ALT, QUAL) # get required INFO & other cols
colnames(X) <- c("Chr", "Pos", "Ref", "Alt", "Qual") # rename

GenomeData <- bind_cols(X, GenomeData) # merge INFO and FORMAT data
GenomeData$Qual %<>% as.numeric()

rm(X)

# Extract P(Best) vs P(NextBest)
GenomeData %<>% rowwise() %>%
    mutate(GLB = GL %>% strsplit(",") %>% unlist()  %>% as.numeric() %>% sort() %>% .[2]) %>%
    ungroup() # split, sort, take the 2nd value (first always 0)

# View col definitions and data to see what each col means
VCF@meta %>% as_tibble() %>% filter(grepl("FORMAT", value))
GenomeData

###########
## Plots ##
###########
## 1 ##
# Relationship between call and QA
Title <- sprintf("Arabidopsis TAI10 Alternate Allele Support by Genotype Call")
GenomeData %>% ggplot(aes(x= GT, y = QA, group = GT, colour = GT)) +
    geom_boxplot() +
    labs(title = Title, x = "GT", y = "QA") +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
          plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5),
          legend.position = "none", legend.title=element_blank())
FileOut <- sprintf("%s%s", Title, ".jpeg")
ggsave(filename = FileOut, path = DirPlot, device = "jpeg")

# Positions called hom ref with strong evidence for alt allele
GenomeData %>% filter(GT == "0/0", QA > 5000)
# Look at 1:15082284

# Positions called 0/2 etc
GenomeData %>% filter(GT %nin% c("0/0", "0/1", "1/1")) #%>% View()
# Look at 1:7822680, 1:8675530, 1:12033496

# Relationship between call and QA truncated
Title <- sprintf("Arabidopsis TAI10 Alternate Allele Support by Genotype Call (Trunc)")
GenomeData %>% filter(QA < 400) %>% ggplot(aes(x= GT, y = QA, group = GT, colour = GT)) +
    geom_boxplot() +
    labs(title = Title, x = "GT", y = "QA") +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
          plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5),
          legend.position = "none", legend.title=element_blank())
FileOut <- sprintf("%s%s", Title, ".jpeg")
ggsave(filename = FileOut, path = DirPlot, device = "jpeg")

## 2 ##
# Relationship between call and GLB
Title <- sprintf("Arabidopsis TAI10 Genotype Call Quality Score by Genotype Box Plot")
GenomeData %>% ggplot(aes(x= GT, y = GLB, group = GT, colour = GT)) +
    geom_boxplot() +
    labs(title = Title, x = "Genotype", y = "GLB") +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
          plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5),
          legend.position = "none", legend.title=element_blank())
FileOut <- sprintf("%s%s", Title, ".jpeg")
ggsave(filename = FileOut, path = DirPlot, device = "jpeg")

Title <- sprintf("Arabidopsis TAI10 Genotype Call Quality Score by Genotype Histogram")
GenomeData %>% filter(GT %in% c("0/0", "0/1", "1/1")) %>%
    ggplot(aes(x = GLB, group = GT, colour = GT)) +
    geom_histogram(colour = "grey", fill = "lightblue") +
    labs(title = Title, x = "GLB", y = "Count") +
    facet_wrap(~GT) +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
          plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5),
          legend.position = "none", legend.title=element_blank())
FileOut <- sprintf("%s%s", Title, ".jpeg")
ggsave(filename = FileOut, path = DirPlot, device = "jpeg")

Title <- sprintf("Arabidopsis TAI10 Reference Score vs Genotype Call Quality Score")
GenomeData %>% filter(GT %in% c("0/0", "0/1", "1/1")) %>%
    ggplot(aes(x = GLB, y = Qual, group = GT, colour = GT)) +
    geom_point(size = 0.3) +
    labs(title = Title, x = "GLB", y = "Qual") +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
          plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5),
          legend.position = "right")
FileOut <- sprintf("%s%s", Title, ".jpeg")
ggsave(filename = FileOut, path = DirPlot, device = "jpeg")

Title <- sprintf("Arabidopsis TAI10 Reference Score vs Genotype Call Quality Score (Trunc)")
GenomeData %>% filter(GT %in% c("0/0", "0/1", "1/1")) %>% filter(GLB > -100, Qual < 1000) %>%
    ggplot(aes(x = GLB, y = Qual, group = GT, colour = GT)) +
    geom_point(size = 0.3) +
    labs(title = Title, x = "GLB", y = "Qual") +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
          plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5),
          legend.position = "right", legend.title=element_blank())
FileOut <- sprintf("%s%s", Title, ".jpeg")
ggsave(filename = FileOut, path = DirPlot, device = "jpeg")


# Subset variants with 0/1 & 1/1 genotypes, with only 1 alt allele
X <- GenomeData %>% filter(GT %in% c("0/1", "1/1")) %>% filter(str_count(AD, ",") == 1) # positions with alt, only 2 alleles

## How many reads support alt genotype calls ##
X$AO %>% hist()
X$AO[X$AO < 100] %>% hist()
X$AO[X$AO < 10] %>% hist()

X$AO %>% summary()

X$AO[X$AO < 10] %>% fsummary()
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 1.00    2.00    4.00   19.91    8.00 5309.00

X %>% filter(AO <= 1) %>% View()
# Look at 3:12247004 - why is the variant shifted?

X %>% filter(AO == 2) %>% arrange(desc(GLB)) # sort by GLB worst to best
X %>% filter(AO == 2) %>% arrange(GLB) # sort by GLB best to worst

X %>% filter(DP >= 20 & DP <= 30) %>% .$GLB %>% hist() # depth 20-30 hist
X %>% filter(DP >= 20 & DP <= 30) %>% arrange(desc(GLB)) # depth 20-30, sort by GLB worst to best
X %>% filter(DP >= 20 & DP <= 30) %>% arrange(GLB) # depth 20-30, sort by GLB best to worst

Title <- sprintf("Arabidopsis TAI10 Subset Read Depth vs Genotype Call Quality Score (Trunc)")
X %>% filter(DP < 30) %>%
    ggplot(aes(x = DP %>% as.factor(), y = GLB, group = GT, colour = GT)) +
    geom_point(size = 0.3) +
    labs(title = Title, x = "DP", y = "GLB") +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
          plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5),
          legend.position = "right", legend.title=element_blank())
FileOut <- sprintf("%s%s", Title, ".jpeg")
ggsave(filename = FileOut, path = DirPlot, device = "jpeg")

Title <- sprintf("Arabidopsis TAI10 Subset Alt Read Depth vs Genotype Call Quality Score (Trunc)")
X %>% filter(AO < 30) %>%
    ggplot(aes(x = AO %>% as.factor(), y = GLB, group = GT, colour = GT)) +
    geom_point(size = 0.3) +
    labs(title = Title, x = "AO", y = "GLB") +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
          plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5),
          legend.position = "right", legend.title=element_blank())
FileOut <- sprintf("%s%s", Title, ".jpeg")
ggsave(filename = FileOut, path = DirPlot, device = "jpeg")






