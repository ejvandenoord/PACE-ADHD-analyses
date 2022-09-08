#Performing QC and preparing covariates for the MWAS
#Code written by Mandy Meijer, Karolina Aberg and Edwin van den Oord, 2020. Center for Biomarker Research and Precision Medicine (BPM), Virginia Commonwealth University
#THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

#This script needs only to be run once
#Only in line 7 and 8 some changes from the user are required
#The expected outcome are a folder 'rw' with a coverage matrix, CpG locations, and chromosome name; a pdf file with a PCA plot of technical probes; outcome file with covariates

dirproject = "/nfs/home2/mmeijer/ramwas/"
base = "/nfs/home2/mmeijer/ramwas/" #Please insert own directory which contains iDAT files

file = "/nfs/home2/mmeijer/ramwas/PACE_ADHD_sample_sheet.csv"  #This is the sample sheet 

#### start

library(ramwas)
library(minfi)

targets = read.csv(file = file,
                    stringsAsFactors = FALSE)
rgSet = read.metharray.exp(
  base = base,
  targets = targets,
  extended = TRUE,
  verbose = TRUE)

work_dir   = paste0(dirproject,"/covariates")
dir.create(work_dir,showWarnings=F,recursive=T)
setwd(work_dir) 

### start function QC and create covariates
if( "IlluminaHumanMethylation450k" %in% rgSet@annotation ){
  host = "http://www.sickkids.ca/MS-Office-Files/Research/Weksberg%20Lab/"
  files = c(
    i450k_ns.xlsx = "48639-non-specific-probes-Illumina450k.xlsx",
    i450k_pl.xlsx = "48640-polymorphic-CpGs-Illumina450k.xlsx")
  
  for( i in seq_along(files) )
    download.file(
      url = paste0(host, files[i]),
      destfile = names(files)[i],
      mode = "wb",
      quiet = TRUE)
  
  library(readxl)
  ex1 = read_excel("i450k_ns.xlsx", sheet = 1)
  ex2 = read_excel("i450k_pl.xlsx", sheet = 1)
  ex3 = read_excel("i450k_pl.xlsx", sheet = 2)
  
  exclude.snp = unique(c(
    ex1$TargetID,
    ex2$PROBE,
    ex3$PROBE[ (ex3$BASE_FROM_SBE < 10) & (ex3$AF > 0.01)]))
  rm(host, files, i, ex1, ex2, ex3)
} else {
  host = "https://static-content.springer.com/esm/art%3A10.1186%2Fs13059-016-1066-1/MediaObjects/"
  files = c(
    S1_cross_reactive.csv     = '13059_2016_1066_MOESM1_ESM.csv',
    S4_snp_cpg.csv            = '13059_2016_1066_MOESM4_ESM.csv',
    S5_snp_base_extension.csv = '13059_2016_1066_MOESM5_ESM.csv',
    S6_snp_body.csv           = '13059_2016_1066_MOESM6_ESM.csv')
  
  for( i in seq_along(files) )
    download.file(
      url = paste0(host, files[i]),
      destfile = names(files)[i],
      mode = "wb",
      quiet = TRUE)
  
  snpcpgs1 = read.csv('S1_cross_reactive.csv', stringsAsFactors = FALSE)
  snpcpgs4 = read.csv('S4_snp_cpg.csv', stringsAsFactors = FALSE)
  snpcpgs5 = read.csv('S5_snp_base_extension.csv', stringsAsFactors = FALSE)
  snpcpgs6 = read.csv('S6_snp_body.csv', stringsAsFactors = FALSE)
  
  exclude.snp = unique(c(
    snpcpgs1$X,
    snpcpgs4$PROBE,
    snpcpgs5$PROBE,
    snpcpgs6$PROBE[
      pmax(snpcpgs6$VARIANT_START - snpcpgs6$MAPINFO, 
           snpcpgs6$MAPINFO - snpcpgs6$VARIANT_END) < 10]))
  rm(host, files, i, snpcpgs1, snpcpgs4, snpcpgs5, snpcpgs6)
}

#Identify probes with low bead count
lb = getNBeads(rgSet) < 3
pi1 = getProbeInfo(rgSet, type = "I")
pi2 = getProbeInfo(rgSet, type = "II")
ex1 = pi1$Name[rowMeans(lb[pi1$AddressA,] | lb[pi1$AddressB,]) > 0.01]
ex2 = pi2$Name[rowMeans(lb[pi2$AddressA,]) > 0.01]
exclude.bds = unique(c(ex1, ex2))
rm(lb, pi1, pi2, ex1, ex2)

#Identify probes and samples with low detection p-values
hp = detectionP(rgSet) > 0.01
exclude.hpv = rownames(hp)[rowMeans(hp) > 0.01]
keep.samples = colMeans(hp) < 0.01
rm(hp)

#Exclusion of low quality samples and probes
rgSet = subsetByLoci(
  rgSet = rgSet[,keep.samples],
  excludeLoci = c(exclude.snp, exclude.bds, exclude.hpv))

#Obtain methylation estimates and save in RaMWAS format
rgSetRaw = fixMethOutliers(preprocessRaw(rgSet))
beta = getBeta(rgSetRaw)

work_dir   = dirproject
setwd(work_dir) 

dir.create('rw', showWarnings = FALSE)

rng = granges(mapToGenome(rgSet))
chr = seqnames(rng)

# Save CpG locations
locs = cbind(chr = as.integer(chr), position = start(rng))
fmloc = fm.create.from.matrix(
  filenamebase = paste0("rw/CpG_locations"),
  mat = locs,
  size = 4)
close(fmloc)
writeLines(con = 'rw/CpG_chromosome_names.txt', text = levels(chr))

# Save estimates
fm = fm.create.from.matrix(
  filenamebase = paste0("rw/Coverage"),
  mat = t(beta))
close(fm)

work_dir   = dirproject
setwd(work_dir) 
work_dir   = paste0(dirproject,"/covariates")
setwd(work_dir) 


### create covariates
# cell types
covariates.cel = estimateCellCounts(
  rgSet = rgSet, 
  compositeCellType = "Blood",
  cellTypes = c("CD8T", "CD4T", "NK", 
                "Bcell", "Mono", "Gran"),
  meanPlot = FALSE)

# Median methylated and unmethylated intensity
covariates.umm = getQC(rgSetRaw)

# Principal components analysis (PCA) on control probes
controlType = unique(getManifest(rgSet)@data$TypeControl$Type)
controlSet = getControlAddress(rgSet, controlType = controlType)
probeData = rbind(getRed(rgSet)[controlSet,], getGreen(rgSet)[controlSet,])

data = probeData - rowMeans(probeData)
covmat = crossprod(data)
eig = eigen(covmat)

pdf("screeplot_techPC_values.pdf") #plot PCA control probes
plotPCvalues(eig$values, n = 20)
plotPCvectors(eig$vectors[,1], i = 1, col = 'blue')
dev.off()

saveRDS(eig$values,paste0(work_dir,"//PCplot_tech.rds"))
        
nPCs = 20
covariates.pca = eig$vectors[,seq_len(nPCs)]
colnames(covariates.pca) = paste0('techPC',seq_len(nPCs))
rm(probeData, data, covmat, eig, nPCs) 

# Phenotypic covariates from the sample sheet
rownames(targets) = targets$Basename;
targets$xSlide = paste0('x',targets$Slide) # Force Slide to be categorical
covariates.phe = targets[
  rownames(covariates.umm),
  colnames(targets)]

covariates = data.frame(
  samples = rownames(covariates.umm),
  covariates.phe,
  covariates.umm,
  covariates.pca,
  covariates.cel
  )

# save covariates to file
write.csv(covariates,"covariates_ramwas_PACE_ADHD.csv",  row.names = FALSE, col.names = TRUE)







