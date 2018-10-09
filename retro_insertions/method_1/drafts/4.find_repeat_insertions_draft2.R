### 4.find_repeat_insertions.R ###

# This script takes discordant reads from svaba output overlaps with a 
# custom DE repeats annotation to find retrotransposon insertions, 
# and filters those lying within breakpoints from CNVnator output within
# following distances: 0, 2, 5, 10, 100, 1000 kb

# module load R/3.5.1 for NCI

library(rtracklayer)
library(GenomicRanges)
library("BiocParallel")
library(org.Hs.eg.db)
library(GEOquery)
library("BSgenome.Hsapiens.UCSC.hg19")


### 0. Define variables/paths ###

project <- "hgsoc_repeats"

#home_dir="/g/data1a/ku3/jt3341/"
home_dir <- "/Users/jamestorpy/nciShare/"
project_dir=paste0(home_dir, "/projects/", project, "/genome/")
results_dir=paste0(project_dir, "/results/")
Robject_dir=paste0(project_dir, "/Robjects/")
system(paste0("mkdir -p ", results_dir))
system(paste0("mkdir -p ", Robject_dir))

ref_dir="/Users/jamestorpy/nciShare/genomes/"
#ref_dir="/g/data1a/ku3/jt3341/genomes/"
svaba_dir=paste0(results_dir, "/svaba/")
cnvnator_dir=paste0(results_dir, "/cnvnator/")
out_dir=paste0(results_dir, "/tables/")
system(paste0("mkdir -p ", out_dir))

#sample_id <- "AOCS-153-1-2_subset"

# define distances from breakpoints in which to look for repeat insertions:
dist <- list(1000, 2000, 5000, 10000, 100000, 1000000)


##############################################################################
### 1. Load in svaba discrepencies and repeat annotation ###
##############################################################################

if ( file.exists(paste0(Robject_dir, "/rep_hg19_gr.rds")) ) {
  
  rep_gr <- readRDS(paste0(Robject_dir, "/rep_hg19_gr.rds"))

} else {
  
  rep_annot <- read.table(paste0(ref_dir, "/repeats/repeats.hg19.gff"))
  
  # convert repeat annotation to GRanges object:
  rep_gr <- GRanges(
    seqnames = Rle(rep_annot$V1),
    ranges = IRanges(start=as.integer(rep_annot$V4), 
                     end = as.integer(rep_annot$V5)),
    strand = Rle("*"),
    gene_id = gsub("ID=", "", rep_annot$V9)
  )
  # remove unwanted chromosomes:
  rep_gr <- rep_gr[grep("M|G|K", seqnames(rep_gr), 
                        invert = T)]
  # remove dust, trf, rRNA:
  rep_gr <- rep_gr[!(rep_gr$gene_id %in% c("dust", "trf", "rRNA"))]
  # save rep_gr as RDS file:
  saveRDS(rep_gr, paste0(Robject_dir, "/rep_hg19_gr.rds"))
  
}

if ( file.exists(paste0(Robject_dir, "/discord_list_hg19.rds")) ) {
  
  discord <- readRDS(paste0(Robject_dir, "/discord_list_hg19.rds"))

} else {
  
  # fetch all svaba discordant reads files from svaba_dir:
  discord_files <- split(
    grep(
      "subset", list.files(svaba_dir, pattern = "discordant", 
        full.names = T), value = T, invert = T
    ), 1:length(discord_files)
  )
  sample_names <- unlist(lapply(discord_files, function(x) {
    return(gsub(".discordant.*$", "", basename(x)))
  }))
  
  discord_list <- lapply(discord_files, function(x) {
    sample_id <- gsub(".discordant.*$", "", basename(x))
    # gunzip svaba discordant reads if needed and load table:
    if ( file.exists(x) ) {
      gunzip(x)
    }
    
    discord <- read.table(paste0(svaba_dir, "/", sample_id, ".discordant.txt"), 
                          sep = "\t", header = T, fill = T)
    
    # convert discord into GRanges object with both donor and acceptor data as main ranges:
    both_seqnames <- c(discord$chr1, discord$chr2)
    both_start <- c(discord$pos1, discord$pos2)
    both_actual_strand <- c(as.character(discord$strand1), as.character(discord$strand2))
    both_match_seqnames <- c(discord$chr2, discord$chr1)
    both_match_pos <- c(discord$pos2, discord$pos1)
    both_match_strand <- c(as.character(discord$strand2), as.character(discord$strand1))
    
    discord_gr <- GRanges(
      seqnames = Rle(paste0("chr",both_seqnames)),
      ranges = IRanges(start=as.integer(both_start), 
                       width=1),
      strand = Rle("*"),
      actual_strand = Rle(both_actual_strand),
      match_seqnames <- Rle(paste0("chr", both_match_seqnames)),
      match_pos <- as.integer(both_match_pos),
      match_strand <- both_match_strand
    )
    colnames(values(discord_gr)) <- c("actual_strand", "match_seqnames", 
                                      "match_pos", "match_strand")
    # remove unwanted chromosomes:
    discord_gr <- discord_gr[grep("M|G|K", as.character(seqnames(discord_gr)), 
                                  invert = T)]
    return(discord_gr)
    
  })
}
names(discord_list) <- sample_names

discord_list <- lapply(discord_list, function(x) {
  x$supporting_reads <- NA
  uniq <- unique(x)
  for ( j in 1:length(uniq) ) {
    print(j)
    uniq$supporting_reads[j] <- length(x[start(ranges(x))==start(ranges(uniq[j]))])
  }
  return(uniq)
})
saveRDS(discord_list, paste0(Robject_dir, "/discord_list_hg19.rds"))

verified_discords <- lapply(discord_list, function(x) return(x[x$supporting_reads > 1]))
saveRDS(verified_discords, paste0(Robject_dir, "/verified_discords_hg19.rds"))


##############################################################################
### 2. Find retrotransposon insertions ###
##############################################################################

if ( file.exists(paste0(Robject_dir, "/rep_insertions_hg19.rds"))) {
  
  rep_insertions <- readRDS(paste0(Robject_dir, "/rep_insertions_hg19.rds"))
  
} else {
  
  rep_insertions <- lapply(discord_list, function(x) {
    # find overlaps of discordant reads with repeat annotation:
    olaps <- findOverlaps(x, rep_gr)
    
    incl_discord <- x[queryHits(olaps)]
    incl_rep <- rep_gr[subjectHits(olaps)]
    
    rep_insert_gr <- GRanges(
      seqnames = seqnames(incl_rep),
      ranges = ranges(incl_rep),
      strand = strand(incl_rep),
      gene_id = incl_rep$gene_id,
      fusion_half_1_seqnames = as.character(seqnames(incl_discord)),
      fusion_half_1_pos = as.integer(start(ranges(incl_discord))),
      fusion_half_1_strand = incl_discord$actual_strand,
      fusion_half_2_seqnames = incl_discord$match_seqnames,
      fusion_half_2_pos = incl_discord$match_pos,
      fusion_half_2_strand = incl_discord$match_strand
    )
    return(rep_insert_gr)
  })
  
  saveRDS(rep_insertions, paste0(Robject_dir, "/rep_insertions_hg19.rds"))
  
}

rep_inserts_pooled <- do.call(getMethod(c, "GenomicRanges"), rep_inserts_pooled)

saveRDS(rep_inserts_pooled, paste0(Robject_dir, "/repeat_insertions_hg19_annotated.rds"))

# make rep insertions gtf with co-ordinate of where repeat inserted:
rep_inserts_pooled_gtf <- GRanges(
  seqnames = Rle(rep_inserts_pooled$fusion_half_2_seqnames),
  ranges = IRanges(start=rep_inserts_pooled$fusion_half_2_pos, width=1),
  strand = "*",
  source = "svaba",
  feature = "repeat",
  score = ".",
  frame = ".",
  repeat_id = rep_inserts_pooled$gene_id
)
rep_inserts_pooled_gtf <- unique(rep_inserts_pooled_gtf)
export(rep_inserts_pooled_gtf, paste0(out_dir, "/repeat_inserts_pooled_hg19.gtf"))


##############################################################################
### 4. Load in CNVnator results and filter out CNVnator breakpoints not 
# supported by at least 1 svaba discordant read ###
##############################################################################

CNV_files <- split(
  grep(
    "subset", list.files(cnvnator_dir, pattern = "CNV_calls", 
                         full.names = T), value = T, invert = T
  ), 1:5
)
CNV_files <- split(
  CNV_files[gsub("_CNV_calls.txt", "", basename(CNV_files)) %in% samplenames], 1:length(CNV_files)
)

# read in CNVnator results:
CNV_list <- lapply(CNV_files, function(x) {
  CNV_breakpoints <- read.table(x, header = F, sep = "\t")
  
  # format data frame:
  CNV_df <- CNV_breakpoints[grep("M|G|K", CNV_breakpoints$V2, invert = T),]
  CNV_df$type <- CNV_df$V1
  split_df <- do.call("rbind", strsplit(as.character(CNV_df$V2), ":"))
  CNV_df$seqnames <- split_df[,1]
  CNV_df$start <- do.call("rbind", strsplit(as.character(split_df[,2]), "-"))[,1]
  CNV_df$end <- do.call("rbind", strsplit(as.character(split_df[,2]), "-"))[,2]
  CNV_df <- subset(CNV_df, select = c("type", "seqnames", "start", "end"))
  
  # convert to GRanges object with all breakpoints as both start and end ranges:
  CNV_gr <- GRanges(
    seqnames = Rle(c(CNV_df$seqnames, CNV_df$seqnames)),
    ranges = IRanges(start=c(as.integer(CNV_df$start), as.integer(CNV_df$end)), 
                     end=c(as.integer(CNV_df$start), as.integer(CNV_df$end))),
    strand = Rle("*"),
    type = c(as.character(CNV_df$type), as.character(CNV_df$type)),
    actual_start = c(as.integer(CNV_df$start), as.integer(CNV_df$start)),
    actual_end = c(as.integer(CNV_df$end), as.integer(CNV_df$end))
  )
  return(CNV_gr)
})

i=1
verified_breakpoints <- lapply(CNV_list, function(x) {
  bp_olaps <- findOverlaps(x, discord_list[[i]])
  i <<- i+1
  return(unique(x[queryHits(bp_olaps)]))
})
pooled_verified_breakpoints <- do.call(getMethod(c, "GenomicRanges"), verified_breakpoints)
  

# create breakpoint gtf file:
pooled_breakpoint_gtf <- GRanges(
  seqnames = seqnames(pooled_verified_breakpoints),
  ranges = ranges(pooled_verified_breakpoints),
  strand = "*",
  source = "cnvnator",
  feature = "breakpoint",
  score = ".",
  frame = ".",
  type = pooled_verified_breakpoints$type
)
export(pooled_breakpoint_gtf, paste0(out_dir, "/breakpoints_pooled.gtf"))


##############################################################################
### 5. Find odds of any repeat insertions at varying distances from verified 
# breakpoints ###
##############################################################################


# convert to list with first element deletions, second duplications:
deletion_gr <- pooled_verified_breakpoints[pooled_verified_breakpoints$type == "deletion",]
duplication_gr <- pooled_verified_breakpoints[pooled_verified_breakpoints$type == "duplication",]

# create function to expand genomic ranges either direction from start/end
# and optionally remove original ranges:
n=1
exp_gr <- function(Length, gr, orig) {
  
  print(paste0("Processing range ", n))
  
  # define lengths of chromosomes:
  seq_lengths <- seqlengths(Hsapiens)[!grepl("_",
                                             names(seqlengths(Hsapiens)))]
  # remove unwanted chromosome constructs:
  gr <- gr[grep("[0-9].[0-9]|MT|G|K", seqnames(gr), 
                      invert = T)]
  
  if ( length(gr) > 0 ) {
    # reduce ranges to ensure no overlaps/double annotations:
    gr <- unique(reduce(gr))

    # assign length of all chromosomes to the variable name 'seq_lengths':
    gr$seq_lengths <- rep(NA, length(gr))
    
    for ( v in names(seq_lengths) ) {
      gr$seq_lengths[as.character(seqnames(gr)) == v] <- seq_lengths[v]
    }
    
    # add length upstream of each range to another gr:
    start_gr <- gr
    # make end position equal to original start position:
    end(ranges(start_gr)) <- start(ranges(start_gr))
    # make start position the original start position - Length if the original start is at least Length
    # away from start of chromosome:
    start(start_gr)[start(ranges(start_gr)) >= Length] <- 
      start(ranges(start_gr))[start(ranges(start_gr)) >= Length] - Length
    
    # add length downstream of each range to another gr:
    end_gr <- gr
    
    # make start position equal to original end position:
    start(ranges(end_gr)) <- end(ranges(end_gr))
    # make end position the original end position + Length if the end is at least Length
    # away from end of chromosome:
    far_enuff <- end(ranges(end_gr)) <= (end_gr$seq_lengths - Length)
    end(ranges(end_gr))[far_enuff] <- end(ranges(end_gr))[far_enuff] + Length
    
    if (orig) {
      # include original sequences in GRanges object:
      start(ranges(gr)) <- start(ranges(start_gr))
      end(ranges(gr)) <- end(ranges(end_gr))
      # remove ranges with end=1:
      gr <- gr[!(end(ranges(gr))==1)]
      n <<- n+1
      return(gr)
    }
    
    # remove co-ordinates that fall after the end of seq_lengths values:
    end_gr <- end_gr[!(end(ranges(end_gr)) > end_gr$seq_lengths)]
    
    # make end position the end of the chromosome if the original end is not at least Length
    # away from end of chromosome:
    not_far_enuff <- end(ranges(end_gr)) > (end_gr$seq_lengths - Length)
    if ( length(not_far_enuff) > 0 ) {
      end(ranges(end_gr))[not_far_enuff] <- end_gr$seq_lengths[not_far_enuff]
    }
    gr <- c(start_gr, end_gr)
    # remove ranges with end=1:
    gr <- gr[!(end(ranges(gr))==1)]
    # return combined start and end annotations:
    n <<- n+1
    return(gr)
  }
}

# Create list of deletions with ranges extending different distances outwards 
# from segments, without original (deleted) segment:
n=1
deletion_borders <- lapply(dist, exp_gr, gr=deletion_gr, orig=FALSE)
non_deletion <- lapply(deletion_borders, gaps)

# Create list of duplications with ranges extending different distances outwards 
# from segments, keeping original segment:
n=1
duplication_borders <- lapply(dist, exp_gr, gr=duplication_gr, orig=FALSE)
non_duplication <- lapply(duplication_borders, gaps)


# create function to calculate odds of group genes are enriched for 
# marks compared with all protein-coding genes
oddsRatio <- function(nonregion, region){
  
  # determine how many repeats in region:
  r_olaps <- findOverlaps(rep_inserts_pooled_gtf, region)
  region_no <- length(unique(queryHits(r_olaps)))
  # determine length of region:
  region_length <- sum(width(reduce(region)))
  
  # determine how many repeats in non-region:
  nr_olaps <- findOverlaps(rep_inserts_pooled_gtf, non_region)
  non_region_no <- length(unique(queryHits(nr_olaps)))
  # determine length of non-region:
  non_region_length <- sum(width(reduce(non_region)))
  
  # calculate odds of region containing a repeat insertion:
  first_ratio <- region_no/region_length
  
  # calculate odds of non-region containing a repeat insertion:
  second_ratio <- non_region_no/non_region_length
  
  # calculate odds ratio:
  odds <- first_ratio/second_ratio
  
  # calculate standard error for odds:
  se <- sqrt(1/region_length + 1/region_no + 1/non_region_length + 
               1/non_region_no)
  print("Odds are:")
  cat(odds)
  cat("\n")
  print("Std error is:")
  cat(se)
  cat("\n")
  
  # calculate p-value for odds:
  database_enrichment <-
    matrix( c(as.numeric(region_no), as.numeric(non_region_no), 
              as.numeric(region_length), as.numeric(non_region_length)), nrow = 2,
            ncol = 2 )
  colnames(database_enrichment) = c("insertions", "length")
  rownames(database_enrichment) = c("region", "non_region")
  
  pval <- chisq.test(database_enrichment)$p.value
  print("P-value is:")
  cat(pval)
  cat("\n")
  
  result <- data.frame(odds, se, pval)
  colnames(result) <- c("odds", "std_error", "pval")
  
  return(result)
}

bp_repeat_enrichment_odds <- lapply()


##############################################################################
### 6. Find odds of specific repeat insertions at varying distances from 
# verified breakpoints ###
##############################################################################

