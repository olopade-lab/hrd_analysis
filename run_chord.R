library(CHORD)
# library(BSgenome.Hsapiens.UCSC.hg38)
library(BSgenome.Hsapiens.UCSC.hg19)
library(mutSigExtractor)
library(optparse)

# NEED TO: check that genomes are both available

options(stringsAsFactors=F)

option_list <- list(
  make_option(c("-v", "--verbose"), action="store_true", 
                     default=FALSE, help="Print extra output"),
  make_option(c("-w", "--overwrite"), action="store_true", 
              default=FALSE, help="Overwrite output files"),
  make_option(c("-s", "--sv"), type="character",  
                     help="SV VCF path pattern",
                     dest="sv_path"),
  make_option(c("-n", "--snv"), type="character",  
                     help="SNV VCF path pattern",
                     dest="snv_path"),
  make_option(c("-i", "--indel"), type="character",  
                     help="Indel VCF path pattern",
                     dest="indel_path"),
  make_option(c("-o", "--output"), type="character",  
                     help="Output directory path")
  # make_option(c("-g", "--genome"), type="character",  
  #                    help="Reference genome (hg19 or hg38)",
  #                    default="hg19")
)


opt <- parse_args(OptionParser(option_list=option_list))

# if (!(opt$genome %in% c("hg19", "hg38"))) {
#   stop("Please specify correct reference genome (hg19 or hg38)")
# }

print(opt)

if (is.null(opt$snv_path)) {
  stop("Please provide SNV VCF path")
}

if (is.null(opt$sv_path)) {
  stop("Please provide SV VCF path")
}

if (is.null(opt$indel_path)) {
  stop("Please provide indel VCF path")
}

if (!dir.exists(opt$output)) {
  stop("Output directory does not exist")
}

if (!(substr(opt$output, nchar(opt$output)-1, nchar(opt$output)) == "/")) {
  opt$output <- paste0(opt$output, "/")
}

contexts_dir <- paste0(opt$output, "contexts/")

if (!dir.exists(contexts_dir)) {
  dir.create(contexts_dir)
}

setwd("/")

snv_vcf_files <- data.frame(
  snv=list.files(dirname(opt$snv_path), pattern=basename(opt$snv_path), full.names=TRUE)
)

print(snv_vcf_files)
snv_vcf_files$sample <- sapply(strsplit(basename(snv_vcf_files$snv),'.', fixed = TRUE),`[`,1)

indel_vcf_files <- data.frame(
  indel=list.files(dirname(opt$indel_path), pattern=basename(opt$indel_path), full.names=TRUE)
)
indel_vcf_files$sample <- sapply(strsplit(basename(indel_vcf_files$indel),'.', fixed = TRUE),`[`,1)

sv_vcf_files <- data.frame(
  sv=list.files(dirname(opt$sv_path), pattern=basename(opt$sv_path), full.names=TRUE)
)
sv_vcf_files$sample <- sapply(strsplit(basename(sv_vcf_files$sv),'.', fixed = TRUE),`[`,1)

vcf_files <- merge(snv_vcf_files, indel_vcf_files, all=FALSE)
vcf_files <- vcf_files[!duplicated(vcf_files),]
vcf_files <- merge(vcf_files, sv_vcf_files, all=FALSE)
vcf_files <- vcf_files[!duplicated(vcf_files),]

print(vcf_files)

if (nrow(vcf_files) == 0) {
  stop("No vcf files match your paths for all three mutation types.")
}

if (opt$verbose) {
  print("List of samples run:")
  print(vcf_files)
}

extract_vcf <- function(sample) {
  # Extract SNV
  df_snv <- readVcfFields(sample$snv, fields=c('CHROM','POS','REF','ALT'))
  colnames(df_snv) <- c('chrom','pos','ref','alt')
  df_snv <- df_snv[which(df_snv$chrom %in% c(1:23, "X", "Y")), ]
  
  # Extract indel
  df_indel <- readVcfFields(sample$indel, fields=c('CHROM','POS','REF','ALT'))
  colnames(df_indel) <- c('chrom','pos','ref','alt')
  df_indel <- df_indel[which(df_indel$chrom %in% c(1:23, "X", "Y")), ]
  
  # Extract SV
  vcf_sv <- readVcfFields(sample$sv, fields='INFO')
  vcf_sv_info <- strsplit(vcf_sv$INFO,';')
  sv_type <- sapply(vcf_sv_info,`[`,4)
  sv_type <- gsub('SVTYPE=','',sv_type)
  sv_len <- sapply(vcf_sv_info,`[`,3)
  sv_len <- gsub('SVLEN=','',sv_len)
  sv_len <- as.integer(sv_len)
  df_sv <- data.frame(sv_type, sv_len)
  df_sv$sv_len <- abs(df_sv$sv_len)
  
  return(list(snv = df_snv, indel = df_indel, sv = df_sv))
}

### TEST ###

# params <- as.list(vcf_files[1,])
# out_path <- paste0(contexts_dir, params$sample, '_contexts.txt')
# list_df <- extract_vcf(params)
# 
# contexts <- extractSigsChord(
#   df.snv = list_df$snv,
#   df.indel = list_df$indel,
#   df.sv = list_df$sv,
#   sample.name='N010911',
#   ref.genome=BSgenome.Hsapiens.UCSC.hg19
# )
# 
# print(contexts)

### END TEST ###

for(i in 1:nrow(vcf_files)) {
  
  params <- as.list(vcf_files[i,])
  out_path <- paste0(contexts_dir, params$sample, '_contexts.txt')
  list_df <- extract_vcf(params)
  print(head(list_df$snv))
  print(head(list_df$sv))
  print(head(list_df$indel))
  
  if(!file.exists(out_path)) {
    extractSigsChord(
      df.snv = list_df$snv,
      df.indel = list_df$indel,
      df.sv = list_df$sv,
      sample.name=params$sample,
      ref.genome=BSgenome.Hsapiens.UCSC.hg19,
      output.path=out_path, verbose=opt$verbose
    )
  } else {
    if (opt$overwrite) {
      file.remove(out_path)
      extractSigsChord(
        df.snv = list_df$snv,
        df.indel = list_df$indel,
        df.sv = list_df$sv,
        sample.name=params$sample,
        ref.genome=BSgenome.Hsapiens.UCSC.hg19,
        output.path=out_path, verbose=opt$verbose
      )
    }
  }
}

context_files <- list.files(contexts_dir, full.names = TRUE)

l_contexts <- lapply(context_files, function(i){
  read.delim(i, check.names = FALSE)
})

merged_contexts <- do.call(rbind, l_contexts)

write.table(merged_contexts, paste0(opt$output, 'merged_contexts.txt'), sep='\t', quote = FALSE)

chord_output <- chordPredict(merged_contexts, do.bootstrap = TRUE, verbose = opt$verbose)
write.table(chord_output, paste0(opt$output, 'chord_pred.txt'), sep='\t', quote = FALSE)