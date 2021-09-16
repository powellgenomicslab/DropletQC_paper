# Load R libraries
library(rtracklayer)
library(GenomicFeatures)

# Avoid scientific notation
options(scipen=999)

# Set base directory where files should be created
base_dir <- "~/Downloads/"

# Define function to simulate an annotation
simulate_annotation <- function(num_genes=100,
                                chrom_name="chr1",
                                chrom_length=1000000,
                                exon_width = 350,
                                gff_path = paste0(base_dir, "test.gff3"),
                                input_seed = 41){

  # Create exon ranges to sample from
  names(chrom_length) <- chrom_name
  potential_exons <- unlist(tileGenome(seqlengths = chrom_length, tilewidth = exon_width))
  strand(potential_exons) <- "+"

  # Select num_genes*2 exons
  set.seed(input_seed)
  selected_exons <- sample(potential_exons, (num_genes*2), replace = FALSE)
  selected_exons <- sort(selected_exons)

  # Group in consecutive twos to form two-exon genes
  my_genes <- split(selected_exons, as.factor(sort(rep(seq_len(num_genes),2))))
  names(my_genes) <- paste0("gene_",seq_along(my_genes))

  export.gff3(asGFF(my_genes, parentType="gene"), gff_path)
}

simulate_annotation()

# Define function to create a simulated BAM file, barcodes file and file with the number of UMIs for testing purposes
simulate_bam <- function(reads = 100000,
                    chromosome_name = "chr1",
                    chr_length = 1000000,
                    number_barcodes = 100,
                    bam_file_path = paste0(base_dir, "possorted_genome_bam"),
                    barcodes_file_path = paste0(base_dir, "barcodes.txt"),
                    annotation_path = paste0(base_dir, "test.gff3"),
                    input_seed = 42){

  # Create a header and write out
  sam_header <- data.frame(x = '@SQ', y = paste0("SN:",chromosome_name), z = paste0("LN:",chr_length))
  write.table(sam_header, bam_file_path, sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)

  # Generate the requested number of cell barcodes and write out
  set.seed(input_seed)
  barcodes <- replicate(number_barcodes, paste(c(sample(c("A", "G", "C","T"), 16, replace = TRUE),"-1"), collapse = ""))
  write.table(data.frame(sort(barcodes)), barcodes_file_path, sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)

  # Get simulated exon and intron intervals
  txdb <- makeTxDbFromGFF(annotation_path)
  exons <- unlist(exonsBy(txdb))
  introns <- unlist(intronsByTranscript(txdb))

  # Tile exons and introns to generate reads
  exon_reads <- unlist(slidingWindows(exons, width = 100))
  intron_reads <- unlist(slidingWindows(introns, width = 100, step=5))
  #export.bed(exon_reads, paste0(base_dir,"exons.bed"))
  #export.bed(intron_reads, paste0(base_dir,"introns.bed"))

  ## Simulate BAM file in three parts; cells (60%), damaged cells (20%) & empty droplets (%20)

  # Cells ~60% exon reads, ~40% intron reads, 80% of the total reads
  cells_bam <- data.frame(
    QNAME = paste0("read", seq_len(reads*0.8)),
    FLAG = 256,
    RNAME = chromosome_name,
    POS = c(sample(start(exon_reads), 0.6*reads*0.8),
            sample(start(intron_reads), 0.4*reads*0.8)),
    MAPQ = 255,
    CIGAR = "6M",
    RNEXT = "*",
    PNEXT = 0,
    TLEN = 6,
    SEQ = "AAAAAA",
    QUAL = "FFFFFF",
    CB = paste0("CB:Z:", sample(x = barcodes[1:c(length(barcodes)*0.6)], size = reads*0.8, replace = TRUE)),
    RE = c(paste0(rep("RE:A:E", 0.6*reads*0.8)), paste0(rep("RE:A:N", 0.4*reads*0.8)))
  )

  # Empty droplets ~90% exon reads, ~10% intron reads, 5% of the total reads
  ed_bam <- data.frame(
    QNAME = paste0("read", seq_len(reads*0.05)),
    FLAG = 256,
    RNAME = chromosome_name,
    POS = c(sample(start(exon_reads), 0.9*reads*0.05),
            sample(start(intron_reads), 0.1*reads*0.05)),
    MAPQ = 255,
    CIGAR = "6M",
    RNEXT = "*",
    PNEXT = 0,
    TLEN = 6,
    SEQ = "AAAAAA",
    QUAL = "FFFFFF",
    CB = paste0("CB:Z:", sample(x = barcodes[c(length(barcodes)*0.6):c(length(barcodes)*0.8)], size = reads*0.05, replace = TRUE)),
    RE = c(paste0(rep("RE:A:E", 0.9*reads*0.05)), paste0(rep("RE:A:N", 0.1*reads*0.05)))
  )

  # Damaged cells  ~90% exon reads, ~10% intron reads, 15% of the total reads
  dc_bam <- data.frame(
    QNAME = paste0("read", seq_len(reads*0.15)),
    FLAG = 256,
    RNAME = chromosome_name,
    POS = c(sample(start(exon_reads), 0.1*reads*0.15),
            sample(start(intron_reads), 0.9*reads*0.15)),
    MAPQ = 255,
    CIGAR = "6M",
    RNEXT = "*",
    PNEXT = 0,
    TLEN = 6,
    SEQ = "AAAAAA",
    QUAL = "FFFFFF",
    CB = paste0("CB:Z:", sample(x = barcodes[c(0.8*length(barcodes)):length(barcodes)], size = reads*0.15, replace = TRUE)),
    RE = c(paste0(rep("RE:A:E", 0.1*reads*0.15)), paste0(rep("RE:A:N", 0.9*reads*0.15)))
  )

  # Combine
  bam <- do.call(rbind, list(cells_bam, ed_bam, dc_bam))
  # Sort
  bam <- bam[order(bam$POS),]

  # Make sure read names are unique
  bam$QNAME <- paste0("read", seq_len(nrow(bam)))

  # Add to header to complete the SAM file
  write.table(bam, bam_file_path, sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE, append = TRUE)

  #Convert the SAM file to a BAM file
  Rsamtools::asBam(bam_file_path, overwrite=TRUE)
  unlink(bam_file_path)

  return("Completed simulating a BAM file for testing purposes")
  }

simulate_bam()
