## Indicating which scaffolds are close for superscaffolding
# 
# need to load GenomicBreaks library first? at least coalesce ting. 

# this function seeks to give information about which contigs/scaffolds in the query genome maybe be part of a larger (super) scaffold by looking at 
# pairwise alignment information 

close_scafs <- function(gr_ob, q_genome, co_tol, ref_dist_tol, edge_tol, contig_len = 2*co_tol){
  
  # initial checks
  if(class(q_genome) == "BSgenome"){
    q_info <- seqinfo(q_genome)
  }
  else if(class(q_genome) == "Seqinfo"){
    q_info <- q_genome
  }
  else{stop("q_genome need be of either BSgenome class or Seqinfo class")}
  
  # coalesce alignments (same contig in both reference and query)
  red_gr_ob1 <- coalesce_contigs(gr_ob, tol = co_tol)
  
  # kick out small alignments
  red_gr_ob2 <- red_gr_ob1[width(ranges(red_gr_ob1)) > contig_len]
  
  # GRanges object of resulting query alignments (conversion from metadata column)
  q_ob <- GRanges(red_gr_ob2$name)
  
  # Are alignment ends close to contig ends?
  q_seqs <- as.vector(seqnames(q_ob))
  genome_seqs <- as.vector(seqnames(q_info))
  seq_lens <- as.numeric(seqlengths(q_info)[match(q_seqs, genome_seqs)])
  

  
    
}



