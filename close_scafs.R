## Indicating which scaffolds are close for superscaffolding
# 
# need to load GenomicBreaks library first? at least coalesce ting. 

# this function seeks to give information about which contigs/scaffolds in the query genome maybe be part of a larger (super) scaffold by looking at 
# pairwise alignment information 

## paramaters
# gr_ob - GRanges object containing pairwise alignment. Query need be in the metadata column "name"
# q_genome - genome package or seqinfo object of query genome
# co_tol - coalescing algorithm tolerance for which, if agreeable in both ref and query, alignments will be coalesced
# ref_dist_tol - gap size on reference genome for which contigs will be considered close
# edge_tol - the end of the alignment must be within this many basepairs of the end of the query scaffold in order for it to be considered 
# contig_len - alignments of length less than this will be kicked out of consideration

# Here I define the left side to be zero and right side to be length(scaffold)

close_scafs <- function(gr_ob, q_genome, co_tol, ref_dist_tol, edge_tol, contig_len = 2*edge_tol){
  
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
  
  red_gr_ob3 <- red_gr_ob2
  # GRanges object of resulting query alignments (conversion from metadata column)
  q_ob <- GRanges(red_gr_ob3$name)
  
  
  # are alignment ends close to contig/scaffold ends? If so L or R?
  q_seqs <- as.vector(seqnames(q_ob))
  genome_seqs <- as.vector(seqnames(q_info))
  seq_lens <- as.numeric(seqlengths(q_info)[match(q_seqs, genome_seqs)])
  
  sat_start_L <- start(q_ob) <= edge_tol
  sat_start_R <- start(q_ob) >= (seq_lens - edge_tol)
  
  sat_end_L <- end(q_ob) <= edge_tol
  sat_end_R <- end(q_ob) >= (seq_lens - edge_tol)
  
  
  # are alignments close on reference? We ask two questions:
  # distance small enough?
  ref_gaps <- abs(end(red_gr_ob2)[1:(length(red_gr_ob2)-1)] -
                    start(red_gr_ob2)[2:length(red_gr_ob2)]) <= ref_dist_tol
  # same reference contig?
  red2_seqs <- as.vector(seqnames(red_gr_ob2))
  sat_ref_seqs <- red2_seqs[1:(length(red2_seqs)-1)] == red2_seqs[2:length(red2_seqs)]
  
  # bring info together
  sat_ref_close <- ref_gaps + sat_ref_seqs == 2
  
  
  # correlating vectors to find valid start(L/R) and end(L/R)
  start_L <- (sat_start_L[2:length(sat_start_L)] + sat_ref_close) == 2
  start_R <- (sat_start_R[2:length(sat_start_R)] + sat_ref_close) == 2
  end_L <- (sat_end_L[1:(length(sat_end_L)-1)] + sat_ref_close) == 2
  end_R <- (sat_end_R[1:(length(sat_end_R)-1)] + sat_ref_close) == 2
  
  # which indices match up?
  #EndL - StartL
  EL_SL <- which((end_L + start_L) == 2)
  #EndL - StartR
  EL_SR <- which((end_L + start_R) == 2)
  #EndR - StartL
  ER_SL <- which((end_R + start_L) == 2)
  #EndR - StartR
  ER_SR <- which((end_R + start_R) == 2)
  
  
  # result for if they are all empty: 
  if(isEmpty(EL_SL) & isEmpty(EL_SR) & isEmpty(ER_SL) & isEmpty(ER_SR)){
    return("No two scaffolds met the requirements for superscaffolding")
  }
  
  
  # Construct results table
  tab_names <- paste(rep(as.vector(seqnames(q_info)), times = 1, each = 2), c("L", "R"))
  res_tab <- matrix(data = 0, ncol = length(tab_names), nrow = length(tab_names))
  colnames(res_tab) <- tab_names
  rownames(res_tab) <- tab_names
  
  
  # constructing vectors of matching scaffolds & fill in table, column and row wise
  #EndL - StartL
  if(!isEmpty(EL_SL)){
  EL_SL_el <- q_seqs[EL_SL]
  EL_SL_sl <- q_seqs[EL_SL + 1]
  
  EL_SL_el_tab <- paste(rep(EL_SL_el, times = 1, each = 2), c("L", "R"))
  EL_SL_sl_tab <- paste(rep(EL_SL_sl, times = 1, each = 2), c("L", "R"))
  
  res_tab[cbind(EL_SL_el_tab, EL_SL_sl_tab)] <- res_tab[cbind(EL_SL_el_tab, EL_SL_sl_tab)] + 1
  }
  #EndL - StartR
  if(!isEmpty(EL_SR)){
  EL_SR_el <- q_seqs[EL_SR]
  EL_SR_sr <- q_seqs[EL_SR + 1]
  
  EL_SR_el_tab <- paste(EL_SR_el, "L")
  EL_SR_sr_tab <- paste(EL_SR_sr, "R")
  
  res_tab[cbind(EL_SR_el_tab, EL_SR_sr_tab)] <- res_tab[cbind(EL_SR_el_tab, EL_SR_sr_tab)] + 1
  }
  #EndR - StartL
  if(!isEmpty(ER_SL)){
  ER_SL_er <- q_seqs[ER_SL]
  ER_SL_sl <- q_seqs[ER_SL + 1]
  
  ER_SL_er_tab <- paste(ER_SL_er, "R")
  ER_SL_sl_tab <- paste(ER_SL_sl, "L")
  
  res_tab[cbind(ER_SL_er_tab, ER_SL_sl_tab)] <- res_tab[cbind(ER_SL_er_tab, ER_SL_sl_tab)] + 1

  }
  #EndR - StartR
  if(!isEmpty(ER_SR)){
  ER_SR_er <- q_seqs[ER_SR]
  ER_SR_sr <- q_seqs[ER_SR + 1]
  
  ER_SR_er_tab <- paste(ER_SR_er, "R")
  ER_SR_sr_tab <- paste(ER_SR_sr, "R")
  
  res_tab[cbind(ER_SR_er_tab, ER_SR_sr_tab)] <- res_tab[cbind(ER_SR_er_tab, ER_SR_sr_tab)] + 1
  }
  
  
  # Condensed version of table
  sat_scafs <- which(res_tab!=0, arr.ind = TRUE)
  scaf_bins <- res_tab[res_tab!=0]
  res_summary <- data.frame(data = cbind(tab_names[as.vector(sat_scafs[,"col"])], tab_names[as.vector(sat_scafs[,"row"])],
                                     as.numeric(scaf_bins)))
  
  # finish filling in results table so that it is symmetric along the diagonal
  if(!isEmpty(EL_SL)){
    res_tab[cbind(EL_SL_sl_tab, EL_SL_el_tab)] <- res_tab[cbind(EL_SL_sl_tab, EL_SL_el_tab)] + 1
  }
  if(!isEmpty(EL_SR)){
    res_tab[cbind(EL_SR_sr_tab, EL_SR_el_tab)] <- res_tab[cbind(EL_SR_sr_tab, EL_SR_el_tab)] + 1
  }
  if(!isEmpty(ER_SL)){
    res_tab[cbind(ER_SL_sl_tab, ER_SL_er_tab)] <- res_tab[cbind(ER_SL_sl_tab, ER_SL_er_tab)] + 1
  }
  if(!isEmpty(ER_SR)){
    res_tab[cbind(ER_SR_sr_tab, ER_SR_er_tab)] <- res_tab[cbind(ER_SR_sr_tab, ER_SR_er_tab)] + 1
  }
  return(list(res_tab, res_summary))
}
  
