#' Algorithm for Coalescing Pairwise Alignments
#'
#' This algorithm take in pairwise alignment, and reduces the number of alignments by coalescing fragments that are within close proximity (user determined). Fragmented alignments can cause artificial breakpoints arising from incorrect basecalls, misassembly and misalignments.
#'
#' @param gr_ob GRanges object of the pairwise alignment, with reference genome as the subject of the GRanges object, and query genome alignment in the metadata column "name".
#' @param tol width of gap that will be bridged in coalescing. The gap must be less than or equal to \code{"tol"} in both the reference and query case.
#' @return new GRanges object of similar structure (GRanges main object is reference genome and metadata is query genome) but with a reduced number of alignment fragments due to coalescion.
#' @export
#' @import GenomicRanges
#' @import IRanges
#' @import GenomeInfoDb
#' @importFrom stats na.omit

# algorithm is vectorized for efficiency

coalesce_contigs <- function(gr_ob, tol){

  # extract query as GRanges object
  gr_q <- GRanges(gr_ob$name)

  # define new GRanges object for output
  gr_ext <- gr_ob
  q_ext <- gr_q

  # Do this in a loop, one ref scaffold at a time
  #######################################################
  ref_scafs <- unique(as.vector(seqnames(gr_ob)))
  ref_gap_sizes_total <- vector(length = length(gr_ob))
  q_gap_sizes_total <- vector(length = length(gr_ob))
  ref_con_met_total <- vector(length = length(gr_ob))
  q_con_met_total <- vector(length = length(gr_ob))
  qscaf_con_met_total <- vector(length = length(gr_ob))

  # end point considerations in query (when same scaffold is not continuous)
  c1 <- end(ranges(gr_q))[1:(length(gr_q)-1)]
  c2 <- start(ranges(gr_q))[2:length(gr_q)]
  qorder_con_met_total <- c(c1<c2, FALSE)

  k = 1 #counter
  for (i in 1:length(ref_scafs)){
    gr_now <- gr_ob[seqnames(gr_ob) == ref_scafs[i]]
    gr_len <- length(gr_now)
    if (gr_len == 1){
      k <- (k + gr_len)
      next()}

    # extract query block (blocks of the same scaffolds in query)
    q_now <- GRanges(gr_now$name)
    q_seq <- as.vector(seqnames(q_now))
    compare1 <- q_seq[2:length(q_seq)]
    compare2 <- q_seq[1:(length(q_seq) - 1)]
    compare_q <- c(compare1 == compare2, FALSE)

    # find where conditions of <tol are met in ref
    ref_starts <- start(ranges(gr_now))[2:length(gr_now)] # define starts
    ref_ends <- end(ranges(gr_now))[1:length(gr_now) - 1] # define ends
    ref_gap_sizes <- c(ref_starts - ref_ends, tol + 1)

    if(any(ref_gap_sizes < 0)){stop("gap sizes should not be negative")}

    ref_con_met <- ref_gap_sizes <= tol

    # find where conditions of <tol are met in query
    q_now <- GRanges(gr_now$name) # define GRanges object of query
    q_starts <- start(ranges(q_now))[2:length(q_now)]
    q_ends <- end(ranges(q_now))[1:length(q_now) - 1]
    q_gap_sizes <- c(q_starts - q_ends, tol + 1)


    q_con_met <- q_gap_sizes <= tol

    ref_gap_sizes_total[(k):(gr_len+k-1)] <- ref_gap_sizes
    q_gap_sizes_total[(k):(gr_len+k-1)] <- q_gap_sizes
    ref_con_met_total[(k):(gr_len+k-1)] <- ref_con_met
    q_con_met_total[(k):(gr_len+k-1)] <- q_con_met
    qscaf_con_met_total[(k):(gr_len+k-1)] <- compare_q

    k <- (k + gr_len)

  }
  #######################################################################

  # find intersection
  con_met_total = (ref_con_met_total + q_con_met_total + qscaf_con_met_total + qorder_con_met_total) == 4

  # apply extension to intersected zone (applying just to end points) (ref only)
  r_add <- ref_gap_sizes_total
  r_add[!con_met_total] = 0

  end(ranges(gr_ext)) <- end(ranges(gr_ext)) + r_add

  # reduce and concatenate
  gr_red <- reduce(gr_ext, min.gapwidth=0)

  # Use match to construct q_red
  ########################################################################
  # initialise vectors for q_red
  qr_starts <- vector(length = length(gr_red))
  qr_ends <- vector(length = length(gr_red))
  qr_seqnames <- vector(length = length(gr_red))
  qr_strand <- vector(length = length(gr_red))

  # 3 comparisons -> starts&ends, starts only, ends only (with scaffold match ofc)
  # start&end
  o_ends <- as.character(end(ranges(gr_ob)))
  o_starts <- as.character(start(ranges(gr_ob)))
  o_scafs <- as.vector(seqnames(gr_ob))
  original_comp_int <- paste(o_scafs, o_starts)
  original_se_comp <- paste(original_comp_int, o_ends)

  r_ends <- as.character(end(ranges(gr_red)))
  r_starts <- as.character(start(ranges(gr_red)))
  r_scafs <- as.vector(seqnames(gr_red))
  red_comp_int <- paste(r_scafs, r_starts)
  red_se_comp <- paste(red_comp_int, r_ends)

  se_match <- match(red_se_comp, original_se_comp)

  # starts
  s_match <- match(red_comp_int, original_comp_int)

  # ends
  original_e_comp <- paste(o_scafs, o_ends)
  red_e_comp <- paste(r_scafs, r_ends)

  e_match <- match(red_e_comp, original_e_comp)

  # construct q_red using match info
  # start&end match
  qr_starts[which(!is.na(se_match))] = start(ranges(q_ext))[na.omit(se_match)]
  qr_ends[which(!is.na(se_match))] = end(ranges(q_ext))[na.omit(se_match)]
  qr_seqnames[which(!is.na(se_match))] = as.vector(seqnames(q_ext)[na.omit(se_match)])
  qr_strand[which(!is.na(se_match))] = as.vector(strand(q_ext)[na.omit(se_match)])

  # start only match
  se_where_1 <- match(se_match, s_match)
  s_only_match <- s_match
  s_only_match[na.omit(se_where_1)] <- NA

  qr_starts[which(!is.na(s_only_match))] = start(ranges(q_ext))[na.omit(s_only_match)]
  qr_seqnames[which(!is.na(s_only_match))] = as.vector(seqnames(q_ext)[na.omit(s_only_match)])
  qr_strand[which(!is.na(s_only_match))] = as.vector(strand(q_ext)[na.omit(s_only_match)])

  # end only match
  se_where_2 <- match(se_match, e_match)
  e_only_match <- e_match
  e_only_match[na.omit(se_where_2)] <- NA

  qr_ends[which(!is.na(e_only_match))] = end(ranges(q_ext))[na.omit(e_only_match)]
  qr_seqnames[which(!is.na(se_match))] = as.vector(seqnames(q_ext)[na.omit(se_match)])
  qr_strand[which(!is.na(se_match))] = as.vector(strand(q_ext)[na.omit(se_match)])

  # final steps
  # construct seqnames Rle


  # construct strand Rle


  qr_red <- GRanges(seqnames = Rle(qr_seqnames), strand = Rle(qr_strand), IRanges(start = qr_starts, end = qr_ends))
  gr_red$name <- qr_red
  return(gr_red)
  ########################################################################

  # for theorectical testing (yay it works)
  ########################################################################
  if (FALSE){
  se_match2 <- as.vector(na.omit(match(original_se_comp, red_se_comp)))
  s_match2 <- as.vector(na.omit(match(original_comp_int, red_comp_int)))
  e_match2 <- as.vector(na.omit(match(original_e_comp, red_e_comp)))
  matches <- sort(unique(c(se_match2, s_match2, e_match2)))
  red_vec <- 1:length(gr_red)
  theory_test <- any(matches!=red_vec)
  a <- 1


  }
  ########################################################################



}
