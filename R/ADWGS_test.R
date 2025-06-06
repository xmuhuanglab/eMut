#
#' Makes mutational signatures
#'
#' @return a dataframe with mutational signatures
.make_mut_signatures = function() {
  nucs = c("A", "C", "G", "T")

  signatures_dfr = data.frame(as.matrix(expand.grid(left=nucs, mid=nucs, right=nucs, alt=nucs)), stringsAsFactors=F)
  signatures_dfr = signatures_dfr[signatures_dfr$mid!=signatures_dfr$alt,]
  signatures_dfr = signatures_dfr[signatures_dfr$mid %in% c("C", "T"),]
  signatures_dfr$tri = paste0(signatures_dfr$left, signatures_dfr$mid, signatures_dfr$right)
  signatures_dfr$signt = paste0(signatures_dfr$tri, ">", signatures_dfr$alt)
  signatures_dfr
}




#' Calculates the number of expected mutations based
#'
#' @param hyp hypothesis to be tested
#' @param select_positions boolean column which indicates which positions are in the element of interest
#' @param dfr a dataframe containing the data to be tested
#' @param colname name of the column which indicates the count of mutations in the positions of interest
#'
#' @return a list of observed mutations and expected mutations
.get_obs_exp = function(hyp, select_positions, dfr, colname) {
  obs_mut = sum(dfr[select_positions, colname])

  exp_probs = hyp$fitted.values[select_positions]
  # derive a fractional estimate of expected mutations as a sum of per-nucleotide probabilties
  exp_mut = sum(exp_probs)

  list(obs_mut, exp_mut)
}

# @import GenomicRanges
# @import GenomeInfoDb
# @import BSgenome.Hsapiens.UCSC.hg19
# @import IRanges
# @import BSgenome
# @import S4Vectors
# @import ActiveDriverWGS

#' ADWGS_test executes the statistical test for ActiveDriverWGS
#'
#' @param id A string used to identify the element of interest. \code{id}
#' corresponds to an element in the id column of the elements file
#' @param gr_element_coords A GenomicRanges object that describes the elements of interest containing the
#' chromosome, start and end coordinates, and an mcols column corresponding to id
#' @param gr_site_coords  A GenomicRanges object that describes the sites of interest which reside
#' in the elements of interest containing the chromosome, start and end coordinates,
#' and an mcols column corresponding to id. Examples of sites include transcription factor binding
#' sites in promoter regions or phosphosites in exons of protein coding genes. An empty GenomicRanges object
#' nullifies the requirement for sites to exist.
#' @param gr_maf A GenomicRanges object that describes the mutations in the dataset containing the chromosome,
#' start and end coordinates, patient id, and trinucleotide context
#' @param win_size An integer indicating the size of the background window in base pairs that is used to establish
#' the expected mutation rate and respective null model. The default is 50000bps
#' @param this_genome The reference genome object of BSgenome, for example BSgenome.Hsapiens.UCSC.hg19::Hsapiens
#' @param detect_depleted_mutations if TRUE, detect elements with significantly fewer than expected mutations. FALSE by default
#' @param openRegions=NULL if set, open Regions as the background region
#'
#' @return A data frame containing the following columns
#' \describe{
#'     \item{id}{A string identifying the element of interest}
#'     \item{pp_element}{The p-value of the element}
#'     \item{element_muts_obs}{The number of patients with a mutation in the element}
#'     \item{element_muts_exp}{The expected number of patients with a mutation in the element with respect to background}
#'     \item{element_enriched}{A boolean indicating whether the element is enriched in mutations}
#'     \item{pp_site}{The p-value of the site}
#'     \item{site_muts_obs}{The number of patients with a mutation in the site}
#'     \item{site_muts_exp}{The expected number of patients with a mutation in the site with respect to element}
#'     \item{site_enriched}{A boolean indicating whether the site is enriched in mutations}
#'     \item{result_number}{A numeric indicator denoting the order in which the results were calculated}
#'     \item{fdr_element}{The FDR corrected p-value of the element}
#'     \item{fdr_site}{The FDR corrected p-value of the site}
#'     \item{has_site_mutations}{A V indicates the presence of site mutations}
#' }
#'
#' @export
#'
#'

ADWGS_test = function(id, gr_element_coords, gr_site_coords, gr_maf, win_size, this_genome, detect_depleted_mutations = FALSE, openRegions=NULL) {

  require(ActiveDriverWGS)
  require(GenomicRanges)

  cat(".")
  null_res = data.frame(id,
                        pp_element = NA, element_muts_obs = NA, element_muts_exp = NA, element_enriched = NA,
                        pp_site = NA, site_muts_obs = NA, site_muts_exp = NA, site_enriched = NA,
                        stringsAsFactors = F)

  ## genome element arithmethics
  gr_elements = gr_element_coords[GenomicRanges::mcols(gr_element_coords)[,1]==id]
  # sites can be missing
  if (length(gr_site_coords) == 0) {
    gr_sites = GenomicRanges::GRanges()
  } else {
    gr_sites = gr_site_coords[GenomicRanges::mcols(gr_site_coords)[,1]==id]
  }

  # return empty result if no mutations found
  if (length(GenomicRanges::findOverlaps(gr_elements, gr_maf)) == 0) {
    return(null_res)
  }

  # take whole sequence around element and account for chromosome ends
  gr_background = .create_background(gr_elements, win_size, this_genome)

  #######  (modified : add openRegions as background)
  if(!is.null(openRegions)){
    gr_background = GenomicRanges::intersect(gr_background, openRegions)
  }

  # take all mutations needed for this analysis (both element and background)
  gr_mutations = gr_maf[S4Vectors::queryHits(GenomicRanges::findOverlaps(gr_maf, gr_background))]
  # separate indels and SNVs for distinct rules
  gr_indel = gr_mutations[GenomicRanges::mcols(gr_mutations)[,2] == "indel>X"]
  gr_snv = gr_mutations[GenomicRanges::mcols(gr_mutations)[,2] != "indel>X"]

  # remove element components from background
  gr_background = GenomicRanges::setdiff(gr_background, gr_elements)

  # make sure sites capture only the element and not the background
  # there is an unexplained err here with chromosomes/levels:
  # when intersect(elements, sites), then err is thrown: seqlevels(seqinfo(x))' and 'levels(seqnames(x))' are not identical
  # if order switched, then err no longer
  gr_sites = GenomicRanges::intersect(gr_elements, gr_sites)

  # remove sites from elements to avoid double counting
  gr_elements = GenomicRanges::setdiff(gr_elements, gr_sites)

  # create set of all mutation signatures (trinucleotide ref + alt nucleotide); AG mapped to CG
  signt_template = .make_mut_signatures()

  # get trinucleotide content of sequences
  trinuc_elements = .seq2signt(gr_elements, this_genome, signt_template)
  trinuc_sites = .seq2signt(gr_sites, this_genome, signt_template)
  trinuc_background = .seq2signt(gr_background, this_genome, signt_template)

  # separate indels and SNVs for distinct analyses
  dfr_snv = NULL; patients_with_SNV_in_element = NULL
  if(length(gr_snv) > 0) {

    # count the mutations of every signt class
    snv_elements = .mut2signt(gr_elements, gr_snv, signt_template,
                              remove_dup_mut_per_patient = TRUE)
    snv_background = .mut2signt(gr_background, gr_snv, signt_template)

    # merge mutation counts and trinucleotide counts
    snv_elements = merge(snv_elements, trinuc_elements, by = "signt")
    snv_background = merge(snv_background, trinuc_background, by = "signt")
    snv_elements$region = "elements"
    snv_background$region = "background"

    # sites are optional - sometimes not provided and
    # sometimes not available for given region
    snv_sites = NULL
    if (length(gr_sites) > 0) {
      snv_sites = .mut2signt(gr_sites, gr_snv, signt_template,
                             remove_dup_mut_per_patient = TRUE)
      snv_sites = merge(snv_sites, trinuc_sites, by = "signt")
      snv_sites$region = "sites"
    }

    # keep track of patients with SNVs, to remove indels later
    gr_snv_el_site = gr_snv[unique(S4Vectors::queryHits(
      GenomicRanges::findOverlaps(gr_snv, c(gr_elements, gr_sites))))]
    patients_with_SNV_in_element =
      unique(GenomicRanges::mcols(gr_snv_el_site)[,"mcols.patient"])

    dfr_snv = rbind(snv_sites, snv_elements, snv_background)
  }

  # separate indels for distinct analyses
  dfr_indel = NULL
  if(length(gr_indel) > 0) {

    # remove multiple indels per patient if these affect the same element/site
    gr_indel_fg = gr_indel[S4Vectors::queryHits(
      GenomicRanges::findOverlaps(gr_indel, c(gr_sites, gr_elements)))]
    gr_indel_fg = gr_indel_fg[!duplicated(
      GenomicRanges::mcols(gr_indel_fg)[,"mcols.patient"])]
    gr_indel_bg = gr_indel[S4Vectors::queryHits(
      GenomicRanges::findOverlaps(gr_indel, gr_background))]
    gr_indel = c(gr_indel_fg, gr_indel_bg)

    # unique indel tag, remove multiple instances
    indel_tag = apply(data.frame(gr_indel)[,c("seqnames", "start", "end", "mcols.patient")],
                      1, paste, collapse = "::")
    gr_indel = gr_indel[!duplicated(indel_tag)]

    # count all indels in regions
    indel_index_sites = unique(S4Vectors::queryHits(
      GenomicRanges::findOverlaps(gr_indel, gr_sites)))
    indel_index_elements = unique(S4Vectors::queryHits(
      GenomicRanges::findOverlaps(gr_indel, gr_elements)))
    indel_index_background = unique(S4Vectors::queryHits(
      GenomicRanges::findOverlaps(gr_indel, gr_background)))

    # to avoid double counting, create hierarchy: sites affected first, then elements, bgrd
    indel_index_elements = setdiff(indel_index_elements, indel_index_sites)
    indel_index_background = setdiff(indel_index_background,
                                     c(indel_index_elements, indel_index_sites))

    # discard element indels that occur in a patient previously seen with SNV there
    which_indel_index_element_dup =
      which(GenomicRanges::mcols(gr_indel[indel_index_elements])[,"mcols.patient"]
            %in% patients_with_SNV_in_element)
    if (length(which_indel_index_element_dup) > 0) {
      indel_index_elements = indel_index_elements[ - which_indel_index_element_dup]
    }
    # discard site indels that occur in a patient previously seen with SNV there
    which_indel_index_sites_dup =
      which(GenomicRanges::mcols(gr_indel[indel_index_sites])[,"mcols.patient"]
            %in% patients_with_SNV_in_element)
    if (length(which_indel_index_sites_dup) > 0) {
      indel_index_sites = indel_index_sites[ - which_indel_index_sites_dup]
    }

    # initiate data frame
    indel_sites = indel_elements = indel_background = data.frame(
      signt = "indel>X", n_mut = NA, tri_nucleotide = "indel",
      n_pos = NA, region = NA, stringsAsFactors = FALSE)

    # set number of mutations in the regions, as defined above hierarchically
    indel_elements$n_mut = length(indel_index_elements)
    indel_background$n_mut = length(indel_index_background)

    # set number of positions as all positions in the regions
    # since indel assumed to affect any site
    indel_elements$n_pos = sum(GenomicRanges::width(gr_elements))
    indel_background$n_pos = sum(GenomicRanges::width(gr_background))

    indel_elements$region = "elements"
    indel_background$region = "background"

    if (length(gr_sites) > 0) {
      indel_sites$n_mut = length(indel_index_sites)
      indel_sites$n_pos = sum(GenomicRanges::width(gr_sites))
      indel_sites$region = "sites"
    } else {
      indel_sites = NULL
    }

    dfr_indel = rbind(indel_sites, indel_elements, indel_background)
  }

  # merge mutations and annotate per region type
  dfr_mut = rbind(dfr_indel, dfr_snv)
  dfr_mut$is_site = 0 + (dfr_mut$region == "sites")
  dfr_mut$is_element = 0 + dfr_mut$region %in% c("sites", "elements")

  # remove nucleotide contexts that are never mutated
  signt_with_muts = names(which(c(by(dfr_mut$n_mut, dfr_mut$signt, sum)) > 0))
  dfr_mut = dfr_mut[dfr_mut$signt %in% signt_with_muts,, drop = FALSE]

  # remove nucleotide contexts that have no positions
  dfr_mut = dfr_mut[dfr_mut$n_pos > 0,, drop = FALSE]

  # apply a simpler formula if only one trinucleotide signature is mutated
  formula_h0 = ifelse(length(signt_with_muts) > 1, "n_mut ~ signt", "n_mut ~ 1")

  # Poisson model testing
  h0 = stats::glm(stats::as.formula(formula_h0), offset = log(dfr_mut$n_pos),
                  family = stats::poisson, data = dfr_mut)
  h1 = stats::update(h0, . ~ . + is_element)
  pp_element = pp_element_2way = stats::anova(h0, h1, test="Chisq")[2,5]

  # coefficient determines enrichment or depletion
  # if significant depletion flip p-value direction
  # this is the default option.
  # if user defines 'allow_negative_selection' then both positively and negatively selected sites will be shown
  coef_element = stats::coef(h1)[['is_element']]
  element_enriched = coef_element > 0
  if (!detect_depleted_mutations & !element_enriched & !is.na(pp_element) & pp_element < 0.5) {
    pp_element = 1 - pp_element_2way
  }
  if (detect_depleted_mutations & element_enriched & !is.na(pp_element) & pp_element < 0.5) {
    pp_element = 1 - pp_element_2way
  }

  # observed and expected values from sampling
  element_stats = .get_obs_exp(h0, dfr_mut$is_element == 1, dfr_mut, "n_mut")
  element_muts_obs = element_stats[[1]]
  element_muts_exp = element_stats[[2]]

  # if region has sites, test second hypothesis on regions
  pp_site = site_muts_obs = site_muts_exp = site_enriched = site_depleted = NA
  if (length(gr_sites) > 0) {

    h2 = stats::update(h1, . ~ . + is_site)
    pp_site = pp_site_2way = stats::anova(h1, h2, test="Chisq")[2,5]

    # coefficient determines enrichment or depletion
    coef_site = stats::coef(h2)[['is_site']]
    site_enriched = coef_site > 0
    if (!detect_depleted_mutations & !site_enriched & !is.na(pp_site) & pp_site < 0.5) {
      pp_site = 1 - pp_site_2way
    }
    if (detect_depleted_mutations & site_enriched & !is.na(pp_site) & pp_site < 0.5) {
      pp_site = 1 - pp_site_2way
    }

    site_stats = .get_obs_exp(h1, dfr_mut$is_site == 1, dfr_mut, "n_mut")
    site_muts_obs = site_stats[[1]]
    site_muts_exp = site_stats[[2]]
  }

  data.frame(id,
             pp_element, element_muts_obs, element_muts_exp, element_enriched,
             pp_site, site_muts_obs, site_muts_exp, site_enriched,
             stringsAsFactors = F)
}



# gr_seq = gr_elements;
# gr_seq = gr_sites
.seq2signt = function(gr_seq, this_genome, signt_template) {

  if (length(gr_seq) == 0) {
    return(NULL)
  }

  # add one bp on either side to capture flanking nucs for trinucleotide
  gr_seq1 = GenomicRanges::GRanges(GenomicRanges::seqnames(gr_seq),
                                   IRanges::IRanges(
                                     GenomicRanges::start(gr_seq) - 1,
                                     GenomicRanges::end(gr_seq) + 1))

  # get sequence trinucleotide
  this_seq = strsplit(as.character(BSgenome::getSeq(this_genome, gr_seq1)), '')

  # extract trinucleotides from sequence
  get_sq_trinucs = function(x) sapply(1:(length(x)-2), function(i) paste0(x[i:(i+2)], collapse=""))
  tri_nucleotide = unlist(lapply(this_seq, get_sq_trinucs))

  # complement trinucleotides for A/G nucleotides
  where_reverse = substr(tri_nucleotide, 2, 2) %in% c("A", "G")
  tri_nucleotide[where_reverse] =
    as.character(Biostrings::reverseComplement(Biostrings::DNAStringSet(
      tri_nucleotide[where_reverse])))

  tri_nucleotide = data.frame(table(tri_nucleotide), stringsAsFactors = FALSE)
  tri_nucleotide[,1] = as.character(tri_nucleotide[,1])

  quad_nucl =  merge(tri_nucleotide, signt_template[,c("tri", "signt")],
                     by.x = "tri_nucleotide", by.y = "tri", all = TRUE)

  # make sure non-represented nucleotides are counted as zeroes
  quad_nucl$Freq[is.na(quad_nucl$Freq)] = 0

  # update column name of frequency
  colnames(quad_nucl)[colnames(quad_nucl) == "Freq"] = "n_pos"

  quad_nucl
}



# gr_seq = gr_elements
# gr_seq = gr_sites
# gr_seq = gr_background
.mut2signt = function(gr_seq, gr_snv, signt_template, remove_dup_mut_per_patient = FALSE) {

  if (length(gr_seq) == 0) {
    return(NULL)
  }

  gr_snv_here =
    gr_snv[S4Vectors::queryHits(GenomicRanges::findOverlaps(gr_snv, gr_seq))]

  # for elements and sites, we need to remove duplicate mutations per patient
  # false positives such as local hypermutation rates
  if (remove_dup_mut_per_patient) {
    gr_snv_here =
      gr_snv_here[!duplicated(GenomicRanges::mcols(gr_snv_here)[,"mcols.patient"])]
  }

  # if there are no mutations, return signts with all zeroes as counts
  if (length(gr_snv_here) == 0 ) {
    return(data.frame(signt = signt_template$signt, n_mut = 0,
                      stringsAsFactors = FALSE))
  }

  dfr_snv = data.frame(table(GenomicRanges::mcols(gr_snv_here)[,"mcols.tag"]),
                       stringsAsFactors = FALSE)
  quad_nucl =  merge(signt_template[,c("signt"), drop = FALSE], dfr_snv,
                     by.x = "signt", by.y = "Var1", all = TRUE)

  # make sure non-represented nucleotides are counted as zeroes
  quad_nucl$Freq[is.na(quad_nucl$Freq)] = 0

  # update column name of frequency
  colnames(quad_nucl)[colnames(quad_nucl) == "Freq"] = "n_mut"

  quad_nucl
}



.create_background = function(gr_elements, win_size, this_genome) {

  # background is plus/minus window around every component of element
  # make sure flanking coordinates stay within chromosomal boundaries
  # add one extra nucleotide of slack on both ends - trinuc tabulation needs that
  bg_starts = GenomicRanges::start(gr_elements) - win_size
  bg_ends = GenomicRanges::end(gr_elements) + win_size
  # max pos is end of chromosome minus one
  max_chr_pos = GenomeInfoDb::seqlengths(this_genome)[as.character(GenomicRanges::seqnames(gr_elements))]
  max_chr_pos = max_chr_pos - 1
  # min pos is 2nd pos of chromosome
  bg_starts[bg_starts < 2] = 2
  bg_ends[bg_ends > max_chr_pos] = max_chr_pos[bg_ends > max_chr_pos]
  gr_background = GenomicRanges::GRanges(GenomicRanges::seqnames(gr_elements),
                                         IRanges::IRanges(bg_starts, bg_ends))
  # take one joined background set to avoid duplicates
  gr_background = GenomicRanges::union(gr_background, gr_background)
  gr_background
}


#' fix_all_results verifies that the results table has the correct format and p-values
#'
#' @param all_results a data frame containing the following columns
#' \describe{
#'     \item{id}{A string identifying the element of interest}
#'     \item{pp_element}{The p-value of the element}
#'     \item{element_muts_obs}{The number of patients with a mutation in the element}
#'     \item{element_muts_exp}{The expected number of patients with a mutation in the element with respect to background}
#'     \item{element_enriched}{A boolean indicating whether the element is enriched in mutations}
#'     \item{pp_site}{The p-value of the element}
#'     \item{site_muts_obs}{The number of patients with a mutation in the site}
#'     \item{site_muts_exp}{The expected number of patients with a mutation in the site with respect to element}
#'     \item{site_enriched}{A boolean indicating whether the site is enriched in mutations}
#'     \item{result_number}{A numeric indicator denoting the order in which the results were calculated}
#' }
#'
#' @return the same data frame
.fix_all_results = function(all_results) {

  resi = data.frame(all_results, stringsAsFactors=F)
  resi$pp_site = as.numeric(resi$pp_site)
  resi$pp_element = as.numeric(resi$pp_element)
  resi$element_muts_obs = as.numeric(resi$element_muts_obs)
  resi$element_muts_exp = as.numeric(resi$element_muts_exp)
  resi$element_enriched = as.logical(gsub("\\s+", "", resi$element_enriched))
  resi$site_muts_obs = as.numeric(resi$site_muts_obs)
  resi$site_muts_exp = as.numeric(resi$site_muts_exp)
  resi$site_enriched = as.logical(gsub("\\s+", "", resi$site_enriched))

  resi[!is.na(resi$pp_element) & resi$pp_element==0,"pp_element"] = 1e-300
  resi[!is.na(resi$pp_site) & resi$pp_site==0,"pp_site"] = 1e-300

  resi_tag = resi[,"id"]
  resi = resi[!duplicated(resi_tag),]
  resi
}

# @import stats

#' Returns significant results
#'
#' @param all_res a data frame containing the following columns
#' \describe{
#'     \item{id}{A string identifying the element of interest}
#'     \item{pp_element}{The p-value of the element}
#'     \item{element_muts_obs}{The number of patients with a mutation in the element}
#'     \item{element_muts_exp}{The expected number of patients with a mutation in the element with respect to background}
#'     \item{element_enriched}{A boolean indicating whether the element is enriched in mutations}
#'     \item{pp_site}{The p-value of the element}
#'     \item{site_muts_obs}{The number of patients with a mutation in the site}
#'     \item{site_muts_exp}{The expected number of patients with a mutation in the site with respect to element}
#'     \item{site_enriched}{A boolean indicating whether the site is enriched in mutations}
#'     \item{result_number}{A numeric indicator denoting the order in which the results were calculated}
#' }
#'
#' @return the same data frame with three addition columns
#' \describe{
#'     \item{fdr_element}{The FDR corrected p-value of the element}
#'     \item{fdr_site}{The FDR corrected p-value of the site}
#'     \item{has_site_mutations}{A V indicates the presence of site mutations}
#' }
#'
.get_signf_results = function(all_res) {
  this_results = all_res
  if (nrow(this_results)==0) {
    return(NULL)
  }
  # this is FDR treating element-level NAs as 1s
  this_results$fdr_element = stats::p.adjust(this_results$pp_element, method="fdr", n=nrow(this_results))

  this_results = this_results[order(this_results$fdr_element),]

  filtered_results = this_results[!is.na(this_results$fdr_element) & this_results$fdr_element<0.05,]
  unsignf_results = this_results[is.na(this_results$fdr_element) | this_results$fdr_element>=0.05,]

  if (nrow(filtered_results) + nrow(unsignf_results) != nrow(this_results)) {
    stop("Error: Something unexpected happened when formatting results")
  }

  if (nrow(filtered_results)!=0) {
    # site-level FDR perform only on elements with pre-selection of FDR<0.05
    filtered_results$fdr_site = stats::p.adjust(filtered_results$pp_site, method="fdr", n=nrow(filtered_results))
    filtered_results$has_site_mutations = !is.na(filtered_results$fdr_site) & filtered_results$fdr_site<0.05
    filtered_results$has_site_mutations = c("","V")[1+c(filtered_results$has_site_mutations)]
  }

  if (nrow(unsignf_results) !=0) {
    unsignf_results$fdr_site = NA
    unsignf_results$has_site_mutations = ""
  }

  final_results = rbind(filtered_results, unsignf_results)
  final_results$pp_element = replace(final_results$pp_element, which(is.na(final_results$pp_element)), 1)
  final_results$pp_site = replace(final_results$pp_site, which(is.na(final_results$pp_site)), 1)
  final_results$fdr_element = replace(final_results$fdr_element, which(is.na(final_results$fdr_element)), 1)
  final_results$fdr_site = replace(final_results$fdr_site, which(is.na(final_results$fdr_site)), 1)
  rownames(final_results) = NULL
  final_results
}
