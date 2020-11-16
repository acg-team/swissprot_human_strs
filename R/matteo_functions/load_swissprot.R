#' Load SwissProt Meta Data
#'
#' @param path Path to .tsv file with UniProt/SwissProt information
#' @param tr_all a data.frame with TR results imported through \code{\link{function}}
#' @export
load_swissprot <- function(path, tr_all){
  # path = paste(local_base_path, path, sep=local_path_separator)
  sp_all = read.csv(path, sep="\t", header=TRUE, quote="", stringsAsFactors = FALSE)
  sp_all = plyr::rename(sp_all, c("Entry"="ID",
                                  "Entry.name" = "prot_name",
                                  "Gene.names" = "gene_names",
                                  "Gene.encoded.by" ="gene_encoded_by",
                                  "Cross.reference..Pfam." = "PFAM_ID",
                                  "Interacts.with" = "interacts_with",
                                  "Active.site" = "active_site",
                                  "Activity.regulation" = "activity_regulation",
                                  "Binding.site" = "binding_site",
                                  "Calcium.binding" = "calcium_binding",
                                  "Catalytic.activity" = "catalytic_activity",
                                  "DNA.binding" = "dna_binding",
                                  "EC.number" = "EC_num",
                                  "Function..CC." = "prot_function",
                                  "Metal.binding" = "metal_binding",
                                  "Nucleotide.binding" = "nucleotide_binding",
                                  "pH.dependence" = "pH_dependence",
                                  "Redox.potential" = "redox_potential",
                                  "Rhea.Ids" = "Rhea_Ids",
                                  "Temperature.dependence" = "temp_dependence",
                                  "Protein.names"="protein_name",
                                  "Virus.hosts"="virus_hosts"))

  # Add boolean to sp_all: has_tr. Later, add: has TR of specific type.
  if(! missing(tr_all)){
    sp_proteins_w_trs = unique(tr_all$ID)
    sp_proteins_w_homorep = unique(tr_all[tr_all$l_effective==1, 'ID'])
    sp_proteins_w_microsats = unique(tr_all[tr_all$l_effective<=3, 'ID'])
    sp_proteins_w_short_trs = unique(tr_all[tr_all$l_effective>3 & tr_all$l_effective<15, 'ID'])
    sp_proteins_w_domain_trs = unique(tr_all[tr_all$l_effective>=15, 'ID'])

    sp_all$has_tr = sp_all$ID %in% sp_proteins_w_trs
    sp_all$has_homo_tr = sp_all$ID %in% sp_proteins_w_homorep
    sp_all$has_micro_tr = sp_all$ID %in% sp_proteins_w_microsats
    sp_all$has_short_tr = sp_all$ID %in% sp_proteins_w_short_trs
    sp_all$has_domain_tr = sp_all$ID %in% sp_proteins_w_domain_trs
  }

  return(sp_all)
}
