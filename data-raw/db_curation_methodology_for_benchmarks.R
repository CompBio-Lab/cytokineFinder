# For the benchmarking we are only interested in CytoSig outputs
# Some cytokines may have multi-subunits 
# Created a mapping of potential gene product subunits that pertain to each cytokine specific to CytoSig

library(cytokineFindeR)

data("dbs_all")

# Cytokine to Gene Mapping - Only for CytoSig list
cytokine_to_gene_mapping <- list(
  "INHBA" = "INHBA",
  "BDNF" = "BDNF", 
  "BMP2" = "BMP2",
  "BMP4" = "BMP4",
  "BMP6" = "BMP6", 
  "BMP7" = "BMP7",
  "CCL2" = "CCL2",
  "CXCL12" = "CXCL12",
  "EGF" = "EGF",
  "EPO" = "EPO",
  "FGF2" = "FGF2",
  "GDF11" = "GDF11", 
  "HGF" = "HGF",
  "HMGB1" = "HMGB1",
  "IFNG" = "IFNG",
  "IGF1" = "IGF1",
  "IL10" = "IL10",
  "IL13" = "IL13",
  "IL15" = "IL15", 
  "IL17A" = "IL17A",
  "IL18" = "IL18",
  "IL1A" = "IL1A",
  "IL1B" = "IL1B", 
  "IL2" = "IL2",
  "IL21" = "IL21",
  "IL22" = "IL22",
  "IL3" = "IL3",
  "IL33" = "IL33",
  "IL4" = "IL4",
  "IL6" = "IL6",
  "INS" = "INS",
  "LIF" = "LIF", 
  "LTA" = "LTA",
  "NRG1" = "NRG1",
  "OSM" = "OSM",
  "PDGFB" = "PDGFB",
  "PDGFD" = "PDGFD",
  "TGFA" = "TGFA",
  "TGFB1" = "TGFB1",
  "TGFB3" = "TGFB3",
  "TSLP" = "TSLP",
  "VEGFA" = "VEGFA", 
  "WNT3A" = "WNT3A",
  "WNT5A" = "WNT5A",
  
  # Single gene conversions (name changes)
  "CD40L" = "CD40LG",
  "GCSF" = "CSF3", 
  "GMCSF" = "CSF2",
  "IFN1" = "IFNA1",
  "IFNL" = "IFNL1",
  "MCSF" = "CSF1",
  "TNFA" = "TNF",
  "TRAIL" = "TNFSF10",
  "TWEAK" = "TNFSF12",
  "PGE2" = "PTGS2",
  # Multi-subunit cytokines (only ones in CytoSig list)
  "IL12" = c("IL12A", "IL12B"),
  "IL27" = c("IL27A", "EBI3"), 
  "IL36" = c("IL36A", "IL36B", "IL36G")
)

# Single function to create CytoSig only features across databases
create_cytokine_dbs <- function(dbs_all, cytokine_mapping) {
  # For each database
  lapply(dbs_all, function(db) {
    # For each cytokine in mapping
    cytokine_features <- lapply(names(cytokine_mapping), function(cytokine) {
      genes <- cytokine_mapping[[cytokine]]
      # Get all receptors for these genes (ignore missing genes)
      all_receptors <- unlist(db[genes[genes %in% names(db)]])
      unique_receptors <- unique(all_receptors[!is.na(all_receptors)])
      # Remove only exact cytokine name matches (keep subunit interactions)
      unique_receptors[unique_receptors != cytokine]
    })
    # Name the list and remove empty features
    names(cytokine_features) <- names(cytokine_mapping)
    cytokine_features[lengths(cytokine_features) > 0]
  })
}

# Apply to create cytokine-only databases
dbs_cytosig <- create_cytokine_dbs(dbs_all, cytokine_to_gene_mapping)

usethis::use_data(dbs_cytosig, overwrite = TRUE)
