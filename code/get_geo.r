library(GEOquery)

# Create a vector of GEO dataset IDs
geo_ids <- c("GSE215039", "GSE226244","GSE160638","GSE220652","GSE235357",
"GSE220972","GSE139940","GSE16879","GSE92415","GSE16879")

# Initialize a list to store the fetched datasets
geo_datasets <- list()

# Loop through the GEO IDs
for (geo_id in geo_ids) {
  geo_data <- getGEO(geo_id, GSEMatrix = TRUE)
  geo_datasets[[geo_id]] <- geo_data
  saveRDS(geo_data, paste0(geo_id, ".rds"))
}