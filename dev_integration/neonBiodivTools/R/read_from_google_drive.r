#' Read file from Google Drive
#'
#' @title read from google drive
#'
#' @author Eric R. Sokol \email{esokol@battelleecology.org}
#'
#' @description read a csv or txt file from google drive
#'
#' @import googledrive readr
#' 
#' 
#' @param file_name_string Data frame you want to write out (dataframe or tibble).
#' @param my_path_to_googledirve_directory Path in Google Drive directory (google id object, or text string).
#' @param keep_local_copy_of_file Do you want to keep file locally (boolean, default = TRUE)?
#' @param ... other arguments passed to readr::read_csv
#'
#'
#' @references License: GNU AFFERO GENERAL PUBLIC LICENSE Version 3, 19 November 2007
#'
#'
#' @export
read_from_google_drive <- function(
  file_name_string = NULL, #can be path or googleid, e.g., for path 'NAQWA_algae_derived_biodiversity_metrics_by_continetnal_US.csv'
  # look at filenames in target directory
  my_path_to_googledirve_directory = NULL,
  keep_local_copy_of_file = TRUE,
  ...){

  #####################
  # download most recent and reading in the raw data file
  
  
  # get list of files
  my_list_of_files <- googledrive::drive_ls(my_path_to_googledirve_directory)
  
  matching_files <- my_list_of_files %>% filter(grepl(file_name_string,name))
  
  if(nrow(matching_files)!=1){
    stop('found ', nrow(matching_files), ' matching files')
  }
  
  my_google_id <- matching_files %>% 
    slice(1) %>% 
    select(id) %>% unlist(use.names = FALSE) %>% 
    googledrive::as_id()
  
  # download file to local working dir
  downloaded_file <- googledrive::drive_download(
    my_google_id, overwrite = TRUE)
  
  # read into R
  my_result <- data.frame()
  my_result <- readr::read_csv(downloaded_file$local_path, ...)
  
  # delete local file 
  if(!keep_local_copy_of_file) file.remove(downloaded_file$local_path)
  
  return(my_result)
  
}