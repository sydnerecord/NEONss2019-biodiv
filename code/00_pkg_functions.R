# Packages used
library(tidyverse)
library(googledrive)
# install_github("EDIorg/ecocomDP", ref = 'development')
library(ecocomDP)
library(neonUtilities)

#' Write file to Google Drive
#' 
#' @param data_to_write Data frame you want to write out.
#' @param write_filename File name to be saved.
#' @param my_path_to_googledirve_directory Path in Google Drive directory.
#' @param keep_file_in_working_dir Do you want to keep file locally?
#' 
write_to_google_drive <- function(
  data_to_write,
  write_filename,
  my_path_to_googledirve_directory, 
  keep_file_in_working_dir = FALSE 
){
  # get list of files
  my_list_of_files <- googledrive::drive_ls(my_path_to_googledirve_directory)
  
  # make a new output filename
  # temp write local
  readr::write_csv(data_to_write, write_filename)
  
  # conditional depending on if we need to overwrite or create new
  if(!write_filename %in% my_list_of_files$name){
    drive_upload(write_filename,
                 path = my_path_to_googledirve_directory,
                 name = write_filename,
                 type = NULL,
                 verbose = TRUE)
    message(paste0('Created ',write_filename, ' in ', my_path_to_googledirve_directory))
  }else{
    google_id <- my_list_of_files %>% filter(name == write_filename) %>% select(id) %>% unlist()
    drive_update(file = as_id(google_id),
                 media = write_filename)
    message(paste0('Updated ',write_filename, ' in ', my_path_to_googledirve_directory))
  }
  
  #remove local file
  if(!keep_file_in_working_dir) file.remove(write_filename)
  
}
