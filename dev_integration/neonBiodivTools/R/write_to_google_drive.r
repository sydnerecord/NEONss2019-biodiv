
#' Write file to Google Drive
#'
#' @title write to google drive
#'
#' @author Eric R. Sokol \email{esokol@battelleecology.org}
#'
#' @description writes to google drive
#'
#' @import googledrive readr
#'
#'
#' @param data_to_write Data frame you want to write out (dataframe or tibble).
#' @param write_filename File name to be saved (character).
#' @param my_path_to_googledirve_directory Path in Google Drive directory (google id object, or text string).
#' @param keep_local_copy_of_file Do you want to keep file locally (boolean)?
#' @param ... Other arguments passed to readr::write_csv
#'
#'
#' @references License: GNU AFFERO GENERAL PUBLIC LICENSE Version 3, 19 November 2007
#'
#'
#' @export
write_to_google_drive <- function(
  data_to_write,
  write_filename = NULL, #can be path or googleid

  # look at filenames in target directory
  my_path_to_googledirve_directory,
  keep_local_copy_of_file = TRUE){

  # get list of files
  my_list_of_files <- googledrive::drive_ls(my_path_to_googledirve_directory)

  # make a new output filename
  # temp write local
  readr::write_csv(
    data_to_write,
    write_filename,
    ...)

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
  if(!keep_local_copy_of_file) file.remove(write_filename)

}
