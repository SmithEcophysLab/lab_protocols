###############################################################################
# clean_licor_file(path, write_to_csv, write_to):
###############################################################################
# Created by: Evan Perkowski (evan.a.perkowski@ttu.edu)
#
# Last edit: April 02, 2024
# Last edit by: Evan Perkowski (evan.a.perkowski@ttu.edu)
#
# Description: helper function that reads an uncleaned LI-6800 source file 
# (compatible with both a UTF-8 or .csv source file) and converts it into
# a clean and tidy file. The function currently is not operational for
# .xslx files as loading LI-6800 .xslx files seems to replace data with 
# zeroes. The helper function is implemented into a second function 
# (`clean_licor_files`) to automate the process of cleaning all LI-6800 
# source files in a working directory.
#
# Arguments:
# - path                = path of uncleaned .txt LI-6800 file
# - write_to_csv        = Boolean operator that dictates whether cleaned script 
#                         gets written as a .csv folder (write_to_csv = TRUE) or 
#                         returned as an object (write_to_csv = FALSE)
# - write_to            = path to write cleaned .txt LI-6800 file. Only used if 
#                         write_to_csv = TRUE
# Returns:
# - a .csv file written to a path of user's choosing (if write_to_csv = TRUE)
# - a data frame (if write_to_csv = FALSE)
clean_licor_file <- function(path = "",
                             write_to_csv = FALSE,
                             write_to = "") {
  
  # Determine whether read-in file is an excel file or UTF-8.
  # Boolean: if is_excel == TRUE, then read file using readxl::read_excel.
  # if is_excel == FALSE, then read file using utils::read.table
  is_excel <- grepl(".xlsx$", path)
  is_csv <- grepl(".csv$", path)
  
  if(is_excel) {
    stop("Function does not currently support .xlsx files. Convert to .csv and try again.")
  }
  
  if(is_csv) {
    data <- read.csv(file = path)
    
    # Remove all rows that occur before SystConst == "obs", but
    # keep all rows after SystConst == "obs"
    data_noheader <- data[cumsum(data$SysConst == "obs") >= 1, ]
    
    # Rename columns and remove extraneous rows that include column name,
    # column units, etc.
    names(data_noheader) <- data_noheader[1, ]
    data_clean <- data_noheader[-c(1:2), ]

  } else if(!is_csv & !is_excel){ # Read file using utils::read.table
    
    # Determine maximum number of columns (needed to designate col.names
    # in utils::read.table)
    col_number <- max(count.fields(file = path,
                                   sep = "\t",
                                   quote = ""))

    # Read file using utils::read.table
    data <- utils::read.delim(file = path, 
                              sep = "\t", 
                              row.names = NULL, 
                              col.names = 1:col_number, 
                              header = FALSE)
    
    # Remove all rows before "[Data]"
    data_noheader <- data[cumsum(data$X1 == "[Data]") >= 1, ]
    
    # Rename columns and remove extraneous rows that include column name,
    # column units, etc.
    names(data_noheader) <- data_noheader[3, ]
    data_clean <- data_noheader[-c(1:4), ]
    
    }
  
  # Conditional return (if write_to_csv = TRUE, then .csv file is written.
  # if write_to_csv = FALSE (default), then cleaned file is returned as 
  # an object)
  if(write_to_csv) { 
    
    # Create file_name path
    file_name <- paste0(basename(path), "_cleaned.csv")
    
    write.csv(data_clean,
              file = file.path(write_to, file_name),
              row.names = FALSE) }
  else{
    return(data_clean)
    }
}

###############################################################################
# clean_licor_files(directory_path, write_directory, return_list):
###############################################################################
# Created by: Evan Perkowski (evan.a.perkowski@ttu.edu)
#
# Last edit: April 02, 2024
# Last edit by: Evan Perkowski (evan.a.perkowski@ttu.edu)
#
# clean_licor_files(directory_path, write_directory, return_list):
#
# Description: helper function that reads all tab-delimited UTF-8 or .xlsx
# files in a folder directory and converts all files into a tidy, cleaned file.
# The function gives the option to write all cleaned files into a separate
# folder, or return a list of data frames with each data frame containing data
# for a single LI-6800 file.
#
# Function inputs:
#
# - directory_path    = folder containing a selection of uncleaned source 
#                       LI-6800 files
# - write_directory   = new folder where cleaned LI-6800 files are stored. Files are
#                       automatically saved as a .csv file with no option for 
#                       writing files as a different extension. This function can
#                       easily be manually altered, however, to write files to a 
#                       different extension
# - return_list       = specifies whether function returns cleaned data files 
#                       as a list of data frames (return_list = TRUE), or writes
#                       
#
# Returns:
# - if return_list = FALSE (default), the function will write a separate .csv file 
#   for each source LI-6800 file into the folder designated as the write_directory. 
#   Note that the function automatically adds "_cleaned" to the end of each .csv, 
#   but this can easily be removed by the user
# - if return_list = TRUE, the function will compile all cleaned LI-6800 files 
#   into a single list of data frames
clean_licor_files <- function(directory_path,
                              write_directory,
                              return_list = FALSE) {
  
  ## Assign directory path
  file_list <- list.files(directory_path, full.names = TRUE)
  
  ## Apply `clean_licor_file` across all files in directory_path 
  cleaned_files_list <-  lapply(file_list, 
                                function(x) 
                                  clean_licor_file(path = x, 
                                                   write_to_csv = FALSE))
  
  ## Change list elements to the basename of the file_list
  names(cleaned_files_list) <- basename(file_list)
  
  if(return_list) {
    
    return(cleaned_files_list)
  } else if(!return_list) {
    
    for(i in seq_along(cleaned_files_list)) {
      write.csv(cleaned_files_list[[i]],
                file = paste0(file.path(write_directory, 
                                        names(cleaned_files_list[i])),
                              "_cleaned.csv"), row.names = FALSE)
    }
    
  }
}