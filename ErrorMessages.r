error_messages <- list(
  no_samples_error = "ERROR: There must be a data input error since there seem to be no samples. Make sure that you input a comma delimited file (.csv file) that can be read using the R command read.csv.sql. Samples correspond to columns in that file.",
  no_probes_error = "ERROR: There must be a data input error since there seem to be zero probes. Make sure that you input a comma delimited file (.csv file) that can be read using the R command read.csv.sql. CpGs correspond to rows.",
  samples_exceed_probes_warning = "MAJOR WARNING: It worries me a lot that there are more samples than CpG probes. Make sure that probes correspond to rows and samples to columns. I wonder whether you want to first transpose the data and then resubmit them? In any event, I will proceed with the analysis.",
  numeric_first_column_error = "Error: The first column does not seem to contain probe identifiers (cg numbers from Illumina) since these entries are numeric values. Make sure that the first column of the file contains probe identifiers such as cg00000292. Instead it contains ",
  non_character_first_column_warning = "Major Warning: The first column does not seem to contain probe identifiers (cg numbers from Illumina) since these entries are numeric values. Make sure that the first column of the file contains CpG probe identifiers such as cg00000292. Instead it contains ",
  input_error = "Input error. Please check the log file for details. Please read the instructions carefully.",
  non_numeric_columns = "\n MAJOR WARNING: Possible input error. The following samples contain non-numeric beta values: ",
  non_numeric_columns_hint = "\n Hint: Maybe you use the wrong symbols for missing data. Make sure to code missing values as NA in the Excel file. To proceed, I will force the entries into numeric values but make sure this makes sense.\n",
  x_missing_probes = "\n \n Comment: There are lots of missing values for X chromosomal probes for some of the samples. This is not a problem when it comes to estimating age but I cannot predict the gender of these samples.\n ",
  missing_probes_1 = "\n \n Input error: You forgot to include the following ",
  missing_probes_2 = " CpG probes (or probe names):\n "
)
