#17/01/21 - new Melissa try for 200 sci-MET and scNMT data
library(Melissa)
library(argparse)

parser <- ArgumentParser()

parser$add_argument("--met_dir", type="character", help="Path to the methylation directory")
parser$add_argument("--anno_file", type="character", help="Path to the annotation file")
parser$add_argument("--no_cores", type="integer", help="Number of cores for parallelisation")

args <- parser$parse_args()

input_melissa <- create_melissa_data_obj(
  met_dir = args$met_dir,
  anno_file = args$anno_file,
  chrom_size_file = NULL,
  chr_discarded = NULL,
  is_centre = TRUE,
  is_window = FALSE,
  upstream = -5000,
  downstream = 5000,
  cov = 5,
  sd_thresh = -1,
  no_cores = 5
) 

print("input_melissa object is created")

input_melissa <- filter_by_cpg_coverage(input_melissa, min_cpgcov = 3)
input_melissa <- filter_by_variability(input_melissa, min_var = 0.05)
input_melissa <- filter_by_coverage_across_cells(input_melissa, min_cell_cov_prcg = 0.05)

print("input_melissa object is filtered")

#run melissa
# Set seed for reproducible results
set.seed(15) 
# Create basis object
basis_obj <- create_rbf_object(M = 4)

print("Melissa algorithm is running")
melissa_obj <- melissa(X = input_melissa$met, K = 5, basis = basis_obj,
                       vb_max_iter = 500, vb_init_nstart = 10, 
                       is_parallel = TRUE, no_cores = args$no_cores, is_verbose = TRUE)






