result_times <- read_tsv("data_track/run_times.tsv")
mean(result_times$nuclear_fraction_run_time/(result_times$reads/100e6))
#[1] 84.06576