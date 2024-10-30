### Installing CmdStan and cmdstanr 
# Pre-installation requirement: C++ tool chain (for Mac e.g. Xcode) 
# If you don't have any c++ tool chain installed, see https://mc-stan.org/docs/cmdstan-guide/installation.html#macos

# Run this in a fresh R session or restarting your current session
install.packages("cmdstanr", repos = c('https://stan-dev.r-universe.dev', getOption("repos")))
# re-start R
library(cmdstanr)
# Check if you have appropriate C++ tool chain
check_cmdstan_toolchain()
# If the version is not specified, the newest available will be installed. 
# I had some issues with the newest version so I'd recommend 2.34.1!
install_cmdstan(version = "2.34.1")
# check path which cmdstan is installed at, by running cmdstan_path()


# Check if cmdstan and cmdstanr is properly installed by running an example
file <- file.path(cmdstan_path(), "examples", "bernoulli", "bernoulli.stan")
mod <- cmdstan_model(file)
mod$print()
mod$exe_file()
# names correspond to the data block in the Stan program
data_list <- list(N = 10, y = c(0,1,0,0,0,0,0,0,0,1))

fit <- mod$sample(
  data = data_list,
  seed = 123,
  chains = 4,
  parallel_chains = 4,
  refresh = 500 # print update every 500 iters
)
# If run without an error, you are good to go! From the second time, we can just run library("cmdstanr")