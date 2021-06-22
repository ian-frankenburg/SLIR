This is a quick readme demonstrating how to replicate the analysis in *A Compartment Model of Human Mobility and Early Covid-19 Dynamics in NYC*, available [here](https://arxiv.org/abs/2102.01821). The entire analysis is contained in the file `analysis.R`.

#### Software preliminaries
The core model code is written in the probabilistic programming language Stan. The model fitting is achieve through the use of CmdStan, a command-line interface avaliable for installation at https://mc-stan.org/docs/2_27/cmdstan-guide/cmdstan-installation.html. The output visuals are build using RStan, available at https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started. Once installed, one run of the `analysis.R` file will fit a vanilla SIR and the extended SLIR to the observed data and provide parameter estimates.

#### Data
The data for this analysis is supplied in the file `boroughs-case-hosp-death.csv` and was initially obtained from https://github.com/nychealth/coronavirus-data. 

#### Fitting the SLIR Model
Below is the code necessary to fit the SLIR compartment model in Stan. The inputs to the Stan model are contained in `model_data`, while the file to be compiled is passed to the `cmdstan_model` function. The SLIR model takes in five arguments: 
1. `n_obs` - length of time series
2. `N` - population size
3. `t0` - starting time
4. `cases` - observed case count data over time
5. `movement_data` - observed cellphone mobility data
```
model_data <- list(n_obs = n_obs, N = N, t0=0,
  cases = df$Cases,
  movement_data=1-df$Apple.Transit.Mobility)
file <- file.path("slir.stan")
slir.mod <- cmdstan_model(file)
slir_mcmc <- slir.mod$sample(
    data = model_data,
    max_treedepth = 10, 
    adapt_delta = .8,
    chains = 4,
    iter_warmup = 500,
    iter_sampling = 500,
    refresh = 100, 
    parallel_chains = getOption("mc.cores", 4)
)
 ```
