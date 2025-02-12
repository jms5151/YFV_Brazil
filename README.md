# README for YFV_Brazil github repository

This repository contains code and data for a research project investigating the role of drought, mosquito dynamics, monkey populations, and human behavior in driving a Yellow Fever Virus (YFV) outbreak in Brazil, the first of its kind in over 100 years. This collaborative project involves researchers from Princeton University and Stanford University. The manuscript has been submitted for peer review under the title "Drought Dynamics Explain Once in a Century Yellow Fever Virus Outbreak in Brazil with Implications for Climate Change".

## Project Overview

The study utilizes dynamical models to simulate the complex interactions between environmental and biological factors influencing YFV transmission. The models assess how drought conditions affect mosquito populations and subsequently impact the spread of the virus among monkeys and humans.

## Repository Structure

```
YFV_Brazil/
├── .gitignore                    # Git ignore file
├── README.md                     # Project documentation
├── YFV_Brazil.Rproj              # R project file
├── assess_results.R              # Script to analyze model results
├── concat_precip.R               # Script for processing precipitation data
├── function_calculate_return_period.R # Function to calculate return periods for climatic events
├── function_extract_evspsbl.R    # Function to extract evapotranspiration and SPEI data
├── function_seasonal_forcing.R   # Function for seasonal forcing in the model
├── init_conditions.csv           # Initial conditions for model simulations
├── model.R                       # Core model script for YFV transmission dynamics
├── model_validation.csv          # Data for model validation
├── parameter_values.csv          # Parameter values used in the model
├── parameters.R                  # Script defining model parameters
├── run_model.R                   # Script to run the model simulations
└── state_variables.R             # Definitions of state variables in the model
```

## Getting Started

### Prerequisites

- **R** (version >= 4.0.0)
- Required R packages:
  - `tidyverse`
  - `data.table`
  - `foreach`
  - `doParallel`
  - `ggplot2`

You can install the necessary packages with:

```R
install.packages(c("tidyverse", "data.table", "foreach", "doParallel", "ggplot2"))
```

### Running the Model

1. Clone the repository:

```bash
git clone https://github.com/jms5151/YFV_Brazil.git
cd YFV_Brazil
```

2. Open `YFV_Brazil.Rproj` in RStudio.

3. Source the parameter and function scripts:

```R
source("parameters.R")
source("function_calculate_return_period.R")
source("function_extract_evspsbl.R")
source("function_seasonal_forcing.R")
```

4. Run the model:

```R
source("run_model.R")
```

5. Assess the results:

```R
source("assess_results.R")
```

## Data

- `init_conditions.csv`: Initial conditions for simulations.
- `parameter_values.csv`: Model parameter settings.
- `model_validation.csv`: Data used for validating model outputs.

## Results

The output includes:

- Time series of YFV cases
- Impact of drought conditions on mosquito and human populations
- Return periods for climatic events
- Visualization of model simulations

## Contributing

Contributions are welcome! Please fork the repository and submit pull requests for any improvements.

## License

This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.

## Acknowledgments

This project was supported by collaborative efforts between Princeton University and Stanford University. We thank all contributors and researchers involved in this study.

---

*For more information, please contact [your email or relevant contact information].*

