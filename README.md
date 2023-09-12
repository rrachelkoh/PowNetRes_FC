[![DOI](https://zenodo.org/badge/665211835.svg)](https://zenodo.org/badge/latestdoi/665211835)

# PowNetRes_FC: Coupled reservoir-power system model with streamflow forecasts

PowNetRes_FC couples a reservoir operation model with a power system model. Reservoir operations are simulated using an ad-hoc python model, and the power system model is [PowNet](https://zenodo.org/record/4688309#.YHc5euhKguU). 

# Data
PowNetRes_FC is implemented for a case study of Cambodia. 
The data used in this repository include:
- Inflow data for the six reservoirs from the Global Flood Awareness System (GloFAS)  ([Harrigan et al., 2021](https://doi.org/10.24381/cds.a4fdd6b9))
- Streamflow forecasts, which consists of an 11-member ensemble from the Global Flood Awareness System (GloFAS) ([Zsoter et al., 2020](https://doi.org/10.24381/cds.a4fdd6b9](https://doi.org/10.24381/cds.2d78664e)))
- Specifications of the power system infrastructure and demand ([EDC, 2016](https://amchamcambodia.net/wp-content/uploads/2019/08/Outlook-of-demand-and-supply-1.pdf); [JICA, 2014](https://openjicareport.jica.go.jp/644/644/644_109_12182697.html)) 


# How to run
### PowNetRes_FC comprises seven python scripts. 

#### 1. Main script: PowNetSolver.py
All other scripts can be called from PowNetSolver.py
Here we initialize the simulation and define the inputs for the models. 
- FC_typ: the type of forecast - perfect / climatology / one of the forecast members ; to run with no forecast - FC_typ = ''
- mem_num: forecast member
- res_reop: input the names of the reservoirs that we would like to re-operate (empty if no re-operation required)
- gen_cost: generation cost of each type of fuel or per-unit import cost
- opt_mtd: if forecasts is used, the optimization method is "mpc_rc", else "none". If running with forecasts, also define p1 (planning horizon) and p2 (forecast horizon)


#### 2. Data setup: DataSetup_Camb.py
Define the data source .csv files here so this script can read and convert all data into a .dat script readable and executable by the solver engine. 

#### 3. Model: PowNetModel_Camb.py
Model that defines the generation and demand nodes, as well as the system constraints. 

#### 4. Solver engine: Solver.py

### Reservoir models
#### 5. Initialization: ResModel1a.py
Initializes the reservoir parameters to be stored in two dictionaries
- res_sysparam: reservoir physical parameters
- res_ops: daily operations

#### 6. Forecast-based reservoir operations: mpc.py 
Reservoir release decisions are made based on a forecast over the next p1 days, and implemented for the next p2 days. 

#### 7. (Optional) Reservoir re-operation model: ResModel2_d.py 
This model will be activated if res_reop in PowNetSolver_FC.py is not empty. 
The goal of this model is to match the hydropower production by the reservoirs and the actual hydropower dispatched by the power system (Refer to [Koh et al., 2022](https://doi.org/10.1016/j.apenergy.2022.119386) for more details about this model). 

