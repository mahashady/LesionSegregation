# Segregating DNA lesions point to high selective advantage of tumor initiating cells
This repository contains the scripts to reproduce analyses described in the manuscript titled: "Segregating DNA lesions point to high selective advantage of tumor initiating cells".
<br>
Vladimir Seplyarskiy, Maha Shady, Maria A Andrianova, Michael Spencer Chapman, Eliezer M Van Allen, Shamil R Sunyaev. Segregating DNA lesions point to high selective advantage of tumor initiating cells. bioRxiv 2025.10.02.680094; doi: https://doi.org/10.1101/2025.10.02.680094. (under review)

## Code
All the code required is available in this repo, and should be run in Python (version >= 3.8.8) or R (code developed and tested with version 4.4.2) or Perl (code developed and tested with version v5.32.1). 

### Python Dependencies
```
import argparse
import glob
import string
import gzip
import re
import os
import subprocess
import pandas as pd # version 2.3.1
import statistics
import numpy as np # version 2.3.2
import matplotlib # version 3.10.3
import scipy # version 1.6.2
import seaborn # version 0.1..1
```
### R dependecies
```
library(HMM)
library(betareg)
library(flexmix)
library(mixtools)
library(readxl)
library(reshape2)
library(stringr)
library(viridis)
library(ggplot2)
library(this.path)
library(dplyr)
library(data.table)
library(purrr)
library(stats4)
library(tidyr)
library(ggridges)

```
## Hardware Requirements
Jobs were executed under the Slurm workload manager with the following maximum resource specifications:
```
CPU: 1 core
Memory: up to 20 GB RAM
Wall time: up to 24 hours
Partition: medium
```
## License
