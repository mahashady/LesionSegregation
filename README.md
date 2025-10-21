# Segregating DNA lesions point to high selective advantage of tumor initiating cells
This repository contains the scripts to reproduce analyses described in the manuscript titled: "Segregating DNA lesions point to high selective advantage of tumor initiating cells".
<br>
Vladimir Seplyarskiy, Maha Shady, Maria A Andrianova, Michael Spencer Chapman, Eliezer M Van Allen, Shamil R Sunyaev. Vladimir Seplyarskiy, Maha Shady, Maria A Andrianova, Michael Spencer Chapman, Eliezer M Van Allen, Shamil R Sunyaev. bioRxiv 2025.10.02.680094; doi: https://doi.org/10.1101/2025.10.02.680094. (under review)

## Code
All the code required is available in this repo, and should be run in Python (version >=) or R (code developed and tested with version 4.4.2) or Perl (code developed and tested with version v5.32.1). 

### Python Dependencies
```
import argparse
import glob
import string
import gzip
import re
import os
import subprocess
import pandas as pd
import statistics
import numpy as np
import matplotlib
import scipy
import seaborn
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
## License
