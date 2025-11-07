

import logging
import sys
from collections.abc import Sequence
from dataclasses import dataclass
from datetime import datetime, timezone, date, timedelta
from pathlib import Path
from time import perf_counter_ns
from typing import Dict, Optional, Union

import numpy as np
from netCDF4 import Dataset, getlibversion, num2date
from numpy.typing import NDArray

from MSU_forward_operator import MSUForwardOperator

import os

print(os.getcwd())

# import process_monthly_msu
# from process_monthly_msu import run_rtm_msu
# from process_monthly_msu import choose_freq_eia

# import RSS_surf_emiss.RSS_surf_emiss as RSS_surf_emiss
# import importlib.resources  # need this to initialize RSS_surf_emiss

from typing import NamedTuple, Optional

import numpy as np
from netCDF4 import Dataset
from numpy.typing import NDArray

import calendar
from rss_plotting.global_map import plot_global_map
import matplotlib.pyplot as plt

from era5 import era5_monthly_files,read_era5_data_monthly

def log_model_files(input_files):
    logging.info(f'land frac file: {input_files["land_frac"]}')
    logging.info(f'orography file: {input_files["orography"]}')
    logging.info(f'surface pressure file: {input_files["surf_pressure"]}')
    logging.info(f'sea ice file: {input_files["sea_ice"]}')
    logging.info(f't2m file: {input_files["t2m"]}')
    logging.info(f'tdew file: {input_files["tdew"]}')
    logging.info(f'tskin file: {input_files["tskin"]}')
    logging.info(f'w10m file: {input_files["w10m"]}')
    logging.info(f'cld file: {input_files["cld"]}')
    logging.info(f'geopot file: {input_files["geopot"]}')
    logging.info(f'q file: {input_files["q"]}')
    logging.info(f't file: {input_files["t"]}')

output_path = Path('/mnt/m/Obs4MIPs/tbs/monthly')
path_to_era5 = Path('/mnt/n/data/model/ERA5/monthly')

year_to_do = 2020

# Set up logging
log_level = logging.INFO
log_fmt = "{asctime} {levelname} {message}"

# Initialize the forward operator
forward_operator_msu2 = MSUForwardOperator('MSU2')
forward_operator_msu3 = MSUForwardOperator('MSU3')

log_datefmt = "%Y-%m-%d %H:%M:%S%z"
logging.basicConfig(style="{", format=log_fmt, datefmt=log_datefmt, level=log_level)

months_to_do = list(range(1, 13))
for month in months_to_do:

    input_files = era5_monthly_files(year_to_do, month, path_to_era5)
    log_model_files(input_files)
    
    output_file = output_path / 'tbs_msu_monthly' / 'era5' / f'y{year_to_do:04d}' / f'm{month:02d}' / f'era5_atmosphere_monthly_{year_to_do:04d}_{month:02d}.msu.ch2.nc'
    output_file.parent.mkdir(parents=True, exist_ok=True)

    logging.info(f'Processing month {month} of year {year_to_do}')
    logging.info(f'Output file: {output_file}')

    print(f'Computing Tbs for {year_to_do}-{month:02d}')
    model_data = read_era5_data_monthly(input_files)

    msu2_tbs = forward_operator_msu2.compute_tbs(model_data)
    tbs_TMT = msu2_tbs['tbs_TMT']
    tbs_TLT = msu2_tbs['tbs_TLT']

    msu3_tbs = forward_operator_msu3.compute_tbs(model_data)
    tbs_TTS = msu3_tbs['tbs_TTS']

    plot_global_map(tbs_TLT[0,:,:], vmin=230, vmax=285,
        title=f'ERA5 Simulated TB TLT {year_to_do}-{month:02d}',plt_colorbar=True,cmap='gist_ncar')
    plot_global_map(tbs_TMT[0,:,:], vmin=210, vmax=265,
        title=f'ERA5 Simulated TB TMT {year_to_do}-{month:02d}',plt_colorbar=True,cmap='gist_ncar')
    plot_global_map(tbs_TTS[0,:,:], vmin=200, vmax=240,
        title=f'ERA5 Simulated TB TTS {year_to_do}-{month:02d}',plt_colorbar=True,cmap='gist_ncar')
    plt.show()
    print()
                

