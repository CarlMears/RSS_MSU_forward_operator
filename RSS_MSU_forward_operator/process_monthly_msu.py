"""Process a monthly file.

"""

import argparse
import logging
import sys
from collections.abc import Sequence
from dataclasses import dataclass
from datetime import datetime, timezone
import calendar
from pathlib import Path
from time import perf_counter_ns
from typing import Optional, Union

import numpy as np
from netCDF4 import Dataset, getlibversion, num2date
from numpy.typing import NDArray

import era5
import rss_atmos_rtm
dir(rss_atmos_rtm)
print(rss_atmos_rtm.__file__)
from rss_atmos_rtm.rss_atmos_rtm import compute_rtm


REF_FREQ_MSU1 = np.array([50.30,50.30,50.30,50.30,50.30,50.30], np.float32)
REF_FREQ_MSU2 = np.array([53.74,53.74,53.74,53.74,53.74,53.74], np.float32)
REF_FREQ_MSU3 = np.array([54.96,54.96,54.96,54.96,54.96,54.96], np.float32)
REF_FREQ_MSU4 = np.array([57.95,57.95,57.95,57.95,57.95,57.95], np.float32)

#incidence angles for all MSU channels are the same
REF_EIA_MSU =  np.array([ 0.00,10.71,21.51,32.51,43.91,56.19], np.float32)
REF_LOOK_ANGLE_MSU = np.array([0.00, 9.47, 18.94, 28.41, 37.88, 47.35], np.float32)

def choose_freq_eia(self, msu_channel: str) -> None:
    """Set the frequency and incidence angle based on the MSU channel."""
    channel_info = {}
    if msu_channel in ['MSU1', 'MSU2', 'MSU3', 'MSU4']:
        channel_info['num_freq'] = 6
        channel_info['num_eia'] = 6
        channel_info['satellite'] = 'MSU'
        channel_info['look_angle'] = REF_LOOK_ANGLE_MSU
        channel_info['incidence'] = REF_EIA_MSU
    if msu_channel == 'MSU1':
        channel_info['frequency'] = REF_FREQ_MSU1
        channel_info['polarization'] = 'V'
    elif msu_channel == 'MSU2':
        channel_info['frequency'] = REF_FREQ_MSU2
        channel_info['polarization'] = 'H'
    elif msu_channel == 'MSU3':
        channel_info['frequency'] = REF_FREQ_MSU3
        channel_info['polarization'] = 'V'
    elif msu_channel == 'MSU4':
        channel_info['frequency'] = REF_FREQ_MSU4
        channel_info['polarization'] = 'H'
    else:
        raise ValueError(f"Unsupported MSU channel: {msu_channel}")
    
    return channel_info

@dataclass
class RtmMonthlyData:
    """Output values after computing the atmospheric RTM.

    This is for a full day of data, so it is more "full" than should normally be
    computed.
    """
    year: int
    month: int

    # Latitude in degrees North, dimensioned as (num_lats, ). They should be in
    # ascending order (e.g., -90 to 90).
    lats: NDArray[np.float32]

    # Longitude in degrees East, dimensioned as (num_lons, ). They should be in
    # ascending order (e.g., 0 to 360).
    lons: NDArray[np.float32]

    # Microwave frequency in GHz, dimensioned as (num_freq, ).
    frequency: NDArray[np.float32]

    # Earth incidence angle in degrees, dimensioned as (num_freq, ).
    incidence: NDArray[np.float32]

    # Atmospheric transmissivity, unitless, dimensioned as (time, lats, lons,
    # freq).
    transmissivity: NDArray[np.float32]

    # Atmospheric upwelling brightness temperature in K, dimensioned as (time,
    # lats, lons, freq).
    tb_up: NDArray[np.float32]

    # Atmospheric downwelling brightness temperature in K, dimensioned as (time,
    # lats, lons, freq).
    tb_down: NDArray[np.float32]



    def write_nc(self, rtm_output: Path) -> None:
        """Write the daily RTM data to a netCDF output file."""
        timestamp = datetime.now(timezone.utc).isoformat(" ", "seconds")
        nc_version = getlibversion().partition(" ")[0]

        time_units = "hours since 1900-01-01 00:00:00Z"
        data_time = (datetime(self.year, self.month, 1, tzinfo=timezone.utc) -
                    datetime(1900, 1, 1, tzinfo=timezone.utc)).total_seconds() / 3600
        
        time_start = datetime(self.year, self.month, 1, tzinfo=timezone.utc).isoformat(" ", "seconds")
        days_in_month = calendar.monthrange(self.year, self.month)[1]
        time_end = datetime(self.year, self.month, days_in_month, 23, 59, 59, tzinfo=timezone.utc).isoformat(" ", "seconds")
        

        with Dataset(rtm_output, "w") as f:
            # ----------
            # Global attributes
            f.setncattr_string("Conventions", "CF-1.9,ACDD-1.3")
            f.setncattr_string("title", "ACCESS RTM output")
            f.setncattr_string("institution", "REMSS")
            f.setncattr_string(
                "history", f"{timestamp} created: {' '.join(sys.argv[1:])}"
            )
            f.setncattr_string("netcdf_version_id", nc_version)
            f.setncattr_string("date_created", timestamp)
            f.setncattr_string("creator_name", "Remote Sensing Systems")
            f.setncattr_string("creator_email", "support@remss.com")
            f.setncattr_string("creator_url", "http://www.remss.com")
            f.setncattr("geospatial_lat_min", np.float32(-90.0))
            f.setncattr("geospatial_lat_max", np.float32(90.0))
            f.setncattr("geospatial_lon_min", np.float32(0.0))
            f.setncattr("geospatial_lon_max", np.float32(360.0))
            f.setncattr_string("time_coverage_start", time_start)
            f.setncattr_string("time_coverage_end", time_end)
            f.setncattr_string("standard_name_vocabulary", "CF Standard Name Table v78")

            # ----------
            # Dimensions
            f.createDimension("time", 1)
            f.createDimension("lat", len(self.lats))
            f.createDimension("lon", len(self.lons))
            f.createDimension("freq", len(self.frequency))
            f.createDimension("eia", len(self.incidence))

            # ----------
            # Coordinate variables
            v = f.createVariable("time", np.int32, ("time",))
            v[:] = data_time
            v.setncattr_string("standard_name", "time")
            v.setncattr_string("axis", "T")
            v.setncattr_string("units", time_units)

            v = f.createVariable("lat", np.float32, ("lat",))
            v[:] = self.lats
            v.setncattr_string("standard_name", "latitude")
            v.setncattr_string("axis", "Y")
            v.setncattr_string("units", "degrees_north")

            v = f.createVariable("lon", np.float32, ("lon",))
            v[:] = self.lons
            v.setncattr_string("standard_name", "longitude")
            v.setncattr_string("axis", "X")
            v.setncattr_string("units", "degrees_east")

            v = f.createVariable("freq", np.float32, ("freq",))
            v[:] = self.frequency
            v.setncattr_string(
                "standard_name", "sensor_band_central_radiation_frequency"
            )
            v.setncattr_string("long_name", "frequency")
            v.setncattr_string("units", "GHz")

            v = f.createVariable("eia", np.float32, ("eia",))
            v[:] = self.incidence
            v.setncattr_string("standard_name", "sensor_zenith_angle")
            v.setncattr_string("long_name", "incidence angle")
            v.setncattr_string("units", "degree")

            # ----------
            # Variables

            v = f.createVariable(
                "tran", np.float32, ("time", "lat", "lon", "eia"), zlib=True
            )
            v[...] = self.transmissivity
            v.setncattr_string("long_name", "atmospheric transmissivity")
            v.setncattr_string("coordinates", "lat lon")

            v = f.createVariable(
                "tb_up", np.float32, ("time", "lat", "lon", "eia"), zlib=True
            )
            v[...] = self.tb_up
            v.setncattr_string("long_name", "upwelling brightness temperature")
            v.setncattr_string("units", "kelvin")
            v.setncattr_string("coordinates", "lat lon")

            v = f.createVariable(
                "tb_down", np.float32, ("time", "lat", "lon", "eia"), zlib=True
            )
            v[...] = self.tb_down
            v.setncattr_string("long_name", "downwelling brightness temperature")
            v.setncattr_string("units", "kelvin")
            v.setncattr_string("coordinates", "lat lon")



def run_rtm_msu(
    era5_data,
    msu_channel: Optional[str] = 'MSU2',
    workers: Optional[int] = None
) -> RtmMonthlyData:
    """Run the RTM on model data."""
    logging.info(f"Running RTM for {msu_channel}")

    channel_info = choose_freq_eia(era5_data, msu_channel)
    freq = channel_info['frequency']
    eia = channel_info['incidence']

    tick = perf_counter_ns()

    # The ERA5 data is organized by lat/lon, but we need to vectorize that down
    # for the RTM and then reshape the output when finished
    shape_4d = era5_data['temperature'].shape
    num_time, num_lat, num_lon, num_levels = shape_4d[0:4]
    num_points = num_time * num_lat * num_lon
    
    
    atmo_results = compute_rtm(
        era5_data['levels'],
        np.reshape(era5_data['temperature'], (num_points, num_levels)),
        np.reshape(era5_data['height'], (num_points, num_levels)),
        np.reshape(era5_data['specific_humidity'], (num_points, num_levels)),
        np.reshape(era5_data['liquid_content'], (num_points, num_levels)),
        np.ravel(era5_data['surface_temperature']),
        np.ravel(era5_data['surface_height']),
        np.ravel(era5_data['surface_dewpoint']),
        np.ravel(era5_data['surface_pressure']),
        eia,
        freq,
        workers,
    )

    tock = perf_counter_ns()
    duration_seconds = (tock - tick) * 1e-9
    logging.info(f"Finished RTM in {duration_seconds:0.2f} s")

    # Now the output values need to be un-vectorized
    num_freq = len(freq)
    
    tran = np.reshape(atmo_results.tran, (num_time, num_lat, num_lon, num_freq))
    tb_up = np.reshape(atmo_results.tb_up, (num_time, num_lat, num_lon, num_freq))
    tb_down = np.reshape(atmo_results.tb_down, (num_time, num_lat, num_lon, num_freq))

    year = era5_data['year']
    month = era5_data['month']

    return RtmMonthlyData(
        year,
        month,
        era5_data['lats'],
        era5_data['lons'],
        freq,
        eia,
        tran,
        tb_up,
        tb_down,
    )


def convert_all_msu(
    era5_surface_input: Path,
    era5_levels_input: Path,
    rtm_output: Path,
    time_indices: Optional[Sequence[int]],
    one_pass: bool,
    msu_channel: Optional[str] = 'MSU2',
    workers: Optional[int] = None,
) -> None:
    """Read the ERA5 profile/surface files and run the RTM and write its output."""
    all_time_indices: Union[Sequence[int], list[int]]
    if time_indices is None:
        all_time_indices = era5.read_time_indices(era5_surface_input, era5_levels_input)
    else:
        all_time_indices = time_indices

    if one_pass:
        era5_data = era5.read_era5_data(
            era5_surface_input, era5_levels_input, all_time_indices
        )
        rtm_data = run_rtm_msu(era5_data, msu_channel, workers)
    else:
        # Read the ERA5 data one hour at a time and process just that much
        hourly_rtm_data = []
        for time_index in all_time_indices:
            era5_data = era5.read_era5_data(
                era5_surface_input, era5_levels_input, [time_index]
            )
            hourly_rtm_data.append(run_rtm_msu(era5_data, msu_channel, workers))

        # Accumulate all the hourly data together
        rtm_data = combine_rtm(hourly_rtm_data)

    logging.info(f"Writing output data: {rtm_output}")
    rtm_data.write_nc(rtm_output)


if __name__ == "__main__":
    main()
