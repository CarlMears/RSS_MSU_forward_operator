"""Read ERA5 daily surface and profile data."""

import logging
from collections.abc import Sequence
from pathlib import Path
from typing import NamedTuple, Optional, Union

import numpy as np
from netCDF4 import Dataset
from numpy.typing import NDArray

class Era5MonthlyData(NamedTuple):
    """The daily data for ERA5 surface and levels data."""

    year: int
    month: int

    # Pressure levels in hPa, with shape (num_levels, ). They should be in
    # descending order (e.g., 1000 to 10).
    levels: NDArray[np.float32]

    # Latitude in degrees North, dimensioned as (num_lats, ). They should be in
    # ascending order (e.g., -90 to 90).
    lats: NDArray[np.float32]

    # Longitude in degrees East, dimensioned as (num_lons, ). They should be in
    # ascending order (e.g., -180 to 180).
    lons: NDArray[np.float32]

    # 2D land fraction, dimensioned as (lats, lons)
    land_fraction: NDArray[np.int32]

    # Profile air temperature in kelvin, dimensioned as (time, lats, lons, levels)
    temperature: NDArray[np.float32]

    # Profile specific humidity in kg/kg, dimensioned as (time, lats, lons,
    # levels)
    specific_humidity: NDArray[np.float32]

    # Geopotential height profile in meters, dimensioned as (time, lats, lons, levels)
    height: NDArray[np.float32]

    # Profile specific liquid water content (from clouds) in kg/kg, dimensioned
    # as (time, lats, lons, levels)
    liquid_content: NDArray[np.float32]

    # Surface pressure in hPa, dimensioned as (time, lats, lons)
    surface_pressure: NDArray[np.float32]

    # 2-meter air temperature in kelvin, dimensioned as (time, lats, lons)
    surface_temperature: NDArray[np.float32]

    # 2-meter dewpoint in kelvin, dimensioned as (time, lats, lons)
    surface_dewpoint: NDArray[np.float32]

    # surface skin temperature in kelvin, dimensioned as (time, lats, lons)
    skin_temperature: NDArray[np.float32]

    # Geopotential height at the surface in meters, dimensioned as (time, lats, lons)
    surface_height: NDArray[np.float32]

    # sea ice fraction, dimensioned as (time, lats, lons)
    sea_ice_fraction: NDArray[np.float32]

    # 10-meter wind speed, dimensioned as (time, lats, lons)
    wind_10m: NDArray[np.float32]

def read_era5_data_monthly(
    input_files: dict,
) -> dict:
    
    """Read ERA5 surface/levels files.
    """

    # The non-coordinate variables are all stored as packed integers and
    # automatically unpacked to float64. To reduce peak memory usage, each one
    # is converted to a float32 array.
    with Dataset(input_files['land_frac'], "r") as f:

        lats = f["latitude"][:]
        lons = f["longitude"][:]
        land_fraction = f["lsm"][:,:].astype(np.float32)

    with Dataset(input_files['orography'], "r") as f:
        surface_height = f["z"][:,:].astype(np.float32)
        # Note: orography is in geopotential, need to convert to meters later
        

    with Dataset(input_files['surf_pressure'], "r") as f:
        surface_pressure = f["sp"][:,:].astype(np.float32)

    with Dataset(input_files['sea_ice'], "r") as f:
        sea_ice_fraction = f["siconc"][:,:].astype(np.float32)

    with Dataset(input_files['t2m'], "r") as f:
        surface_temperature = f["t2m"][:,:].astype(np.float32)

    with Dataset(input_files['tdew'], "r") as f:
        surface_dewpoint = f["d2m"][:,:].astype(np.float32)

    with Dataset(input_files['tskin'], "r") as f:
        skin_temperature = f["skt"][:,:].astype(np.float32)

    with Dataset(input_files['w10m'], "r") as f:
        wind_10m = f["si10"][:,:].astype(np.float32)

    with Dataset(input_files['cld'], "r") as f:
        levels = f["level"][:].astype(np.float32)
        liquid_content = f["clwc"][:,:,:].astype(np.float32)

    with Dataset(input_files['geopot'], "r") as f:
        height = f["z"][:,:,:].astype(np.float32)
        # Convert geopotential to geopotential height later
        

    with Dataset(input_files['q'], "r") as f:
        specific_humidity = f["q"][:,:,:].astype(np.float32)

    with Dataset(input_files['t'], "r") as f:
        temperature = f["t"][:,:,:].astype(np.float32)

    logging.info(f"Post-processing ERA5 data ({len(levels)} pressure levels)")

    # By default, netCDF4 returns masked arrays for all the variables above.
    # However, there shouldn't be any values that are actually masked in the
    # ERA5 data. Check that assumption and then convert everything to "vanilla"
    # ndarrays.


    if any(
        np.ma.count_masked(a) > 0
        for a in (
            lats,
            lons,
            land_fraction,
            surface_pressure,
            surface_temperature,
            surface_dewpoint,
            surface_height,
            skin_temperature,
            #sea_ice_fraction, # sea ice can have masked values
            wind_10m,
            levels,
            temperature,
            specific_humidity,
            height,
            liquid_content,
        )
    ):
        raise Exception("Masked input values detected")
    lats = np.ma.getdata(lats)
    lons = np.ma.getdata(lons)
    land_fraction = np.ma.getdata(land_fraction)
    surface_pressure = np.ma.getdata(surface_pressure)
    surface_temperature = np.ma.getdata(surface_temperature)
    surface_dewpoint = np.ma.getdata(surface_dewpoint)
    surface_height = np.ma.getdata(surface_height)
    skin_temperature = np.ma.getdata(skin_temperature)
    wind_10m = np.ma.getdata(wind_10m)
    levels = np.ma.getdata(levels)
    temperature = np.ma.getdata(temperature)
    specific_humidity = np.ma.getdata(specific_humidity)
    height = np.ma.getdata(height)
    liquid_content = np.ma.getdata(liquid_content)

    # sea ice has masked values where there is no ice, fill those with 0
    sea_ice_fraction = sea_ice_fraction.filled(fill_value=0)

    # The reciprocal of the standard gravity, in units of s^2 / m
    # https://en.wikipedia.org/wiki/Standard_gravity
    INV_STANDARD_GRAVITY = 1 / 9.80665

    # Convert geopotential to geopotential height
    # (https://apps.ecmwf.int/codes/grib/param-db?id=129)
    height *= INV_STANDARD_GRAVITY
    surface_height *= INV_STANDARD_GRAVITY

    # Convert surface pressure from Pa to hPa
    surface_pressure *= 1e-2

    # The 4d arrays should be reordered from (time, levels, lat, lon) to (time,
    # lat, lon, levels)
    temperature = np.moveaxis(temperature, 1, -1)
    specific_humidity = np.moveaxis(specific_humidity, 1, -1)
    height = np.moveaxis(height, 1, -1)
    liquid_content = np.moveaxis(liquid_content, 1, -1)

    # The latitudes need to be adjusted. In the ERA5 files, the latitudes are in
    # *descending* order from 90 to -90, and the longitudes are in ascending
    # order from 0 to 360. Leave the longitudes alone, but flip the latitudes so
    # they go from -90 to 90.
    lats = -lats

    temperature = np.flip(temperature, 1).copy()
    specific_humidity = np.flip(specific_humidity, 1).copy()
    height = np.flip(height, 1).copy()
    liquid_content = np.flip(liquid_content, 1).copy()
    surface_pressure = np.flip(surface_pressure, 1).copy()
    surface_temperature = np.flip(surface_temperature, 1).copy()
    surface_dewpoint = np.flip(surface_dewpoint, 1).copy()
    surface_height = np.flip(surface_height, 1).copy()
    skin_temperature = np.flip(skin_temperature, 1).copy()
    wind_10m = np.flip(wind_10m, 1).copy()
    sea_ice_fraction = np.flip(sea_ice_fraction, 1).copy()
    land_fraction = np.flip(land_fraction, 1).copy()

    # the pressure levels are upside down, from 1,000 hPa to 1 hPa, need to reverse them
    levels = levels[::-1].copy()
    temperature = np.flip(temperature,3).copy()
    specific_humidity = np.flip(specific_humidity,3).copy()
    height = np.flip(height,3).copy()
    liquid_content = np.flip(liquid_content,3).copy()
    liquid_content[liquid_content < 0] = 0.0  # ensure no negative cloud liquid water content

    year = input_files['year']
    month = input_files['month']

    return_dict = {
        'year': year,
        'month': month,
        'levels': levels,
        'lats': lats,
        'lons': lons,
        'land_fraction': land_fraction,
        'temperature': temperature,
        'specific_humidity': specific_humidity,
        'height': height,
        'liquid_content': liquid_content,
        'surface_pressure': surface_pressure,
        'surface_temperature': surface_temperature,
        'surface_dewpoint': surface_dewpoint,
        'skin_temperature': skin_temperature,
        'surface_height': surface_height,
        'sea_ice_fraction': sea_ice_fraction,
        'wind_10m': wind_10m
    }

    return return_dict

def era5_monthly_files(year_to_do : int, month_to_do : int, path_to_rtm : Path) -> dict:

    print(f'Processing month {month_to_do} of year {year_to_do}')
    land_frac_file = path_to_rtm / '2D' / f'{year_to_do:04d}' / f'era5_land_frac_2D_360_181_{year_to_do:04d}_{month_to_do:02d}.nc'
    orography_file = path_to_rtm / '2D'  / f'{year_to_do:04d}' / f'era5_GEO_2D_360_181_{year_to_do:04d}_{month_to_do:02d}.nc'
    surf_pressure_file = path_to_rtm / '2D'  / f'{year_to_do:04d}' / f'era5_PS_2D_360_181_{year_to_do:04d}_{month_to_do:02d}.nc'
    sea_ice_file = path_to_rtm / '2D'  / f'{year_to_do:04d}' / f'era5_SeaIce_2D_360_181_{year_to_do:04d}_{month_to_do:02d}.nc'
    t2m_file = path_to_rtm / '2D'  / f'{year_to_do:04d}' / f'era5_T2m_2D_360_181_{year_to_do:04d}_{month_to_do:02d}.nc'
    tdew_file = path_to_rtm / '2D'  / f'{year_to_do:04d}' / f'era5_TDew_2D_360_181_{year_to_do:04d}_{month_to_do:02d}.nc'
    tskin_file = path_to_rtm / '2D'  / f'{year_to_do:04d}' / f'era5_TSkin_2D_360_181_{year_to_do:04d}_{month_to_do:02d}.nc'
    w10m_file = path_to_rtm / '2D'  / f'{year_to_do:04d}' / f'era5_W10_2D_360_181_{year_to_do:04d}_{month_to_do:02d}.nc'

    cld_file = path_to_rtm / '3D' / f'{year_to_do:04d}' / f'era5_CLD_3D_360_181_{year_to_do:04d}_{month_to_do:02d}.nc'
    geopot_file = path_to_rtm / '3D' / f'{year_to_do:04d}' / f'era5_G_3D_360_181_{year_to_do:04d}_{month_to_do:02d}.nc'
    q_file = path_to_rtm / '3D' / f'{year_to_do:04d}' / f'era5_Q_3D_360_181_{year_to_do:04d}_{month_to_do:02d}.nc'
    t_file = path_to_rtm / '3D' / f'{year_to_do:04d}' / f'era5_T_3D_360_181_{year_to_do:04d}_{month_to_do:02d}.nc'

    input_files = {
        'year': year_to_do,
        'month': month_to_do,
        'land_frac': land_frac_file,
        'orography': orography_file,
        'surf_pressure': surf_pressure_file,
        'sea_ice': sea_ice_file,
        't2m': t2m_file,
        'tdew': tdew_file,
        'tskin': tskin_file,
        'w10m': w10m_file,
        'cld': cld_file,
        'geopot': geopot_file,
        'q': q_file,
        't': t_file
    }

    return input_files



class Era5DailyData(NamedTuple):
    """The daily data for ERA5 surface and levels data."""

    # Pressure levels in hPa, with shape (num_levels, ). They should be in
    # descending order (e.g., 1000 to 10).
    levels: NDArray[np.float32]

    # Latitude in degrees North, dimensioned as (num_lats, ). They should be in
    # ascending order (e.g., -90 to 90).
    lats: NDArray[np.float32]

    # Longitude in degrees East, dimensioned as (num_lons, ). They should be in
    # ascending order (e.g., -180 to 180).
    lons: NDArray[np.float32]

    # Hours since 1900-01-01, dimensioned as (num_time, ).
    time: NDArray[np.int32]

    # Profile air temperature in kelvin, dimensioned as (time, lats, lons, levels)
    temperature: NDArray[np.float32]

    # Profile specific humidity in kg/kg, dimensioned as (time, lats, lons,
    # levels)
    specific_humidity: NDArray[np.float32]

    # Geopotential height profile in meters, dimensioned as (time, lats, lons, levels)
    height: NDArray[np.float32]

    # Profile specific liquid water content (from clouds) in kg/kg, dimensioned
    # as (time, lats, lons, levels)
    liquid_content: NDArray[np.float32]

    # Surface pressure in hPa, dimensioned as (time, lats, lons)
    surface_pressure: NDArray[np.float32]

    # 2-meter air temperature in kelvin, dimensioned as (time, lats, lons)
    surface_temperature: NDArray[np.float32]

    # 2-meter dewpoint in kelvin, dimensioned as (time, lats, lons)
    surface_dewpoint: NDArray[np.float32]

    # Geopotential height at the surface in meters, dimensioned as (time, lats, lons)
    surface_height: NDArray[np.float32]

    # Total column water vapor in kg/m^2, dimensioned as (time, lats, lons)
    columnar_water_vapor: NDArray[np.float32]

    # Total column cloud liquid water in kg/m^2, dimensioned as (time, lats, lons)
    columnar_cloud_liquid: NDArray[np.float32]


def buck_vap(temperature: NDArray[np.float32]) -> NDArray[np.float32]:
    """Buck equation.

    Use the Buck equation to convert temperature in kelvin into water vapor
    saturation pressure in hPa. The equation is from [1], which cites Buck 1996.

    To convert to water vapor partial pressure, multiply the result by the
    relative humidity.

    [1] https://en.wikipedia.org/wiki/Arden_Buck_equation
    """
    # Temperature in degrees Celsius
    temp_c: NDArray[np.float32] = temperature - 273.15
    return 6.1121 * np.exp((18.678 - temp_c / 234.5) * (temp_c / (257.14 + temp_c)))


def read_time_indices(surface_file: Path, levels_file: Path) -> list[int]:
    """Return the available time indices in the two ERA5 files."""
    with Dataset(surface_file, "r") as f:
        surface_times = f["time"][:]
    with Dataset(levels_file, "r") as f:
        level_times = f["time"][:]

    # The "time" coordinate variable is an integer so an exact comparison works
    if not np.array_equal(surface_times, level_times):
        raise Exception("ERA5 surface/level files have mismatched times")

    # The actual values for the time variable don't matter at this point, only
    # the available indices, which go from 0 to N-1.
    return list(range(len(surface_times)))


def read_era5_data(
    surface_file: Path,
    levels_file: Path,
    time_subset: Optional[Sequence[int]] = None,
) -> Era5DailyData:
    """Read the pair of ERA5 surface/levels files.

    Optionally, a subset of the time values can be read.
    """
    logging.info(f"Reading surface data: {surface_file}")
    if time_subset is not None:
        logging.info(f"Subsetting hour indices to: {time_subset}")

    times: Union[slice, Sequence[int]]
    if time_subset is None:
        times = slice(None)
    else:
        times = time_subset

    # The non-coordinate variables are all stored as packed integers and
    # automatically unpacked to float64. To reduce peak memory usage, each one
    # is converted to a float32 array.
    with Dataset(surface_file, "r") as f:
        lats = f["latitude"][:]
        lons = f["longitude"][:]
        time = f["time"][times]
        surface_pressure = f["sp"][times, :, :].astype(np.float32)
        surface_temperature = f["t2m"][times, :, :].astype(np.float32)
        surface_dewpoint = f["d2m"][times, :, :].astype(np.float32)
        surface_height = f["z"][times, :, :].astype(np.float32)
        columnar_water_vapor = f["tcwv"][times, :, :].astype(np.float32)
        columnar_cloud_liquid = f["tclw"][times, :, :].astype(np.float32)

    logging.info(f"Reading profiles data: {levels_file}")
    with Dataset(levels_file, "r") as f:
        # TODO: Probably the lats/lons/time coordinate variables should be
        # checked to ensure the levels file matches the surface file...but for
        # now we'll just be really trusting
        levels = f["level"][:].astype(np.float32)
        temperature = f["t"][times, :, :, :].astype(np.float32)
        specific_humidity = f["q"][times, :, :, :].astype(np.float32)
        height = f["z"][times, :, :, :].astype(np.float32)
        liquid_content = f["clwc"][times, :, :, :].astype(np.float32)

    logging.info(f"Post-processing ERA5 data ({len(levels)} pressure levels)")

    # By default, netCDF4 returns masked arrays for all the variables above.
    # However, there shouldn't be any values that are actually masked in the
    # ERA5 data. Check that assumption and then convert everything to "vanilla"
    # ndarrays.
    if any(
        np.ma.count_masked(a) > 0
        for a in (
            lats,
            lons,
            time,
            surface_pressure,
            surface_temperature,
            surface_dewpoint,
            surface_height,
            columnar_water_vapor,
            columnar_cloud_liquid,
            levels,
            temperature,
            specific_humidity,
            height,
            liquid_content,
        )
    ):
        raise Exception("Masked input values detected")
    lats = np.ma.getdata(lats)
    lons = np.ma.getdata(lons)
    time = np.ma.getdata(time)
    surface_pressure = np.ma.getdata(surface_pressure)
    surface_temperature = np.ma.getdata(surface_temperature)
    surface_dewpoint = np.ma.getdata(surface_dewpoint)
    surface_height = np.ma.getdata(surface_height)
    columnar_water_vapor = np.ma.getdata(columnar_water_vapor)
    columnar_cloud_liquid = np.ma.getdata(columnar_cloud_liquid)
    levels = np.ma.getdata(levels)
    temperature = np.ma.getdata(temperature)
    specific_humidity = np.ma.getdata(specific_humidity)
    height = np.ma.getdata(height)
    liquid_content = np.ma.getdata(liquid_content)

    # The reciprocal of the standard gravity, in units of s^2 / m
    # https://en.wikipedia.org/wiki/Standard_gravity
    INV_STANDARD_GRAVITY = 1 / 9.80665

    # Convert geopotential to geopotential height
    # (https://apps.ecmwf.int/codes/grib/param-db?id=129)
    height *= INV_STANDARD_GRAVITY
    surface_height *= INV_STANDARD_GRAVITY

    # Convert surface pressure from Pa to hPa
    surface_pressure *= 1e-2

    # The 4d arrays should be reordered from (time, levels, lat, lon) to (time,
    # lat, lon, levels)
    temperature = np.moveaxis(temperature, 1, -1)
    specific_humidity = np.moveaxis(specific_humidity, 1, -1)
    height = np.moveaxis(height, 1, -1)
    liquid_content = np.moveaxis(liquid_content, 1, -1)

    # The latitudes need to be adjusted. In the ERA5 files, the latitudes are in
    # *descending* order from 90 to -90, and the longitudes are in ascending
    # order from 0 to 360. Leave the longitudes alone, but flip the latitudes so
    # they go from -90 to 90.
    lats = -lats

    temperature = np.flip(temperature, 1)
    specific_humidity = np.flip(specific_humidity, 1)
    height = np.flip(height, 1)
    liquid_content = np.flip(liquid_content, 1)
    surface_pressure = np.flip(surface_pressure, 1)
    surface_temperature = np.flip(surface_temperature, 1)
    surface_dewpoint = np.flip(surface_dewpoint, 1)
    surface_height = np.flip(surface_height, 1)
    columnar_water_vapor = np.flip(columnar_water_vapor, 1)
    columnar_cloud_liquid = np.flip(columnar_cloud_liquid, 1)

    return Era5DailyData(
        levels,
        lats,
        lons,
        time,
        temperature,
        specific_humidity,
        height,
        liquid_content,
        surface_pressure,
        surface_temperature,
        surface_dewpoint,
        surface_height,
        columnar_water_vapor,
        columnar_cloud_liquid,
    )
