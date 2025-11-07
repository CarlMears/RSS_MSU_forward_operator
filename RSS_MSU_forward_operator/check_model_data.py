import numpy as np
import matplotlib.pyplot as plt

# This could moved to a config file later
min_max_dict = {
        'year': (0, 3000),
        'month': (1,12),
        'levels': (0,1000.0),
        'lats': (-90.0,90.0),
        'lons': (-180.0,360.0),
        'land_fraction': (0.0,1.0),
        'temperature': (150.0,350.0),
        'specific_humidity': (0.0,0.03),
        'height': (-200.0,100000.0),
        'liquid_content': (0,0.01),
        'surface_pressure': (400.0,1100.0),
        'surface_temperature': (200.0,350.0),
        'surface_dewpoint': (150.0, 350.0),
        'skin_temperature': (200.0, 350.0),
        'surface_height': (-500.0, 9000.0),
        'sea_ice_fraction': (0.0,1.0),
        'wind_10m': (0.0, 100.0),
    }


def check_min_max(arr, min_val, max_val):
    """Check if all elements in arr are within [min_val, max_val]."""
    arr_max = np.nanmax(arr)
    arr_min = np.nanmin(arr)
    return arr_min >= min_val and arr_max <= max_val

def check_model_data(model_data: dict,verbose: bool = False) -> dict:
    """Check if model data meets specified criteria."""
    
    results = {}
    for key, (min_val, max_val) in min_max_dict.items():
        if key in model_data:
            arr = model_data[key]
            is_valid = check_min_max(arr, min_val, max_val)
            if is_valid:
                results[key] = 'In Range'
                if verbose:
                    print(f"{key} is within the range [{min_val}, {max_val}]")
            else:
                results[key] = 'Out of Range'
                if verbose:
                    print(f"Warning: {key} has values outside the range [{min_val}, {max_val}]")
                    print(f"Values found: max {np.nanmax(arr)}, min {np.nanmin(arr)}")

        else:
            results[key] = 'Missing'  # Key not found in model data
            if verbose:
                print(f"Warning: {key} is missing from model data")
    return results
