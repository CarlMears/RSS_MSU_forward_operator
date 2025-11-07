import geomod10.geomod10 as geomod10
import importlib.resources
import numpy as np


path_to_data = importlib.resources.path('geomod10').args[0] / 'data'
print(dir(geomod10))
print(path_to_data)



freq = 19.0
tht = 53.0
sst = np.array([10.0, 10.0, 10.0, 10.0], dtype=np.float32)
wind = np.array([5.0, 7.0, 9.0, 11.0], dtype=np.float32)
phir = np.array([-999.0, -999.0, -999.0, -999.0], dtype=np.float32)

emiss = geomod10.wind_emiss(freq, tht, sst, wind, phir)
print(emiss)
geomod10.init(str(path_to_data)+'/')
print()

emiss = geomod10.wind_emiss(freq, tht, sst, wind, phir)
print(emiss)

expected_emiss = np.array(
    [
        [0.58449423, 0.28792182],
        [0.5853587, 0.29340607],
        [0.588001, 0.30069447],
        [0.59222317, 0.31007487],
    ],
    dtype=np.float32,
)
try:
    assert np.array_equal(emiss, expected_emiss)
    print("Test passed!")
except AssertionError:
    print("Test failed!")

