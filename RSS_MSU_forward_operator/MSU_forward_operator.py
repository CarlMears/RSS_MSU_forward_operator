
from process_monthly_msu import run_rtm_msu
from process_monthly_msu import choose_freq_eia

import RSS_surf_emiss.RSS_surf_emiss as RSS_surf_emiss
import importlib.resources  # need this to initialize RSS_surf_emiss

from typing import Dict
from numpy.typing import NDArray
import numpy as np

class MSUForwardOperator:

    def __init__(self, msu_channel):

        self.msu_channel = msu_channel

        #initialize RSS_surf_emiss (ocean emissivity model)
        path_to_data = importlib.resources.path('RSS_surf_emiss').args[0] / 'data'
        RSS_surf_emiss.init(str(path_to_data)+'/')

    def compute_tbs(self, model_data: Dict[str, NDArray[np.float32]]) -> Dict[str, NDArray[np.float32]]:

        workers=1
        rtm_data = run_rtm_msu(model_data, self.msu_channel, workers)

        rtm_shape = rtm_data.transmissivity.shape
        emiss_ocean = np.full((rtm_shape[0], rtm_shape[1], rtm_shape[2], rtm_shape[3], 2), np.nan, dtype=np.float32)  # last dim is pol H/V
        channel_info = choose_freq_eia(model_data, self.msu_channel)
        eia_array = channel_info['incidence']
        freq_array = channel_info['frequency']
    
        for eia_index,eia in enumerate(eia_array):
            freq = freq_array[0]
        
            surface_temp_c = np.ravel(model_data['skin_temperature']-273.15)
            wind = np.ravel(model_data['wind_10m'])
            emiss_this_eia = RSS_surf_emiss.wind_emiss(freq, 
                            eia, 
                            surface_temp_c, 
                            wind, 
                            np.full_like(wind,-999.9, dtype=np.float32))
            
            z = emiss_this_eia.reshape(rtm_shape[0], rtm_shape[1], rtm_shape[2], 2)
            emiss_ocean[:,:,:,eia_index,:] = emiss_this_eia.reshape(rtm_shape[0], rtm_shape[1], rtm_shape[2], 2)

        # combine ocean and land emissivity  Land emissivity is assumed to be 0.9
        emiss_combined = np.full_like(emiss_ocean, np.nan, dtype=np.float32)
        for i in range(rtm_shape[0]):
            for j in range(rtm_shape[3]):
                for k in range(2):
                    emiss_combined[i,:,:,j,k] = emiss_ocean[i,:,:,j,k]*(1 - model_data['land_fraction']) + 0.9*(model_data['land_fraction'])
                    
        # compute TOA brightness temperature for each eia and pol
        tbs_by_eia_pol = np.full_like(emiss_combined, np.nan, dtype=np.float32)
        for eia_index in range(rtm_shape[3]):
            for pol_index in range(2):
                tbs_by_eia_pol[:,:,:,eia_index,pol_index] = (rtm_data.tb_up[:,:,:,eia_index] + \
                    (rtm_data.transmissivity[:,:,:,eia_index] * emiss_combined[:,:,:,eia_index,pol_index] * model_data['skin_temperature'][:, :, :]) +        
                    (1.0 - emiss_combined[:,:,:,eia_index,pol_index]) * rtm_data.transmissivity[:,:,:,eia_index] * rtm_data.tb_down[:,:,:,eia_index])   # downwelling reflected by surface
                            

        # combine the two polarizations into final TB     
        look_angle = channel_info['look_angle']
        cosd_sqr_look_angle = np.cos(np.deg2rad(look_angle))**2
        
        tbs_by_eia = np.full((rtm_shape[0], rtm_shape[1], rtm_shape[2], rtm_shape[3]), np.nan, dtype=np.float32)
        for eia_index in range(rtm_shape[3]):
            if channel_info['polarization'] == 'V':
                tbs_by_eia[:,:,:,eia_index] = tbs_by_eia_pol[:,:,:,eia_index,1]*cosd_sqr_look_angle[eia_index] + tbs_by_eia_pol[:,:,:,eia_index,0]*(1 - cosd_sqr_look_angle[eia_index])  #V-POL at NADIR
            elif channel_info['polarization'] == 'H':
                tbs_by_eia[:,:,:,eia_index] = tbs_by_eia_pol[:,:,:,eia_index,0]*cosd_sqr_look_angle[eia_index] + tbs_by_eia_pol[:,:,:,eia_index,1]*(1 - cosd_sqr_look_angle[eia_index])  #H-POL at NADIR
            else:
                raise ValueError(f"Unsupported polarization: {channel_info['polarization']}")
        
        output_dict = {}
        if self.msu_channel == 'MSU2':
            output_dict['tbs_TMT'] = 0.2*tbs_by_eia[:,:,:,0] + 0.4*tbs_by_eia[:,:,:,1] + 0.4*tbs_by_eia[:,:,:,2]
            output_dict['tbs_TLT'] = 2.0*(tbs_by_eia[:,:,:,2] + tbs_by_eia[:,:,:,3]) - 1.5*(tbs_by_eia[:,:,:,4] + tbs_by_eia[:,:,:,5])
        elif self.msu_channel == 'MSU3':
            output_dict['tbs_TTS'] = 0.2*tbs_by_eia[:,:,:,0] + 0.4*tbs_by_eia[:,:,:,1] + 0.4*tbs_by_eia[:,:,:,2]
        elif self.msu_channel == 'MSU4':
            output_dict['tbs_TLS'] = 0.2*tbs_by_eia[:,:,:,0] + 0.4*tbs_by_eia[:,:,:,1] + 0.4*tbs_by_eia[:,:,:,2]
        else:
            raise ValueError(f'Unsupported msu_channel: {self.msu_channel}')
        
        return output_dict