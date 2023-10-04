import os
os.environ["DEVELOPMENT"] = '1'

from aquacrop import AquaCropModel, Soil, Crop, InitialWaterContent,IrrigationManagement
from aquacrop.utils import prepare_weather, get_filepath
import pandas as pd
import numpy as np

# Initial setup for calibration run:
filepath=get_filepath('tunis_climate.txt')
weather_data = prepare_weather(filepath)
sandy_loam = Soil(soil_type='SandyLoam')
InitWC = InitialWaterContent(value=['FC'])
irr_mngt = IrrigationManagement(irrigation_method=0)

Rel=60
RedaCCx=77
CCx=96
wheat = Crop('WheatGDDAQTEST', planting_date='11/01',need_calib=1,
             RelativeBio=Rel/100,Ksccx_in=RedaCCx/CCx,fcdecline_in=1,
             sfertstress=0.57)

# Calibration run
model1 = AquaCropModel(sim_start_time=f'{1979}/11/01',
                      sim_end_time=f'{1980}/05/05',
                      weather_df=weather_data,
                      soil=sandy_loam,
                      crop=wheat,
                      irrigation_management=irr_mngt,
                      initial_water_content=InitWC)
# model1._initialize()

# Store soil fertility stress calibrated outputs
# sf_es=model1.crop.sf_es
# Ksexpf_es=model1.crop.Ksexpf_es
# fcdecline_es=model1.crop.fcdecline_es
# Kswp_es=model1.crop.Kswp_es
# Ksccx_es=model1.crop.Ksccx_es
# relbio_es=model1.crop.relbio_es

# # Calculate additional soil fertility stress parameters
# stress=1-0.6
# IrrMethod=0
# loc_=np.argmin(np.abs(sf_es[0:100]-stress))

# Ksccx=Ksccx_es[loc_]
# Ksexpf=Ksexpf_es[loc_]
# Kswp=Kswp_es[loc_]
# fcdecline=fcdecline_es[loc_]

# ccx_=(1-Ksccx)*100
# cgc_=(1-Ksexpf)*100
# dcc_=fcdecline*10000/100
# wp_=(1-Kswp)*100

# wheat = Crop('WheatGDDAQTEST', planting_date='11/01',Ksccx=1-ccx_/100,Ksexpf=1-cgc_/100,Kswp=1-wp_/100,fcdecline=dcc_/100,\
#                  sfertstress=stress,sf_es=sf_es,Ksexpf_es=Ksexpf_es,fcdecline_es=fcdecline_es,Kswp_es=Kswp_es,\
#                 Ksccx_es=Ksccx_es,relbio_es=relbio_es)


# irrMet=IrrigationManagement(irrigation_method=IrrMethod)

# model_os = AquaCropModel(
#             sim_start_time=f"{1979}/11/01",
#             sim_end_time=f"{1980}/05/05",
#             weather_df=weather_data,
#             soil=Soil(soil_type='SandyLoam'),
#             crop=wheat,
#             initial_water_content=InitialWaterContent(value=['FC']),
#             irrigation_management=irrMet
#         )

model1.run_model(till_termination=True)
model_results = model1.get_crop_growth()

aq=pd.read_table('AquaCropV61Nr02052018\SIMUL\Crop.OUT',skiprows=4, delim_whitespace=True,encoding="latin1")
aq=aq.drop([0])

aq.to_csv('aq_temp.csv')

aq=pd.read_csv('aq_temp.csv')

import matplotlib.pyplot as plt

aq

time=np.array(range(1,len(aq.index)))
plt.figure(figsize=(20,10))
plt.rcParams.update({'font.size': 22})
plt.plot(time,aq.loc[time,'CC'],label='Win')
plt.plot(time,model_results.loc[time,'canopy_cover']*100,label='PY')
plt.plot(time,aq.loc[time,'CC']-model_results.loc[time,'canopy_cover']*100,label='Diff')
plt.legend()
# plt.savefig('test_figs\wheat_CCx_irrmethod{0}_stress{1}.png'.format(IrrMethod,stress))
plt.show()

plt.figure()
plt.figure(figsize=(20,10))
plt.plot(time,aq.loc[time,'Biomass'],label='Win')
plt.plot(time,model_results.loc[time,'biomass']/100,label='PY')
plt.legend()
# plt.savefig('test_figs\wheat_Biomass_irrmethod{0}_stress{1}.png'.format(IrrMethod,stress))
plt.show()

plt.figure()
plt.figure(figsize=(20,10))
plt.plot(time,aq.loc[time,'Tr'],label='Win')
plt.plot(time,model_results.loc[time,'Tr'],label='PY')
plt.legend()
plt.show()
