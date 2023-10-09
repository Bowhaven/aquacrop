import os
os.environ["DEVELOPMENT"] = '1'

from aquacrop import AquaCropModel, Soil, Crop, InitialWaterContent,IrrigationManagement
from aquacrop.utils import prepare_weather, get_filepath
import pandas as pd
import numpy as np

# Initial setup for calibration run:
filepath=get_filepath('cordaba_climate.txt')
weather_data = prepare_weather(filepath)
sandy_loam = Soil(soil_type='SandyLoam')
InitWC = InitialWaterContent(value=['FC'])
irr_mngt = IrrigationManagement(irrigation_method=0)

Rel=60
RedaCCx=77
CCx=96
maize = Crop('MaizeGDDAQTEST', planting_date='05/01',soil_fert_stress=1,
             RelativeBio=Rel/100,Ksccx_in=RedaCCx/CCx,fcdecline_in=1,
             sfertstress=0.72)

# Calibration run
model1 = AquaCropModel(sim_start_time=f'{1986}/05/01',
                      sim_end_time=f'{1986}/08/30',
                      weather_df=weather_data,
                      soil=sandy_loam,
                      crop=maize,
                      irrigation_management=irr_mngt,
                      initial_water_content=InitWC)

model1.run_model(till_termination=True)
model_results = model1.get_crop_growth()

aq=pd.read_table('AquaCropV61Nr02052018\SIMUL\Crop.OUT',skiprows=4, delim_whitespace=True,encoding="latin1")
aq=aq.drop([0])

aq.to_csv('aq_temp.csv')

aq=pd.read_csv('aq_temp.csv')

import matplotlib.pyplot as plt

aq

time=np.array(range(1,len(aq.index)))

# CCx
plt.figure(figsize=(20,10))
plt.rcParams.update({'font.size': 22})
plt.plot(time,aq.loc[time,'CC'],label='Win')
plt.plot(time,model_results.loc[time,'canopy_cover']*100,label='PY')
plt.plot(time,aq.loc[time,'CC']-model_results.loc[time,'canopy_cover']*100,label='Diff')
plt.legend()
plt.show()

# Biomass
plt.figure()
plt.figure(figsize=(20,10))
plt.plot(time,aq.loc[time,'Biomass'],label='Win')
plt.plot(time,model_results.loc[time,'biomass']/100,label='PY')
plt.legend()
plt.show()

# Cumulative Tr
plt.figure()
plt.figure(figsize=(20,10))
plt.plot(time,aq.loc[time,'Tr'],label='Win')
plt.plot(time,model_results.loc[time,'Tr'],label='PY')
plt.legend()
plt.figure()
plt.figure(figsize=(20,10))
plt.plot(time,aq.loc[time,'Tr'],label='Win')
plt.plot(time,model_results.loc[time,'Tr'],label='PY')
plt.legend()
import itertools
plt.figure()
plt.figure(figsize=(20,10))
plt.plot(time,list(itertools.accumulate(aq.loc[time,'Tr'])),label='Win')
plt.plot(time,list(itertools.accumulate(model_results.loc[time,'Tr'])),label='PY')
plt.legend()
plt.show()