import numpy as np
import pandas as pd

from ..entities.modelConstants import ModelConstants
from typing import TYPE_CHECKING

from .calculate_HIGC import calculate_HIGC
from .calculate_HI_linear import calculate_HI_linear

if TYPE_CHECKING:
    # Important: classes are only imported when types are checked, not in production.
    from aquacrop.entities.crop import Crop
    from pandas import DatetimeIndex, DataFrame


def compute_crop_calendar(
    crop: "Crop",
    clock_struct_planting_dates: "DatetimeIndex",
    clock_struct_simulation_start_date: str,
    clock_struct_time_span: "DatetimeIndex",
    weather_df: "DataFrame",
    ParamStruct:"ParamStruct"):
    """
    Function to compute additional parameters needed to define crop phenological calendar

    <a href="https://www.fao.org/3/BR248E/br248e.pdf#page=28" target="_blank">Reference Manual</a> (pg. 19-20)


    Arguments:

        crop (Crop):  Crop object containing crop paramaters

        clock_struct_planting_dates (DatetimeIndex):  list of planting dates

        clock_struct_simulation_start_date (str):  sim start date

        clock_struct_time_span (DatetimeIndex):  all dates between sim start and end dates

        weather_df (DataFrame):  weather data for simulation period


    Returns:

        crop (Crop): updated Crop object



    """

    if len(clock_struct_planting_dates) == 0:
        plant_year = pd.DatetimeIndex([clock_struct_simulation_start_date]).year[0]
        if (
            pd.to_datetime(str(plant_year) + "/" + crop.planting_date)
            < clock_struct_simulation_start_date
        ):
            pl_date = str(plant_year + 1) + "/" + crop.planting_date
        else:
            pl_date = str(plant_year) + "/" + crop.planting_date
    else:
        pl_date = clock_struct_planting_dates[0]

    # Define crop calendar mode
    Mode = crop.CalendarType

    # Calculate variables %%
    if Mode == 1:  # Growth in calendar days

        # Time from sowing to end of vegatative growth period
        if crop.Determinant == 1:
            crop.CanopyDevEndCD = round(crop.HIstartCD + (crop.FloweringCD / 2))
        else:
            crop.CanopyDevEndCD = crop.SenescenceCD

        # Time from sowing to 10% canopy cover (non-stressed conditions)
        crop.Canopy10PctCD = round(
            crop.EmergenceCD + (np.log(0.1 / crop.CC0) / crop.CGC_CD)
        )

        # Time from sowing to maximum canopy cover (non-stressed conditions)
        crop.MaxCanopyCD = round(
            crop.EmergenceCD
            + (
                np.log(
                    (0.25 * crop.CCx * crop.CCx / crop.CC0)
                    / (crop.CCx - (0.98 * crop.CCx))
                )
                / crop.CGC_CD
            )
        )

        # Time from sowing to end of yield_ formation
        crop.HIendCD = crop.HIstartCD + crop.YldFormCD

        # Duplicate calendar values (needed to minimise if
        # statements when switching between gdd and CD runs)
        crop.Emergence = crop.EmergenceCD
        crop.Canopy10Pct = crop.Canopy10PctCD
        crop.MaxRooting = crop.MaxRootingCD
        crop.Senescence = crop.SenescenceCD
        crop.Maturity = crop.MaturityCD
        crop.MaxCanopy = crop.MaxCanopyCD
        crop.CanopyDevEnd = crop.CanopyDevEndCD
        crop.HIstart = crop.HIstartCD
        crop.HIend = crop.HIendCD
        crop.YldForm = crop.YldFormCD
        if crop.CropType == 3:
            crop.FloweringEndCD = crop.HIstartCD + crop.FloweringCD
            # crop.FloweringEndCD = crop.FloweringEnd
            # crop.FloweringCD = crop.Flowering
        else:
            crop.FloweringEnd = ModelConstants.NO_VALUE
            crop.FloweringEndCD = ModelConstants.NO_VALUE
            crop.FloweringCD = ModelConstants.NO_VALUE

        # Check if converting crop calendar to gdd mode
        if crop.SwitchGDD == 1:
            #             # Extract weather data for first growing season that crop is planted
            #             for i,n in enumerate(ParamStruct.CropChoices):
            #                 if n == crop.Name:
            #                     idx = i
            #                     break
            #                 else:
            #                     idx = -1
            #             assert idx > -1

            date_range = pd.date_range(pl_date, clock_struct_time_span[-1])
            weather_df = weather_df.copy()
            weather_df.index = weather_df.Date
            weather_df = weather_df.loc[date_range]
            temp_min = weather_df.MinTemp
            temp_max = weather_df.MaxTemp

            # Calculate gdd's
            if crop.GDDmethod == 1:

                Tmean = (temp_max + temp_min) / 2
                Tmean = Tmean.clip(lower=crop.Tbase, upper=crop.Tupp)
                gdd = Tmean - crop.Tbase

            elif crop.GDDmethod == 2:

                temp_max = temp_max.clip(lower=crop.Tbase, upper=crop.Tupp)
                temp_min = temp_min.clip(lower=crop.Tbase, upper=crop.Tupp)
                Tmean = (temp_max + temp_min) / 2
                gdd = Tmean - crop.Tbase

            elif crop.GDDmethod == 3:

                temp_max = temp_max.clip(lower=crop.Tbase, upper=crop.Tupp)
                temp_min = temp_min.clip(upper=crop.Tupp)
                Tmean = (temp_max + temp_min) / 2
                Tmean = Tmean.clip(lower=crop.Tbase)
                gdd = Tmean - crop.Tbase

            gdd_cum = np.cumsum(gdd)
            # Find gdd equivalent for each crop calendar variable
            # 1. gdd's from sowing to emergence
            crop.Emergence = gdd_cum.iloc[int(crop.EmergenceCD)]
            # 2. gdd's from sowing to 10# canopy cover
            crop.Canopy10Pct = gdd_cum.iloc[int(crop.Canopy10PctCD)]
            # 3. gdd's from sowing to maximum rooting
            crop.MaxRooting = gdd_cum.iloc[int(crop.MaxRootingCD)]
            # 4. gdd's from sowing to maximum canopy cover
            crop.MaxCanopy = gdd_cum.iloc[int(crop.MaxCanopyCD)]
            # 5. gdd's from sowing to end of vegetative growth
            crop.CanopyDevEnd = gdd_cum.iloc[int(crop.CanopyDevEndCD)]
            # 6. gdd's from sowing to senescence
            crop.Senescence = gdd_cum.iloc[int(crop.SenescenceCD)]
            # 7. gdd's from sowing to maturity
            crop.Maturity = gdd_cum.iloc[int(crop.MaturityCD)]
            # 8. gdd's from sowing to start of yield_ formation
            crop.HIstart = gdd_cum.iloc[int(crop.HIstartCD)]
            # 9. gdd's from sowing to start of yield_ formation
            crop.HIend = gdd_cum.iloc[int(crop.HIendCD)]
            # 10. Duration of yield_ formation (gdd's)
            crop.YldForm = crop.HIend - crop.HIstart

            # 11. Duration of flowering (gdd's) - (fruit/grain crops only)
            if crop.CropType == 3:
                # gdd's from sowing to end of flowering
                crop.FloweringEnd = gdd_cum.iloc[int(crop.FloweringEndCD)]
                # Duration of flowering (gdd's)
                crop.Flowering = crop.FloweringEnd - crop.HIstart

            # Convert CGC to gdd mode
            # crop.CGC_CD = crop.CGC
            crop.CGC = (
                np.log(
                    (((0.98 * crop.CCx) - crop.CCx) * crop.CC0)
                    / (-0.25 * (crop.CCx**2))
                )
            ) / (-(crop.MaxCanopy - crop.Emergence))

            # Convert CDC to gdd mode
            # crop.CDC_CD = crop.CDC
            tCD = crop.MaturityCD - crop.SenescenceCD
            if tCD <= 0:
                tCD = 1

            CCi = crop.CCx * (1 - 0.05 * (np.exp((crop.CDC_CD / crop.CCx) * tCD) - 1))
            if CCi < 0:
                CCi = 0

            tGDD = crop.Maturity - crop.Senescence
            if tGDD <= 0:
                tGDD = 5

            crop.CDC = (crop.CCx / tGDD) * np.log(1 + ((1 - CCi / crop.CCx) / 0.05))
            # Set calendar type to gdd mode
            crop.CalendarType = 2

        else:
            crop.CDC = crop.CDC_CD
            crop.CGC = crop.CGC_CD

        # print(crop.__dict__)
    elif Mode == 2:
        # Growth in growing degree days
        # Time from sowing to end of vegatative growth period
        if crop.Determinant == 1:
            crop.CanopyDevEnd = round(crop.HIstart + (crop.Flowering / 2))
        else:
            crop.CanopyDevEnd = crop.Senescence

        # Time from sowing to 10# canopy cover (non-stressed conditions)
        crop.Canopy10Pct = round(crop.Emergence + (np.log(0.1 / crop.CC0) / crop.CGC))

        # Time from sowing to maximum canopy cover (non-stressed conditions)
        crop.MaxCanopy = round(crop.Emergence+(np.log((0.25*crop.CCx*crop.Ksccx*crop.CCx*crop.Ksccx/crop.CC0)
                                                                    /(crop.CCx*crop.Ksccx-(0.98*crop.CCx*crop.Ksccx)))/crop.CGC/crop.Ksexpf))

        # Time from sowing to end of yield_ formation
        crop.HIend = crop.HIstart + crop.YldForm

        # Time from sowing to end of flowering (if fruit/grain crop)
        if crop.CropType == 3:
            crop.FloweringEnd = crop.HIstart + crop.Flowering

        # Extract weather data for first growing season that crop is planted
        #         for i,n in enumerate(ParamStruct.CropChoices):
        #             if n == crop.Name:
        #                 idx = i
        #                 break
        #             else:
        #                 idx = -1
        #         assert idx> -1
        date_range = pd.date_range(pl_date, clock_struct_time_span[-1])
        weather_df = weather_df.copy()
        weather_df.index = weather_df.Date

        weather_df = weather_df.loc[date_range]
        temp_min = weather_df.MinTemp
        temp_max = weather_df.MaxTemp

        # Calculate gdd's
        if crop.GDDmethod == 1:

            Tmean = (temp_max + temp_min) / 2
            Tmean = Tmean.clip(lower=crop.Tbase, upper=crop.Tupp)
            gdd = Tmean - crop.Tbase

        elif crop.GDDmethod == 2:

            temp_max = temp_max.clip(lower=crop.Tbase, upper=crop.Tupp)
            temp_min = temp_min.clip(lower=crop.Tbase, upper=crop.Tupp)
            Tmean = (temp_max + temp_min) / 2
            gdd = Tmean - crop.Tbase

        elif crop.GDDmethod == 3:

            temp_max = temp_max.clip(lower=crop.Tbase, upper=crop.Tupp)
            temp_min = temp_min.clip(upper=crop.Tupp)
            Tmean = (temp_max + temp_min) / 2
            Tmean = Tmean.clip(lower=crop.Tbase)
            gdd = Tmean - crop.Tbase

        gdd_cum = np.cumsum(gdd).reset_index(drop=True)

        assert (
            gdd_cum.values[-1] > crop.Maturity
        ), f"not enough growing degree days in simulation ({gdd_cum.values[-1]}) to reach maturity ({crop.Maturity})"

        crop.MaturityCD = (gdd_cum > crop.Maturity).idxmax() + 1

        assert crop.MaturityCD < 365, "crop will take longer than 1 year to mature"

        # 1. gdd's from sowing to maximum canopy cover
        crop.MaxCanopyCD = (gdd_cum > crop.MaxCanopy).idxmax() + 1
        # 2. gdd's from sowing to end of vegetative growth
        crop.CanopyDevEndCD = (gdd_cum > crop.CanopyDevEnd).idxmax() + 1
        # 3. Calendar days from sowing to start of yield_ formation
        crop.HIstartCD = (gdd_cum > crop.HIstart).idxmax() + 1
        # 4. Calendar days from sowing to end of yield_ formation
        crop.HIendCD = (gdd_cum > crop.HIend).idxmax() + 1
        # 5. Duration of yield_ formation in calendar days
        crop.YldFormCD = crop.HIendCD - crop.HIstartCD
        
        crop.SenescenceCD = (gdd_cum>crop.Senescence).idxmax()+1
        crop.EmergenceCD = (gdd_cum>crop.Emergence).idxmax()+1

        #look unnecessary, but cannot get the same results with AquaCrop-win without it
        if crop.Ksccx<1 or crop.Ksexpf<1:
            if crop.CGC_CD==-1:
                crop.CGC_CD=crop.MaxCanopy/crop.MaxCanopyCD*crop.CGC
                
            crop.MaxCanopyCD = round(crop.EmergenceCD+(np.log((0.25*crop.CCx*crop.Ksccx*crop.CCx*crop.Ksccx/crop.CC0)
                                                                        /(crop.CCx*crop.Ksccx-(0.98*crop.CCx*crop.Ksccx)))/crop.CGC_CD/crop.Ksexpf))

            if crop.MaxCanopyCD>crop.CanopyDevEndCD:
                while crop.MaxCanopyCD>crop.CanopyDevEndCD and crop.Ksexpf<1:
                    crop.Ksexpf=crop.Ksexpf+0.01
                    crop.MaxCanopyCD = round(crop.EmergenceCD+(np.log((0.25*crop.CCx*crop.Ksccx*crop.CCx*crop.Ksccx/crop.CC0)
                                                                    /(crop.CCx*crop.Ksccx-(0.98*crop.CCx*crop.Ksccx)))/crop.CGC_CD/crop.Ksexpf))
                while crop.MaxCanopyCD>crop.CanopyDevEndCD and crop.CCx*crop.Ksccx>0.1 and crop.Ksccx>0.5:
                    crop.Ksccx=crop.Ksccx-0.01
                    crop.MaxCanopyCD = round(crop.EmergenceCD+(np.log((0.25*crop.CCx*crop.Ksccx*crop.CCx*crop.Ksccx/crop.CC0)
                                                                    /(crop.CCx*crop.Ksccx-(0.98*crop.CCx*crop.Ksccx)))/crop.CGC_CD/crop.Ksexpf))
            crop.MaxCanopy=gdd_cum.values[crop.MaxCanopyCD-1]  
        
        if crop.CropType == 3:
            # 1. Calendar days from sowing to end of flowering
            FloweringEnd = (gdd_cum > crop.FloweringEnd).idxmax() + 1
            # 2. Duration of flowering in calendar days
            crop.FloweringCD = FloweringEnd - crop.HIstartCD
        else:
            crop.FloweringCD = ModelConstants.NO_VALUE

        #calculate the normalized Tr for soil fertility stress, it is a theritically one, could be derived without simulation 
        CCx=crop.CCx#*crop.Ksccx
        CGC=crop.CGC#*crop.Ksexpf
        
        Half_CCx = round(crop.Emergence+(np.log(0.5*CCx/crop.CC0)/CGC))
            
            
        Full_CCx = round(crop.Emergence+(np.log((0.25*CCx*CCx/crop.CC0)
                                                                    /(CCx-(0.98*CCx)))/CGC))
        
        if gdd_cum.values[-1] < crop.Maturity:
            Half_CCx = Half_CCx*gdd_cum.values[-1]/crop.Maturity
            Full_CCx = Full_CCx*gdd_cum.values[-1]/crop.Maturity

        Half_CCxCD = (gdd_cum>Half_CCx).idxmax()+1 # isn't used, so why is it defined?
        Full_CCxCD = (gdd_cum>Full_CCx).idxmax()+1
        #Full_CCxCD=crop.MaxCanopyCD

        Ks_Tr=[]# cold stress 
        Kc_Tr=[]# crop transpiration coefficient with soil fertility stress 
        Ksc_Total=[]
        max_cc=0

        for day_ in range(1,np.min([crop.MaturityCD+1,len(gdd_cum)])):
            #cold stress
            GDD_=gdd_cum.values[day_]-gdd_cum.values[day_-1]
            
            #copy from solution.py
            if crop.TrColdStress == 0:
            # Cold temperature stress does not affect transpiration
                KsCold = 1
            elif crop.TrColdStress == 1:
                # Transpiration can be affected by cold temperature stress
                if GDD_ >= crop.GDD_up:
                    # No cold temperature stress
                    KsCold = 1
                elif GDD_ <= crop.GDD_lo:
                    # Transpiration fully inhibited by cold temperature stress
                    KsCold = 0
                else:
                    # Transpiration partially inhibited by cold temperature stress
                    # Get parameters for logistic curve
                    KsTr_up = 1
                    KsTr_lo = 0.02
                    fshapeb = (-1) * (
                        np.log(((KsTr_lo * KsTr_up) - 0.98 * KsTr_lo) / (0.98 * (KsTr_up - KsTr_lo)))
                    )
                    # Calculate cold stress level
                    GDDrel = (GDD_ - crop.GDD_lo) / (crop.GDD_up - crop.GDD_lo)
                    KsCold = (KsTr_up * KsTr_lo) / (
                        KsTr_lo + (KsTr_up - KsTr_lo) * np.exp(-fshapeb * GDDrel)
                    )
                    KsCold = KsCold - KsTr_lo * (1 - GDDrel)
            #record cold stress 
            Ks_Tr.append(KsCold)

            if gdd_cum.values[day_]<crop.Emergence:
                CC=0
                Kctr=crop.Kcb
            
            elif gdd_cum.values[day_] <= Half_CCx:
                CC=crop.CC0*np.exp((gdd_cum.values[day_]-crop.Emergence)*CGC)
                if CC>CCx/2:
                    CC=CCx-0.25*CCx*CCx/crop.CC0*np.exp(-(gdd_cum.values[day_]-crop.Emergence)*CGC)
                Kctr=crop.Kcb
                
                max_cc=CC

            elif gdd_cum.values[day_] > Half_CCx and gdd_cum.values[day_] <= Full_CCx:
                CC=CCx-0.25*CCx*CCx/crop.CC0*np.exp(-(gdd_cum.values[day_]-crop.Emergence)*CGC)
                Kctr=crop.Kcb
                
                max_cc=CC
                    
            elif gdd_cum.values[day_] > Full_CCx and gdd_cum.values[day_-5] <= Full_CCx:
            
                if gdd_cum.values[day_]<crop.CanopyDevEnd:
                    CC=CCx-0.25*CCx*CCx/crop.CC0*np.exp(-(gdd_cum.values[day_]-crop.Emergence)*CGC)
                    max_cc=CC
                else:
                    CC=max_cc
                    
                Kctr=crop.Kcb

            elif gdd_cum.values[day_-5] > Full_CCx and gdd_cum.values[day_] <= crop.Senescence:
                
                if gdd_cum.values[day_]<crop.CanopyDevEnd:
                    CC=CCx-0.25*CCx*CCx/crop.CC0*np.exp(-(gdd_cum.values[day_]-crop.Emergence)*CGC)
                    max_cc=CC
                else:
                    CC=max_cc

                Kctr=crop.Kcb-(day_-Full_CCxCD-5)*(crop.fage / 100)*max_cc

            elif gdd_cum.values[day_] > crop.Senescence and gdd_cum.values[day_] <= crop.Maturity:
                CC_adj=max_cc
                
                CDC = crop.CDC*((CC_adj+2.29)/(crop.CCx+2.29))
                CC=CC_adj*(1-0.05*(np.exp(3.33*CDC*(gdd_cum.values[day_]-crop.Senescence)/(CC_adj+2.29))-1))
                Kctr=(crop.Kcb-(day_-Full_CCxCD-5)*(crop.fage / 100)*max_cc)*(CC/max_cc)**crop.a_Tr
                
            if CC<=0:
                CC=0
            CC_star=1.72*CC-CC*CC+0.3*CC*CC*CC

            #print(CC)
            #print(ParamStruct.CO2.CurrentConc)
            try:
                CO2conc=ParamStruct.CO2.current_concentration 
            except:
                CO2conc=ParamStruct.CO2.co2_data_processed.iloc[0]
            Kc_TrCo2=1            
            if CO2conc>369.41:
                Kc_TrCo2=1-0.05*(CO2conc-369.41)/(550-369.41)
            
            Kc_Tr_=CC_star*Kctr*Kc_TrCo2
            #print(day_)
            #print(CC_star)
            #print(Kctr)
            #print(Kc_TrCo2)
            #print(Kc_Tr_)
            Kc_Tr.append(Kc_Tr_)

            Ksc_Total.append(Kc_Tr_*KsCold)
        
        crop.TR_ET0_fertstress=np.sum(Ksc_Total[:])
        #print(Ksc_Total[:])
        #print(crop.MaturityCD)
        
        crop.HIGC = calculate_HIGC(
            crop.YldFormCD,
            crop.HI0,
            crop.HIini,
        )
        if crop.CropType == 3:
            # Determine linear switch point and HIGC rate for fruit/grain crops
            crop.tLinSwitch, crop.dHILinear = calculate_HI_linear(
                crop.YldFormCD, crop.HIini, crop.HI0, crop.HIGC
            )
        else:
            # No linear switch for leafy vegetable or root/tiber crops
            crop.tLinSwitch = 0
            crop.dHILinear = 0.
        
        Bio_top=0
        for i_ in range(len(Ksc_Total)):
            #print(Bio_top)
            crop.Bio_top[i_]=Bio_top
            if i_>=crop.HIstartCD:

                if ((crop.CropType == 2) or (crop.CropType == 3)):
                
                    if crop.CropType == 2:
                        PctLagPhase_=100
                    else:
                        if (i_-crop.HIstartCD) < crop.tLinSwitch:
                            PctLagPhase_ = 100*((i_-crop.HIstartCD)/crop.tLinSwitch)
                        else:
                            PctLagPhase_=100

                    # Adjust WP for reproductive stage
                    if crop.Determinant == 1:
                        fswitch = PctLagPhase_/100
                    else:
                        if (i_-crop.HIstartCD) < (crop.YldFormCD/3):
                            fswitch = (i_-crop.HIstartCD)/(crop.YldFormCD/3)
                        else:
                            fswitch = 1

                    if fswitch>1:
                        fswitch=1
                        
                    Bio_top+= Ksc_Total[i_]*(1-(1-crop.WPy/100)*fswitch)
                else:
                    Bio_top+= Ksc_Total[i_]
            else:
                Bio_top+= Ksc_Total[i_]
        crop.Bio_top[len(Ksc_Total):len(Ksc_Total)+100]=Bio_top

    if crop.soil_fert_stress == 1:
        return crop, gdd_cum, Ksc_Total
    else:
        return crop
