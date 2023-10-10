import numpy as np
import pandas as pd

from ..entities.paramStruct import ParamStruct
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    # Important: classes are only imported when types are checked, not in production.
    from aquacrop.entities.crop import Crop
    from pandas import DatetimeIndex, DataFrame

from .calculate_HIGC import calculate_HIGC
from .calculate_HI_linear import calculate_HI_linear

def calibrate_soil_fert_stress(
    crop: "Crop",
    gdd_cum: "pd.Series",
    ParamStruct: "ParamStruct",
) -> "Crop":
    """
    Function to compute additional parameters for soil fertility stress


    Arguments:

        crop (Crop):  Crop object containing crop paramaters

        gdd_cum (pd.Series):  cumulative GDD values throughout season

        param_struct (ParamStruct):  Contains model crop and soil paramaters


    Returns:

        crop (Crop): updated Crop object



    """

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

#soil fertility stress initialization and calibration, basically calculate the whole crop growth under ideal conditions (no water stress)  

        #This part can be optimized(shorten) later for the final version
        
        #print(crop.EmergenceCD)
        #print(crop.MaxCanopy)
        #print(crop.MaxCanopyCD)
        #print(crop.CanopyDevEndCD)

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
    
    
    def Biomas_ini_es(Bio_mul,TopStress,Ksccx_temp,Ksexpf_temp,fcdecline_temp,Kswp_temp):
        CCx=crop.CCx*Ksccx_temp
    
        CGC=crop.CGC*Ksexpf_temp
        
        Half_CCx = round(crop.Emergence+(np.log(0.5*CCx/crop.CC0)/CGC))
        
        
        Full_CCx = round(crop.Emergence+(np.log((0.25*CCx*CCx/crop.CC0)
                                                                    /(CCx-(0.98*CCx)))/CGC))
        
        if gdd_cum.values[-1] < crop.Maturity:
            Half_CCx = Half_CCx*gdd_cum.values[-1]/crop.Maturity
            Full_CCx = Full_CCx*gdd_cum.values[-1]/crop.Maturity

        Half_CCxCD = (gdd_cum>Half_CCx).idxmax()+1 # isn't used, so why is it defined?
        Full_CCxCD = (gdd_cum>Full_CCx).idxmax()+1

        Kc_Tr_es=[]# crop transpiration coefficient with soil fertility stress 
        max_cc=0

        for day_ in range(1,np.min([crop.MaturityCD+1,len(gdd_cum)])):

            # crop transpiration coefficient
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

                if crop.SenescenceCD>Full_CCxCD:
                    CC=max_cc-fcdecline_temp*(day_-Full_CCxCD)*(day_-Full_CCxCD)/(crop.SenescenceCD-Full_CCxCD)
                Kctr=crop.Kcb
                if CC<0:
                    CC=0

            elif gdd_cum.values[day_-5] > Full_CCx and gdd_cum.values[day_] <= crop.Senescence:
                
                if gdd_cum.values[day_]<crop.CanopyDevEnd:
                    CC=CCx-0.25*CCx*CCx/crop.CC0*np.exp(-(gdd_cum.values[day_]-crop.Emergence)*CGC)
                    max_cc=CC
                else:
                    CC=max_cc


                if crop.SenescenceCD>Full_CCxCD:
                    CC=max_cc-fcdecline_temp*(day_-Full_CCxCD)*(day_-Full_CCxCD)/(crop.SenescenceCD-Full_CCxCD)
                Kctr=crop.Kcb-(day_-Full_CCxCD-5)*(crop.fage / 100)*max_cc
                if CC<0:
                    CC=0

            elif gdd_cum.values[day_] > crop.Senescence and gdd_cum.values[day_] <= crop.Maturity:
                if crop.SenescenceCD>Full_CCxCD:
                    CC_fs=max_cc-fcdecline_temp*(day_-Full_CCxCD)*(day_-Full_CCxCD)/(crop.SenescenceCD-Full_CCxCD)

                    CC_adj=max_cc-fcdecline_temp*(crop.SenescenceCD-Full_CCxCD)
                else:
                    CC_fs=max_cc
                    CC_adj=max_cc
                CDC = crop.CDC*((CC_adj+2.29)/(crop.CCx+2.29))
                CC=CC_adj*(1-0.05*(np.exp(3.33*CDC*(gdd_cum.values[day_]-crop.Senescence)/(CC_adj+2.29))-1))
                
                #if CC_fs<CC:#not in AquaCrop v6
                #    CC=CC_fs
                if CC<0:
                    CC=0
                
                Kctr=(crop.Kcb-(day_-Full_CCxCD-5)*(crop.fage / 100)*max_cc)*(CC/max_cc)**crop.a_Tr
            
            #print(Full_CCxCD)
            #print(Kc_Tr_es)
            #print(Kc_Tr)
            
            if CC<=0:
                CC=0
            CC_star=1.72*CC-CC*CC+0.3*CC*CC*CC
            
            #print(ParamStruct.CO2.CurrentConc)
            try:
                CO2conc=ParamStruct.CO2.current_concentration 
            except:
                CO2conc=ParamStruct.CO2.co2_data_processed.iloc[0]
            Kc_TrCo2=1            
            if CO2conc>369.41:
                Kc_TrCo2=1-0.05*(CO2conc-369.41)/(550-369.41)
            
            Kc_Tr_=CC_star*Kctr*Kc_TrCo2
            #if Kc_Tr_==0:
            #    print("Zero"+str(day_))
            #    print(CC_star)
            #    print(Kctr)
            Kc_Tr_es.append(Kc_Tr_)

        Bio_cur=0
        Curstress=0
        #print(len(Kc_Tr_es))
        #print(len(Kc_Tr))
        #print(len(Bio_mul))
        for i in range(len(Kc_Tr_es)):
            if Curstress<TopStress:
                Curstress+=Kc_Tr_es[i]*Ks_Tr[i]
                Kswp_=1-(1-Kswp_temp)*(Curstress/TopStress)*(Curstress/TopStress)
            else:
                Kswp_=1-(1-Kswp_temp)
            Bio_cur+=Kswp_*Kc_Tr_es[i]*Ks_Tr[i]*Bio_mul[i]
        
        return Bio_cur
    
    
    if crop.soil_fert_stress==1:
        
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
            # No linear switch for leafy vegetable or root/tuber crops
            crop.tLinSwitch = 0
            crop.dHILinear = 0.

        TopStress=crop.TR_ET0_fertstress*crop.RelativeBio
        
        Bio_top=0
        Bio_mul=[]
        #print('TOtal')
        #print(len(Ksc_Total))
        for i_ in range(len(Ksc_Total)):
            #print(i_)
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
                    Bio_mul.append(1-(1-crop.WPy/100)*fswitch)
                    
                else:
                    Bio_top+= Ksc_Total[i_]
                    Bio_mul.append(1)
            else:
                Bio_top+= Ksc_Total[i_]
                Bio_mul.append(1)
        
        #initialze parameters
        adj_flag=True
        Ksccx_temp=crop.Ksccx_in
        L123=crop.SenescenceCD
        L12=crop.MaxCanopyCD
        kk_=0#avoid periodic loop
        while adj_flag and kk_<100:
            kk_+=1

            CCx_temp=crop.CCx*Ksccx_temp # not used
            if crop.fcdecline_in==0:
                CCxfinal=0.92*crop.CCx
            elif crop.fcdecline_in==1:
                CCxfinal=0.85*crop.CCx
            else:
                CCxfinal=0.77*crop.CCx
                
            if L123>L12:
                fcdecline_temp=(1-Ksccx_temp)*(crop.CCx-CCxfinal)/(L123-L12)
                if fcdecline_temp>0.01:
                    fcdecline_temp=0.009
            else:
                fcdecline_temp=0
            
            Ksexpf_temp=np.round((1-(1-Ksccx_temp)*0.6)*100)/100
            Kswp_temp=np.round((1-(1-Ksccx_temp)*0.5)*100)/100 
            
            full_temp = round(crop.EmergenceCD+(np.log((0.25*crop.CCx*Ksccx_temp*crop.CCx*Ksccx_temp/crop.CC0)
                                                                /(crop.CCx*Ksccx_temp-(0.98*crop.CCx*Ksccx_temp)))/crop.CGC_CD/Ksexpf_temp))
            

            if full_temp>crop.CanopyDevEndCD:

                while Ksexpf_temp<0.99 and full_temp>crop.CanopyDevEndCD:
                    Ksexpf_temp=Ksexpf_temp+0.01
                    full_temp = round(crop.EmergenceCD+(np.log((0.25*crop.CCx*Ksccx_temp*crop.CCx*Ksccx_temp/crop.CC0)
                                                                /(crop.CCx*Ksccx_temp-(0.98*crop.CCx*Ksccx_temp)))/crop.CGC_CD/Ksexpf_temp))
                                                        
                while full_temp>crop.CanopyDevEndCD and crop.CCx*Ksccx_temp>0.1 and Ksccx_temp>0.5:
                    Ksccx_temp=Ksccx_temp-0.01
                    full_temp = round(crop.EmergenceCD+(np.log((0.25*crop.CCx*Ksccx_temp*crop.CCx*Ksccx_temp/crop.CC0)
                                                                /(crop.CCx*Ksccx_temp-(0.98*crop.CCx*Ksccx_temp)))/crop.CGC_CD/Ksexpf_temp))

            adj_wp=True
            k_=0#avoid periodic loop
            while adj_wp and k_<100:
                #print('wp')
                k_+=1
                bio_cur=Biomas_ini_es(Bio_mul,TopStress,Ksccx_temp,Ksexpf_temp,fcdecline_temp,Kswp_temp)
                #print(Kswp_temp)
                #print(int(bio_cur/Bio_top*100))
                if int(bio_cur/Bio_top*100)==int(crop.RelativeBio*100):
                    adj_flag=False
                    adj_wp=False
                elif bio_cur/Bio_top<crop.RelativeBio:
                    Kswp_temp=Kswp_temp+0.01
                else:
                    Kswp_temp=Kswp_temp-0.01
                
                if Kswp_temp<0.3-0.001:
                    Kswp_temp=Kswp_temp+0.01
                    adj_wp=False
                if Kswp_temp>0.99+0.001:
                    Kswp_temp=Kswp_temp-0.01
                    adj_wp=False
                    
                #print(Ksccx_temp)
                #print(Ksexpf_temp)
                #print(Kswp_temp)
                #print(int(bio_cur/Bio_top*100))
                #print(bio_cur)
                #print(Bio_top)
                #print(int(crop.RelativeBio*100))
            
            if adj_flag:
                adj_Ksexpf=True
                k_=0
                while adj_Ksexpf and k_<100:
                    k_+=1
                    bio_cur=Biomas_ini_es(Bio_mul,TopStress,Ksccx_temp,Ksexpf_temp,fcdecline_temp,Kswp_temp)
                    if int(bio_cur/Bio_top*100)==int(crop.RelativeBio*100):
                        adj_flag=False
                        adj_Ksexpf=False
                    elif bio_cur/Bio_top<crop.RelativeBio:
                        Ksexpf_temp=Ksexpf_temp+0.01
                    else:
                        Ksexpf_temp=Ksexpf_temp-0.01
                    
                    full_temp = round(crop.EmergenceCD+(np.log((0.25*crop.CCx*Ksccx_temp*crop.CCx*Ksccx_temp/crop.CC0)
                                                                /(crop.CCx*Ksccx_temp-(0.98*crop.CCx*Ksccx_temp)))/crop.CGC_CD/Ksexpf_temp))
            
                    if full_temp>crop.CanopyDevEndCD:
                        while Ksexpf_temp<0.99 and full_temp>crop.CanopyDevEndCD:
                            Ksexpf_temp=Ksexpf_temp+0.01
                            full_temp = round(crop.EmergenceCD+(np.log((0.25*crop.CCx*Ksccx_temp*crop.CCx*Ksccx_temp/crop.CC0)
                                                                /(crop.CCx*Ksccx_temp-(0.98*crop.CCx*Ksccx_temp)))/crop.CGC_CD/Ksexpf_temp))
                        while full_temp>crop.CanopyDevEndCD and crop.CCx*Ksccx_temp>0.1 and Ksccx_temp>0.5:
                            Ksccx_temp=Ksccx_temp-0.01
                            full_temp = round(crop.EmergenceCD+(np.log((0.25*crop.CCx*Ksccx_temp*crop.CCx*Ksccx_temp/crop.CC0)
                                                                /(crop.CCx*Ksccx_temp-(0.98*crop.CCx*Ksccx_temp)))/crop.CGC_CD/Ksexpf_temp))

                        adj_Ksexpf=False
                    
                    if (1-Ksexpf_temp)<0.1*(1-Ksccx_temp)-0.001:
                        Ksexpf_temp=Ksexpf_temp-0.01
                        adj_Ksexpf=False
                    if (1-Ksexpf_temp)>0.8*(1-Ksccx_temp)+0.001:
                        Ksexpf_temp=Ksexpf_temp+0.01
                        adj_Ksexpf=False
                        
            if adj_flag:
                adj_fcdecline=True
                k_=0
                while adj_fcdecline and k_<100:
                    k_+=1
                    #print('fcdecline')
                    bio_cur=Biomas_ini_es(Bio_mul,TopStress,Ksccx_temp,Ksexpf_temp,fcdecline_temp,Kswp_temp)

                    if int(bio_cur/Bio_top*100)==int(crop.RelativeBio*100):
                        adj_flag=False
                        adj_fcdecline=False
                    elif bio_cur/Bio_top<crop.RelativeBio:
                        fcdecline_temp=fcdecline_temp-0.0001
                    else:
                        fcdecline_temp=fcdecline_temp+0.0001
                    
                    if fcdecline_temp<0.0001-0.00001:
                        fcdecline_temp=fcdecline_temp+0.0001
                        adj_fcdecline=False
                    if fcdecline_temp>0.009+0.00001:
                        fcdecline_temp=fcdecline_temp-0.0001
                        adj_fcdecline=False
                    
            if adj_flag:

                if bio_cur/Bio_top<crop.RelativeBio:
                    Ksccx_temp=Ksccx_temp+0.01
                else:
                    Ksccx_temp=Ksccx_temp-0.01
                    
                if Ksccx_temp<0:
                    Ksccx_temp=Ksccx_temp+0.01
                    adj_flag=False
                if Ksccx_temp>1:
                    Ksccx_temp=Ksccx_temp-0.01
                    adj_flag=False
            
            #print(int(bio_cur/Bio_top*100))
            #print(int(crop.RelativeBio*100))
            #print(adj_flag)
            CCx=crop.CCx*Ksccx_temp
    
            CGC=crop.CGC*Ksexpf_temp

            Full_CCx = round(crop.Emergence+(np.log((0.25*CCx*CCx/crop.CC0)
                                                                        /(CCx-(0.98*CCx)))/CGC))
            if gdd_cum.values[-1] < crop.Maturity:

                Full_CCx = Full_CCx*gdd_cum.values[-1]/crop.Maturity
    
            Full_CCxCD = (gdd_cum>Full_CCx).idxmax()+1
            L12=Full_CCxCD
            
            #if L12>L123:
            #    print("L12>L123")


        if (int(bio_cur/Bio_top*100)==int(crop.RelativeBio*100)) or (kk_>=100 and np.abs(int(bio_cur/Bio_top*100)-int(crop.RelativeBio*100))<5):
        #accept some errors when exit because of maximum loop, in theory, the above can ensure ==

            crop.Ksccx_in=Ksccx_temp
            crop.Ksccx_es[0]=Ksccx_temp
            crop.Ksexpf_es[0]=Ksexpf_temp
            crop.fcdecline_es[0]=fcdecline_temp
            crop.Kswp_es[0]=Kswp_temp
            crop.sf_es[0]=1-bio_cur/Bio_top
            crop.relbio_es[0]=bio_cur/Bio_top
            
            def get_shape(ks,sf):
                if ks>1-sf:
                    up_=50
                    low_=0
                else:
                    up_=0
                    low_=-50
                flag_=True
                while flag_:
                    shape_=(up_+low_)/2
                    ks_=1-(np.exp(sf*shape_)-1)/(np.exp(shape_)-1)
                    
                    if np.abs(ks_-ks)<0.0001 or np.abs(up_-low_)<0.0001:
                        flag_=False
                    elif ks_<ks:
                        low_=shape_
                    else:
                        up_=shape_
                return shape_
            
            Ksccx_shape=get_shape(Ksccx_temp,crop.sf_es[0])
            Ksexpf_shape=get_shape(Ksexpf_temp,crop.sf_es[0])
            Kswp_shape=get_shape(Kswp_temp,crop.sf_es[0])
            fcdecline_shape=get_shape(1-fcdecline_temp*100,crop.sf_es[0])
            
            for i in range(0,100):
                crop.sf_es[i+1]=i/100
                Ksccx_temp=1-(np.exp(crop.sf_es[i+1]*Ksccx_shape)-1)/(np.exp(Ksccx_shape)-1)
                Ksexpf_temp=1-(np.exp(crop.sf_es[i+1]*Ksexpf_shape)-1)/(np.exp(Ksexpf_shape)-1)
                Kswp_temp=1-(np.exp(crop.sf_es[i+1]*Kswp_shape)-1)/(np.exp(Kswp_shape)-1)
                fcdecline_temp=(np.exp(crop.sf_es[i+1]*fcdecline_shape)-1)/(np.exp(fcdecline_shape)-1)/100
                bio_temp=Biomas_ini_es(Bio_mul,crop.TR_ET0_fertstress*(1-crop.sf_es[i+1]),Ksccx_temp,Ksexpf_temp,fcdecline_temp,Kswp_temp)

                crop.Ksccx_es[i+1]=Ksccx_temp
                crop.Ksexpf_es[i+1]=Ksexpf_temp
                crop.fcdecline_es[i+1]=fcdecline_temp
                crop.Kswp_es[i+1]=Kswp_temp
                crop.relbio_es[i+1]=bio_temp/Bio_top


    if crop.soil_fert_stress==2:
        
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

        TopStress=crop.TR_ET0_fertstress*crop.RelativeBio
        
        possible_flag=1
        crop.Ksccx_in=crop.Ksccx_in-0.01
        ind_=0
        
        while possible_flag:
            crop.Ksccx_in=crop.Ksccx_in+0.01
            CCx=crop.CCx*crop.Ksccx_in
            
            Bio_top=0
            Bio_mul=[]
            for i_ in range(len(Ksc_Total)):
                #print(Bio_top)
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
                        Bio_mul.append(1-(1-crop.WPy/100)*fswitch)
                    else:
                        Bio_top+= Ksc_Total[i_]
                        Bio_mul.append(1)
                else:
                    Bio_top+= Ksc_Total[i_]
                    Bio_mul.append(1)


            for Ksexpf in np.arange(1,101,2)/100:
                for fcdecline in np.arange(0,100,2)/100/100:
                    
                    #Ksexpf=0.79
                    #fcdecline=0.0006
                
                    
                    CGC=crop.CGC*Ksexpf
                    Half_CCx = round(crop.Emergence+(np.log(0.5*CCx/crop.CC0)/CGC))
        
        
                    Full_CCx = round(crop.Emergence+(np.log((0.25*CCx*CCx/crop.CC0)
                                                                                /(CCx-(0.98*CCx)))/CGC))
                    
                    if gdd_cum.values[-1] < crop.Maturity:
                        Half_CCx = Half_CCx*gdd_cum.values[-1]/crop.Maturity
                        Full_CCx = Full_CCx*gdd_cum.values[-1]/crop.Maturity
            
                    Half_CCxCD = (gdd_cum>Half_CCx).idxmax()+1 # isn't used, so why is it defined?
                    Full_CCxCD = (gdd_cum>Full_CCx).idxmax()+1
                    
                    if Full_CCxCD>=crop.SenescenceCD:
                        break
                    
                    #Full_CCxCD=crop.MaxCanopyCD

                    Kc_Tr_es=[]# crop transpiration coefficient with soil fertility stress 
                    max_cc=0

                    for day_ in range(1,np.min([crop.MaturityCD+1,len(gdd_cum)])):

                        # crop transpiration coefficient
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
                        
                            if crop.SenescenceCD>Full_CCxCD:
                                CC=max_cc-fcdecline_temp*(day_-Full_CCxCD)*(day_-Full_CCxCD)/(crop.SenescenceCD-Full_CCxCD)
                            Kctr=crop.Kcb
                            if CC<0:
                                CC=0
        
                        elif gdd_cum.values[day_-5] > Full_CCx and gdd_cum.values[day_] <= crop.Senescence:
                        
                            if gdd_cum.values[day_]<crop.CanopyDevEnd:
                                CC=CCx-0.25*CCx*CCx/crop.CC0*np.exp(-(gdd_cum.values[day_]-crop.Emergence)*CGC)
                                max_cc=CC
                            else:
                                CC=max_cc
                        
                            if crop.SenescenceCD>Full_CCxCD:
                                CC=max_cc-fcdecline_temp*(day_-Full_CCxCD)*(day_-Full_CCxCD)/(crop.SenescenceCD-Full_CCxCD)
                            Kctr=crop.Kcb-(day_-Full_CCxCD-5)*(crop.fage / 100)*max_cc
                            if CC<0:
                                CC=0
        
                        elif gdd_cum.values[day_] > crop.Senescence and gdd_cum.values[day_] <= crop.Maturity:
                            if crop.SenescenceCD>Full_CCxCD:
                                CC_fs=max_cc-fcdecline_temp*(day_-Full_CCxCD)*(day_-Full_CCxCD)/(crop.SenescenceCD-Full_CCxCD)
        
                                CC_adj=max_cc-fcdecline_temp*(crop.SenescenceCD-Full_CCxCD)
                            else:
                                CC_fs=max_cc
                                CC_adj=max_cc
                            CDC = crop.CDC*((CC_adj+2.29)/(crop.CCx+2.29))
                            CC=CC_adj*(1-0.05*(np.exp(3.33*CDC*(gdd_cum.values[day_]-crop.Senescence)/(CC_adj+2.29))-1))
                            
                            #if CC_fs<CC:
                            #    CC=CC_fs
                            if CC<0:
                                CC=0
                            
                            Kctr=(crop.Kcb-(day_-Full_CCxCD-5)*(crop.fage / 100)*max_cc)*(CC/max_cc)**crop.a_Tr
        
                        if CC<=0:
                            CC=0
                        CC_star=1.72*CC-CC*CC+0.3*CC*CC*CC
                        
                        #print(ParamStruct.CO2.CurrentConc)
                        try:
                            CO2conc=ParamStruct.CO2.current_concentration 
                        except:
                            CO2conc=ParamStruct.CO2.co2_data_processed.iloc[0]
                        Kc_TrCo2=1            
                        if CO2conc>369.41:
                            Kc_TrCo2=1-0.05*(CO2conc-369.41)/(550-369.41)
                        
                        Kc_Tr_=CC_star*Kctr*Kc_TrCo2
                        #if Kc_Tr_==0:
                        #    print("Zero"+str(day_))
                        #    print(CC_star)
                        #    print(Kctr)
                        Kc_Tr_es.append(Kc_Tr_)
                    
                    Bio_up=0
                    for i in range(len(Kc_Tr_es)):
                        Bio_up+=Kc_Tr_es[i]*Ks_Tr[i]*Bio_mul[i]
                        #print(Bio_up)
                        
                    if int(Bio_up/Bio_top*100)<int(crop.RelativeBio*100):
                        #pass
                        break
                        
                    Kswp_l=0
                    Kswp_u=1
                    con_=True
                    
                    while con_:
                        
                        Kswp=(Kswp_l+Kswp_u)/2
                        #Kswp=0.6

                        Bio_cur=0
                        Curstress=0
                        
                        for i in range(len(Kc_Tr_es)):
                            if Curstress<TopStress:
                                Curstress+=Kc_Tr_es[i]*Ks_Tr[i]
                                Kswp_=1-(1-Kswp)*(Curstress/TopStress)*(Curstress/TopStress)
                            else:
                                Kswp_=1-(1-Kswp)
                            Bio_cur+=Kswp_*Kc_Tr_es[i]*Ks_Tr[i]*Bio_mul[i]
                            #print(Bio_cur)
                            

                        if int(Bio_cur/Bio_top*100)==crop.RelativeBio*100:
                            con_=False
                            possible_flag=0
                            crop.Ksexpf_es[ind_]=Ksexpf
                            crop.fcdecline_es[ind_]=fcdecline
                            crop.Kswp_es[ind_]=Kswp

                            ind_+=1
                            
                            break
                        else:
                            if Bio_cur/Bio_top>crop.RelativeBio:
                                Kswp_u=Kswp
                            else:
                                Kswp_l=Kswp

    

    return crop