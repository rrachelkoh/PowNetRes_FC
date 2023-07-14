# import datetime
# yr_start = 2000
# yr_end = 2003
res_ops = {}
# =============================================================================
# #c: runoff coeff, A: basin area in m2, V: reservoir storage capacity, 
#Hmax/min: max/min res lvl, hydhead: hydr head (dist from max WL to turbine), 
#DayMin/Max: day with lowest/highest water level
# Hmin_Kamchay=592; Hmax_Kamchay=602; DayMin_Kamchay=180; DayMax_Kamchay=298

res_sysparam = { 'KAMCHAY': { 'c':0.62, 'A':710000000, 'V':432000000, 
                          'Hmin':500, 'Hmax':610, 'hydhead':122,
                          'DayMin':180, 'DayMax':298, 'Release':163.5,
                          'plant_cap':194.1}, #, 's_ini':311599190
                'ATAY': { 'c':0.75, 'A':1157000000, 'V':443800000, 
                          'Hmin':510, 'Hmax':545, 'hydhead':216,
                          'DayMin':158, 'DayMax':311, 'Release':125,
                          'plant_cap':240}, #, 's_ini':310660000
                'LRCHRUM': { 'c':0.51, 'A':1550000000, 'V':62000000, 
                          'Hmin':500, 'Hmax':568, 'hydhead':132,
                          'DayMin':150, 'DayMax':310, 'Release':300,
                          'plant_cap':338}, #, 's_ini':43400000
                'TATAY': { 'c':0.71, 'A':1073000000, 'V':322000000, 
                          'Hmin':500, 'Hmax':577, 'hydhead':188,
                          'DayMin':150, 'DayMax':310, 'Release':150,
                          'plant_cap':246}, #, 's_ini':225400000
                'KIRIROM1': { 'c':0.48, 'A':99000000, 'V':30000000, 
                          'Hmin':500, 'Hmax':534, 'hydhead':373.5,
                          'DayMin':150, 'DayMax':310, 'Release':20,
                          'plant_cap':12}, #, 's_ini':21000000
                'KIRIROM3': { 'c':0.48, 'A':105000000, 'V':30000000, 
                          'Hmin':500, 'Hmax':540, 'hydhead':271,
                          'DayMin':150, 'DayMax':310, 'Release':40,
                          'plant_cap':18}}  #, 's_ini':21000000



cascade = { 'Upper':{'1': ['ATAY']}, 'Lower':{'1': 'LRCHRUM'} } 
# cascade = { 'Upper':{'1': ['ATAY','KAMCHAY'],'2':['KIRIROM1']}, 'Lower':{'1': 'LRCHRUM','2':'KIRIROM3'} } 
#all length units in [m]

# Change model parameters when appropriate. 
# Note: if you change Hmax/Hmin, you have to change V
# If change Hmax, Hmin, Tmax, Tmin to an appropriate value, you can use dams
# to simulate the impact of soil storage so no need to redistribute rainfall

#rainfall_all = pd.read_csv('precip_1976to2007.csv')
# Q = pd.read_csv('Q_in_2016.csv')

# if yr%4 == 0:
#     Q = Q.drop([Q.index[59]])
#-------------------------------------------------

# Constants
# name = 'KAMCHAY'
# yr_start = 2000
# yr_end = 2000
# yr = yr_start
# turbine_factor = 0.9
# APH_cor_factor = 1.26 # APHRODITE correction factor 
# because APH underestimates real rainfall --> x correction factor

# res_ops = dict.fromkeys(res_sysparam.keys(), {})
# res_ops = {}

def Model1a(p1):
    # if num_yrs == 5:
    #     Q  = pd.read_csv('../Data/data_camb_Q_12to16.csv')
    # elif num_yrs == 37:
    #     Q  = pd.read_csv('../Data/data_camb_Q_80to16.csv')
    
    # else:
    #     Q  = pd.read_csv('../Data/Q_1000/data_camb_Q_syn_run'+str(run_no)+'.csv')
    # Q = pd.concat([Q] * num_yrs, ignore_index=True)
    

    Q_MAF = pd.read_csv('../Data/data_camb_Q_MAF.csv')
    Q_MAF = pd.concat([Q_MAF] * (num_yrs+1), ignore_index=True)
    # if yr%4 == 0:
    #     Q = Q.drop([Q.index[59]])
    
    for name in res_sysparam.keys(): 
        # Q = pd.read_csv('glofas_'+name+'_20to22_obs.csv')
        Q = pd.read_csv('../Data/Qglofas_2000to2019.csv')
        
        Q['Date'] = pd.to_datetime(Q['Date'],dayfirst=True)
        Q = Q[~((Q.Date.dt.month == 2) & (Q.Date.dt.day == 29))]
                
        date_start = pd.to_datetime(str(yr_start)+'-1-1')
        date_end = pd.to_datetime(str(yr_end)+'-12-31') + timedelta(days=p1)

        Q = Q[(Q['Date'] >= date_start) & (Q['Date'] <= date_end)]
        Q  = Q.drop(['Date'],axis=1)
    
        ### *** DIRECT INFLOW INPUTS *** ###
        # if FC_typ == 'perfect' or FC_typ == 'clim':
            # Q_in = Q.iloc[:,0]
        Q_in = np.array(Q[name])
        Q_in = np.append(np.nan,Q_in) #arbitrary value for day 0
        Q_raw = list(Q_in)
        
        Q_MAF_ = np.array(Q_MAF[name])
        Q_MAF_ = np.append(np.nan,Q_MAF_) #arbitrary value for day 0
        Q_MAF_ = list(Q_MAF_)
        
        if FC_typ == 'FC':
            Q_FC = pd.read_csv('../Data/Q_FC/data_camb_Q_FCmember'+str(mem_num)+'_'+name+'.csv')
            Q_FC['Date'] = pd.to_datetime(Q_FC['Date'])
            Q_FC = Q_FC[~((Q_FC.Date.dt.month == 2) & (Q_FC.Date.dt.day == 29))]
            Q_FC = Q_FC[(Q_FC['Date'] >= date_start) & (Q_FC['Date'] <= date_end)]
            Q_FC  = Q_FC.drop(['Date'],axis=1)
        
        # Design model parameters
        DayMin = res_sysparam[name]['DayMin']
        DayMax = res_sysparam[name]['DayMax']
        H = res_sysparam[name]['hydhead']
        Hmin = res_sysparam[name]['Hmin']
        Hmax = res_sysparam[name]['Hmax']
        V = res_sysparam[name]['V']
        # s_ini = V*.7
        Qmax = res_sysparam[name]['Release']
        length = num_yrs*365+1#len(Q_MAF_)
        
        #-------------------------------------------------
        #Reservoir operation
        s = np.nan * np.ones(shape = (length,1)) #storage from opt operations
        Hd_m = np.nan * np.ones(shape = (length,1))
        r = np.nan * np.ones(shape = (length,1)) #release from opt operations
        Qspil = np.nan * np.ones(shape = (length,1)) #
        Q_MEF = np.nan * np.ones(shape = (length,1)) #
        availhydro = np.nan * np.ones(shape = (length,1)) #available hydropower
        day_reop = np.zeros(shape = (length,1)) #which day re-operated
        
        #-------------------------------------------------
        # Initialize for different operating methods #
        #-------------------------------------------------
        if opt_mtd=='none' or opt_mtd=='mpc_rc':
        #Operating rule of reservoir
            Hdesign=[] #Design water level
            Vdesign=[] #Design storage
            # DayNum=[]
                    
            for t in range(1,366):
                DayNo = t%366
                if DayMin<DayNo<DayMax :
                    Hdes = (DayNo-DayMin)/(DayMax-DayMin)*(Hmax-Hmin)  + Hmin
                elif DayNo>=DayMax:
                    Hdes = (365-DayNo+DayMin)/(365-DayMax+DayMin)*(Hmax-Hmin) + Hmin
                elif DayNo<=DayMin:
                    Hdes = (DayMin-DayNo)/(365-DayMax+DayMin)*(Hmax-Hmin) + Hmin
                Vdes = (Hdes-Hmin)*V/(Hmax-Hmin)    
                Hdesign.append(Hdes)
                Vdesign.append(Vdes)
            Vdesign = np.append(np.nan,np.tile(Vdesign,num_yrs))
            Hdesign = np.append(np.nan,np.tile(Hdesign,num_yrs))
            
            s[0] = Vdesign[1]
            Hd_m[0] = s[0]/V*(Hmax-Hmin) + Hmin
            
        if opt_mtd=='mpc_rc':    
            r[0] = Q_in[1]
            
            if name in cascade['Lower'].values():
                cascade_no = [k for k in cascade['Lower'].keys() if cascade['Lower'][k] == name][0]
                # Find inflow from upper reservoir(s)
                for res_up in cascade['Upper'][cascade_no]:
                    Q_us = res_ops[res_up]['r'][0]
                r[0] = Q_in[1] + Q_us
            
        if opt_mtd=='mpc':
            Vdesign = np.nan * np.ones(shape = (len(Q_in),1)) #storage from forecasts
            Hdesign = np.nan * np.ones(shape = (len(Q_in),1)) #storage from forecasts
            
            s[0] = 0.7*V
            Hd_m[0] = s[0]/V*(Hmax-Hmin) + Hmin
            r[0] = Q_in[1]
            
            
        
        #------------------------------------------------------------------------------
        
        #for iterating reservoir op (feedback loop)
        #dictionary for 365 days, 24hrs/day
        # r_search = { i: {j:[] for j in range(1,25)} for i in range(1,len(Q_in)+1) }
        # Hd_m_search = { i: {j:[] for j in range(1,25)} for i in range(1,len(Q_in)+1) }
        # s_search = { i: {j:[] for j in range(1,25)} for i in range(1,len(Q_in)+1) }
        
        Qspil_hr = {i:[] for i in range(1,length)}
        
        # Minimum environmental flow
        for n in range(1,length):
            if Q_in[n] <= 0.4*Q_MAF_[n]:
                Q_MEF[n] = 0.6*Q_in[n]
            elif Q_in[n] > 0.8*Q_MAF_[n]:
                Q_MEF[n] = 0.3*Q_in[n]
            else: 
                Q_MEF[n] = 0.45*Q_in[n]
        
        res_ops[name] = {'Hdesign':Hdesign, 'Vdesign':Vdesign, 'Q_in':Q_raw, #'Q':Q,
                         'Q_MAF':Q_MAF_, 'Q_MEF':Q_MEF,
                          's':s, 'r':r, 'Hd_m':Hd_m, 'Qspil':Qspil, 
                          'availhydro': availhydro,
                          # 's_rc':s_rc, 'r_rc':r_rc, 'Hd_m_rc':Hd_m_rc, 'Qspil_rc':Qspil_rc, 
                          # 'availhydro_rc': availhydro_rc,                           
                          'Qspil_hr':Qspil_hr,'day_reop':day_reop,'numdays_reop':0} #'r_search':r_search, 
                            #'dispatchedhydro_final':dispatchedhydro_final,, 
                            #'dispatchedhydro':dispatchedhydro,'Hd_m_search':Hd_m_search, 's_search':s_search, 
        res_ops['cascade'] = cascade
        res_ops[name]['goal_HP'] = np.nan * np.ones(shape = (length,1))
        res_ops[name]['goal_finalV'] = np.nan * np.ones(shape = (length,1))
        
        if FC_typ == 'FC':
            res_ops[name]['Q_FC'] = Q_FC
            
Model1a(p1)        
