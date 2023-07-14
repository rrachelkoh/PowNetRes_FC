import os
# os.chdir('E:\\SUTD2021Research\\CoSimulation\\Camb0')

# import pyomo.environ as pyo
from pyomo.opt import SolverFactory, SolverStatus
# from pyomo.core import Var
# from pyomo.core import Param
# from operator import itemgetter
import pandas as pd
from datetime import datetime, timedelta
# import pyomo.environ as pyo
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.dates as mdates

mem_num = 1

run_no = 0
FC_typ = 'FC' #'perfect','clim','FC'



# t_start = datetime.now()

# from ResModel1b import ResModel1b
from Solver import Solver
# =============================================================================
# Specify optimization method (or none)
# =============================================================================
opt_mtd = 'mpc_rc' #'mpc','none','mpc_rc'

###Segment C.1
yr_start = 2000
yr_end = 2018
num_yrs = yr_end-yr_start+1

FC_name = FC_typ+'_'+str(num_yrs)+'yrs'

# #Define the reservoirs that should be reoperated
res_reop = [] #"KAMCHAY","ATAY","LRCHRUM","TATAY","KIRIROM1","KIRIROM3"
cascade = { "Upper":{"1": {"ATAY"}}, "Lower":{"1": {"LRCHRUM"}} } 


# =============================================================================


# =============================================================================
# Ensure folder exists to store results
# =============================================================================

if len(res_reop) == 0:
    path = 'results/'+FC_name+'/mem'+str(mem_num)+'/'
else:
    path = 'results/'+FC_name+'/mem'+str(mem_num)+'_withFB/'


def assure_path_exists(path_):
    '''Creates directory if it doesn't exist'''
    dir = os.path.dirname(path_)
    if not os.path.exists(dir):
        os.makedirs(dir)
assure_path_exists(os.getcwd() + '/'+ path)

# =============================================================================

# yr=2016

# res_ops = {}
turbine_factor = 0.9
#Unit cost of generation / import of each fuel type
gen_cost = {'coal_st':5.2, 'oil_ic':6.0, 'oil_st':6.0,
             'imp_viet':65, 'imp_thai':66, 'slack':1000}
#Unit cost of hydro import
h_import_cost = 48

if opt_mtd == 'mpc' or opt_mtd == 'mpc_rc':
    # For MPC #
    p1=30 # Planning horizon
    p2=1 # Forecast horizon
else:
    p1 = 0
    
exec(open('ResModel1a.py').read())
# exec(open('mpc_setup.py').read())

# =============================================================================
###Segment C.2
# from .. import PowNetModel_Camb
from PowNetModel_Camb import model


###Segment C.3
if len(res_reop) == 0:
    dat_name = 'pownet_data_Camb_'+FC_name+'_mem'+str(mem_num)
else: 
    dat_name = 'pownet_data_Camb_'+FC_name+'_mem'+str(mem_num)+'_FB'


exec(open('DataSetup_Camb.py').read())
instance = model.create_instance(dat_name+'.dat')


###Segment C.4
opt = SolverFactory("gurobi") ##SolverFactory("cplex")
opt.options["threads"] = 1
opt.options["TimeLimit"] = 15 ##in seconds

H = instance.HorizonHours
K=range(1,H+1)
start = 1 ##1 to 364
end  = 365 * num_yrs



###Segment C.5
#Space to store results
on=[]
switch=[]

mwh=[]
hydro=[]
# solar=[]
# wind=[]

hydro_import=[]

srsv=[]
nrsv=[]

vlt_angle=[]

opt_status=[]
# =============================================================================

for day in range(start,end+1):
    if opt_mtd == 'mpc' or opt_mtd == 'mpc_rc':
        if day%p2==1 or p2==1:
            if (end-day+1) > p2:
                exec(open('mpc_try.py').read())

    
    for z in instance.d_nodes:
      #load Demand and Reserve time series data
        for i in K:
            instance.HorizonDemand[z,i] = instance.SimDemand[z,((day-1)*24)%8760+i]
            instance.HorizonReserves[i] = instance.SimReserves[((day-1)*24)%8760+i] 
    # print('Node ',z,' hr',i, 'Demand:',instance.SimDemand[z,((day-1)*24)%8760+i])
    
    # for z in instance.s_nodes:
    #   # load Solar time series data
    #     for i in K:
    #         instance.HorizonSolar[z,i] = instance.SimSolar[z,(day-1)*24+i]
    
    # for name in instance.h_nodes:
    # exec(open('ResModel1b.py').read())
    

    
    for z in instance.h_nodes:
     #load Hydropower time series data
        instance.plant_cap[z] = res_sysparam[z]['plant_cap']*turbine_factor
        instance.HorizonHydro[z] = res_ops[z]['availhydro'][day]
 
            
#     for z in instance.w_nodes:
#      #load Wind time series data
#         for i in K:
#             instance.HorizonWind[z,i] = instance.SimWind[z,(day-1)*24+i]
            
    for z in instance.h_imports:
      #load Hydropower time series data
        for i in K:
            instance.HorizonHydroImport[z,i] = instance.SimHydroImport[z,(day-1)*24+i]     

    hydro_,hydro_import_,vlt_angle_,switch_,on_,mwh_,nrsv_,srsv_ = Solver(instance,day)  #,solar_
# =============================================================================
#     result = opt.solve(instance) ##,tee=True to check number of variables
#     # instance.display()
#     if result.solver.status == SolverStatus.aborted: #max time limit reached #v1.3
#         result.solver.status = SolverStatus.warning #change status so that results can be loaded
#         opt_status.append(day)
#         
#     instance.solutions.load_from(result) 
#     # opt_obj.append(getattr(instance, 'SystemCost')()) 
#  
# #  #The following section is for storing and sorting results
#     for v in instance.component_objects(Var, active=True):
#         varobject = getattr(instance, str(v))
#         a=str(v)
#         if a=='hydro':   
#              hydro_ = []
#              for index in varobject:
#                  if int(index[1]>0 and index[1]<25):
#                     if index[0] in instance.h_nodes:
#                         # hydro.append((index[0],index[1]+((day-1)*24),varobject[index].value)) 
#                         hydro_.append((index[0],index[1]+((day-1)*24),varobject[index].value)) 
# 
#         if a=='solar':
#               solar_ = []          
#               for index in varobject:
#                   if int(index[1]>0 and index[1]<25):
#                     if index[0] in instance.s_nodes:
#                         # solar.append((index[0],index[1]+((day-1)*24),varobject[index].value)) 
#                         solar_.append((index[0],index[1]+((day-1)*24),varobject[index].value)) 
# 
# #         if a=='wind':
#       
# #              for index in varobject:
# #                  if int(index[1]>0 and index[1]<25):
# #                     if index[0] in instance.w_nodes:
# #                         wind.append((index[0],index[1]+((day-1)*24),varobject[index].value))   
#                             
# 
#         if a=='hydro_import':     
#               hydro_import_ = []
#               for index in varobject:
#                   if int(index[1]>0 and index[1]<25):
#                     # if index[0] in instance.h_imports:
#                         # hydro_import.append((index[0],index[1]+((day-1)*24),varobject[index].value))  
#                         hydro_import_.append((index[0],index[1]+((day-1)*24),varobject[index].value)) 
# 
#         if a=='vlt_angle':
#              vlt_angle_ = []
#              for index in varobject:
#                  if int(index[1]>0 and index[1]<25):
#                     if index[0] in instance.nodes:
#                         # vlt_angle.append((index[0],index[1]+((day-1)*24),varobject[index].value))   
#                         vlt_angle_.append((index[0],index[1]+((day-1)*24),varobject[index].value))
# 
#         if a=='mwh':  
#             mwh_ = []
#             ini_mwh_ = {} #v1.3
#             for index in varobject:
#                 if int(index[1]>0 and index[1]<25):
#                     # mwh.append((index[0],index[1]+((day-1)*24),varobject[index].value))  
#                     mwh_.append((index[0],index[1]+((day-1)*24),varobject[index].value))  
#                 if int(index[1])==24:
#                     ini_mwh_[index[0]] = varobject[index].value                            
# 
#                         
#         if a=='on':   
#             on_ = []
#             ini_on_ = {}  
#             for index in varobject:
#                 if int(index[1]>0 and index[1]<25):
#                     # on.append((index[0],index[1]+((day-1)*24),varobject[index].value))
#                     on_.append((index[0],index[1]+((day-1)*24),varobject[index].value))
#                 if int(index[1])==24:
#                     ini_on_[index[0]] = varobject[index].value    
# 
#         if a=='switch':  
#             switch_ = []
#             for index in varobject:
#                 if int(index[1]>0 and index[1]<25):
#                     # switch.append((index[0],index[1]+((day-1)*24),varobject[index].value))
#                     switch_.append((index[0],index[1]+((day-1)*24),varobject[index].value))
# 
#         if a=='srsv':   
#             srsv_ = []
#             for index in varobject:
#                 if int(index[1]>0 and index[1]<25):
#                     # srsv.append((index[0],index[1]+((day-1)*24),varobject[index].value))
#                     srsv_.append((index[0],index[1]+((day-1)*24),varobject[index].value))
#                         
#         if a=='nrsv':   
#             nrsv_ = []
#             for index in varobject:
#                 if int(index[1]>0 and index[1]<25):
#                     # nrsv.append((index[0],index[1]+((day-1)*24),varobject[index].value)) 
#                     nrsv_.append((index[0],index[1]+((day-1)*24),varobject[index].value))           
#     
#     # Update initialization values for "on" 
#     for z in instance.Generators:
#         instance.ini_on[z] = round(ini_on_[z])
#         instance.ini_mwh[z] = max(ini_mwh_[z],0) #v1.3
# =============================================================================
    
    # for name in instance.h_nodes:
    # repeat = False
    if len(res_reop) != 0:
        exec(open('ResModel2_d.py').read())
    
    # if repeat == True:
    #     continue
    
    
    # solar.extend(solar_)
    hydro.extend(hydro_)
    hydro_import.extend(hydro_import_)
    on.extend(on_)
    mwh.extend(mwh_)
    vlt_angle.extend(vlt_angle_)
    srsv.extend(srsv_)
    nrsv.extend(nrsv_)
    switch.extend(switch_)
    
    print('Member ', mem_num,'  complete day ',day)
    # day = day+1
    
    print(str(datetime.now()))
        
hydro_pd=pd.DataFrame(hydro,columns=('Node','Time','Value'))
hydro_import_pd=pd.DataFrame(hydro_import,columns=('Node','Time','Value'))
# solar_pd=pd.DataFrame(solar,columns=('Node','Time','Value'))
#     wind_pd=pd.DataFrame(wind,columns=('Node','Time','Value'))

vlt_angle_pd=pd.DataFrame(vlt_angle,columns=('Node','Time','Value'))

mwh_pd=pd.DataFrame(mwh,columns=('Generator','Time','Value'))    
on_pd=pd.DataFrame(on,columns=('Generator','Time','Value'))
switch_pd=pd.DataFrame(switch,columns=('Generator','Time','Value'))
srsv_pd=pd.DataFrame(srsv,columns=('Generator','Time','Value'))
nrsv_pd=pd.DataFrame(nrsv,columns=('Generator','Time','Value'))


# t_end = datetime.now()

# del a, v
###to save outputs
# out_path = 'results/'+FC_name+'/R'+str(run_no)+'/out_R'+str(run_no)
out_path = path+'out_mem'+str(mem_num)

mwh_pd.to_csv(out_path+'_mwh.csv')

hydro_pd.to_csv(out_path+'_hydro.csv')
hydro_import_pd.to_csv(out_path+'_hydro_import.csv')
# solar_pd.to_csv(out_path+'_solar.csv')
#     wind_pd.to_csv('out_syn_R'+str(run_no)+'_'+str(yr)+'_wind.csv')

vlt_angle_pd.to_csv(out_path+'_vlt_angle.csv')

on_pd.to_csv(out_path+'_on.csv')
switch_pd.to_csv(out_path+'_switch.csv')
srsv_pd.to_csv(out_path+'_srsv.csv')
nrsv_pd.to_csv(out_path+'_nrsv.csv')


##############################
# exec(open('ana_all.py').read())

def postprocess():
    #Calculate Co2 emissions
    
    num_days = num_yrs * 365
    day_range = np.arange(num_days)
    hr_range = np.arange(num_days*24)
    
    
    # =============================================================================
    #Total annual generation by fuel type
    # gen={'hydro':sum(hydro_pd['Value'])/1000, 'solar':sum(solar_pd['Value'])/1000}
    # df_gen_gwh = (mwh_pd.groupby(['type'])['Value'].sum())/1000
    # gen.update(df_gen_gwh)
    # print(gen)
    
    
    # Dictionary for dispatchable plant parameters
    gen_param = df_gen.set_index('name')
    gen_param = gen_param.to_dict()
    # Unit of gen_all: MWh
    gen_all = mwh_pd[['Generator','Time']].copy()
    
    for z in list(gen_param.keys()):
        gen_all[z] = gen_all['Generator'].map(gen_param[z])
    gen_all['mwh'] = mwh_pd.Value
    gen_all['on'] = on_pd.Value
    gen_all['switch'] = switch_pd.Value
    
    gen_all['fx_cost'] = gen_all.maxcap*gen_all.fix_om*gen_all.on
    gen_all['start_cost'] = gen_all.maxcap*gen_all.st_cost*gen_all.switch
    gen_all['vr_cost'] = gen_all.mwh*(gen_all.heat_rate*gen_all.gen_cost+gen_all.var_om)
    gen_all['total_cost'] = gen_all.fx_cost + gen_all.start_cost + gen_all.vr_cost #
    
    # =============================================================================
    #Generation mix
    # ==============
    # types = list(set(df_gen.typ)) #unique generation types
    
    ##---Hourly---##
    #Generations by fuel type and time (GWh)
    disp_hourly = gen_all[['Time','typ','mwh']].copy()
    disp_hourly = disp_hourly.groupby(['typ','Time'])['mwh'].sum()/1000
    ##Hourly TS of dispatchables
    disp_hourly = disp_hourly.unstack().T
    ##Hourly TS of hydro
    hydro_hourly = hydro_pd.groupby(['Time'])['Value'].sum()/1000
    hydro_import_hourly = hydro_import_pd.groupby(['Time'])['Value'].sum()/1000
    ##Hourly TS of solar
    # solar_hourly = solar_pd.groupby(['Time'])['Value'].sum()/1000
    
    
    ##---Daily---##
    ##Daily TS of dispatchables (GWh)
    disp_daily = disp_hourly.groupby(hr_range // 24).sum()
    ##Daily TS of hydro (in GWh)
    hydro_daily = hydro_hourly.groupby(hr_range // 24).sum()
    hydro_import_daily = hydro_import_hourly.groupby(hr_range // 24).sum()
    # solar_daily = solar_hourly.groupby(hr_range // 24).sum()
    
    genmix_daily = disp_daily.copy()
    genmix_daily['hydro'] = hydro_daily
    # genmix_daily['solar'] = solar_daily
    genmix_daily['imp_hydro'] = hydro_import_daily
    genmix_daily['coal'] = disp_daily['coal_st']
    genmix_daily['oil'] = disp_daily['oil_st']+disp_daily['oil_ic']
    genmix_daily['solar'] = np.zeros(len(genmix_daily))
    
    #####--- PLot daily gen mix ---- #
    # =============================================================================
    fig, (ax1) = plt.subplots(nrows=1,ncols=1)
    ax1.stackplot(day_range,genmix_daily['hydro'],genmix_daily['solar'],genmix_daily['coal_st'],genmix_daily['imp_hydro'],
                      genmix_daily['imp_thai'],genmix_daily['imp_viet'],genmix_daily['oil'],genmix_daily['slack'],
                      labels=['hydro','solar','coal','import(Laos)','import(Thailand)','import(Vietnam)','oil','slack'])
    ax1.set_ylabel('Electricity generation (GWh/day)')
    ax1.legend(loc=(0.04,1.01),ncol=3)
    
    # Set the locator
    locator = mdates.MonthLocator()  # every month
    # Specify the format - %b gives us Jan, Feb...
    fmt = mdates.DateFormatter('%b')
    X = plt.gca().xaxis
    X.set_major_locator(locator)
    # Specify formatter
    X.set_major_formatter(fmt)
    # =============================================================================
    
    # plt.savefig('dailygenmix'+str(run_no)+'.tiff', format='tiff',dpi=1200,bbox_inches = "tight")
    # plt.savefig('results/'+FC_name+'/R'+str(run_no)+'/dailygenmix_R'+str(run_no)+'.png', bbox_inches='tight', transparent=True, dpi=1200)
    
    # =============================================================================
    #Annual costs
    # ==============
    
    
    # global sys_cost (thermal plants + imports)
    sys_cost_byplant = gen_all.groupby('Generator')['total_cost'].sum() ##Unit: MWh
    sys_cost = (sys_cost_byplant.sum() - sys_cost_byplant.SLACK1)/(10**6) \
        + hydro_import_daily.sum()*h_import_cost / (10**3) ##Unit of hydro_import_daily: GWh
    
    # print('Total sys cost ($M/yr): ',sys_cost)
    
    # =============================================================================
    #CO2 emissions
    # ==============
    
    co2_coal_rate = 1.04  ##tonne/mwh
    co2_oil_rate = 0.73  ##tonne/mwh
    
    # ann_coal_mwh = gen['coal']*(10**3)
    # ann_oil_mwh = gen['oil']*(10**3)
    
    ann_coal_mwh = genmix_daily.coal_st.sum() * 1000
    ann_oil_mwh = ( genmix_daily.oil_st.sum() + genmix_daily.oil_ic.sum() ) * 1000
    
    ann_co2_coal = ann_coal_mwh * co2_coal_rate/(10**6) ##in megatonnes
    ann_co2_oil = ann_oil_mwh * co2_oil_rate/(10**6) ##in megatonnes
    
    ann_co2_Mtonnes = ann_co2_coal + ann_co2_oil
    
    # print('coal emi:',ann_co2_coal,', oil emi:',ann_co2_oil,', total emi (Mton):',ann_co2_Mtonnes)
    
    
    
    # =============================================================================
    # Operating costs
    # ==============
    
    #Annual gen mix
    
    ann_hydro = genmix_daily['hydro'].sum()
    ann_solar = genmix_daily['solar'].sum()
    ann_coal = genmix_daily['coal_st'].sum()
    ann_imprt_Laos = genmix_daily['imp_hydro'].sum() 
    ann_imprt_Thai = genmix_daily['imp_thai'].sum()
    ann_imprt_Viet = genmix_daily['imp_viet'].sum()
    ann_oil = genmix_daily['oil'].sum()
    ann_slack = genmix_daily['slack'].sum()
     
    # # =============================================================================
    # # Plots
    # # =============================================================================
    
    
    # =============================================================================
    # # #Analyse hydro dispatch for each reservoir
    hy_hr = pd.DataFrame()
    hy_d = pd.DataFrame()
    for n in hydro_pd.Node.unique():
        hy_hr[n] = hydro_pd[hydro_pd.Node==n].set_index('Time').Value
        hy_d[n] = hydro_pd[hydro_pd.Node==n].reset_index().Value.groupby(hr_range//24).sum()   
    # # =============================================================================
    
    # #Hydro dispatch by HP plant
    h_all = pd.DataFrame()
    for h in instance.h_nodes:
        h_all[h] = hydro_pd[hydro_pd.Node==h].reset_index().Value
    h_all['total'] = h_all.sum(axis=1)
    # h_all_d = h_all.groupby(hr_range // 24).sum()
    
    availhydro_all = pd.DataFrame()
    for h in instance.h_nodes:
        availhydro_all[h] = res_ops[h]['availhydro'][1:].flatten()
    availhydro_all = availhydro_all
    
    
    for name in instance.h_nodes:
        s = res_ops[name]['s']; #s_rc = res_ops[name]['s_rc']
        r = res_ops[name]['r']; #r_rc = res_ops[name]['r_rc']
        availhydro = res_ops[name]['availhydro']; #availhydro_rc = res_ops[name]['availhydro_rc']
        # hyd_dispatched = h_all_d[name]
        Vdesign = res_ops[name]['Vdesign']
        Q_in = res_ops[name]['Q_in']
        Q_MEF = res_ops[name]['Q_MEF']
        Qspil = res_ops[name]['Qspil']
        day_reop = res_ops[name]['day_reop']
        
        res_all = [Vdesign[1:].flatten()/(10**6),s[1:].flatten()/(10**6),
                    Q_in[1:len(r)],Q_MEF[1:].flatten(),r[1:].flatten(),
                    availhydro[1:].flatten(),np.array(hy_d[name]).flatten(),
                    Qspil[1:].flatten(),day_reop[1:].flatten()]#
        
        res_all = pd.DataFrame(res_all).T
        res_all.columns=['Vdesign','s','Qin','Q_MEF','r','availhydro','dispatched','spill','day_reop']#
        # res_all.columns=['s_rc','s','r_rc','r','availhydro_rc','availhydro','dispatched']
        
        
        #WRITE RESERVOIR RESULTS TO EXCEL#
        # if run_no ==1:
        res_all.to_csv(path+name+'_mem'+str(mem_num)+'.csv',index=False)
        
        
    
    # =============================================================================
    # Check mindown constraint for coal plants
    # coal_plants = gen_all[gen_all.typ.str[:3]=='coa']
        
    # coal_mwh = pd.DataFrame()
    # for n in coal_plants.Generator.unique():
    #     coal_mwh[n] = coal_plants[coal_plants.Generator==n].reset_index().mwh
    #     coal_ana = coal_mwh[n]
    #     # coal_ana[coal_ana<0.001] = 0
    #     len_zeros = []
    #     len_nonzeros = []
    #     for k, v in coal_ana[coal_ana == 0].groupby((coal_ana != 0).cumsum()):
    #         len_zeros.append(len(v))
    #     if len_zeros == []:
    #         print('min down time check:', n, 'NA')
    #     else:
    #         print('min down time check:', n, min(len_zeros), ',occurrence:',k)
    #     for k, v in coal_ana[coal_ana != 0].groupby((coal_ana == 0).cumsum()):
    #         len_nonzeros.append(len(v))
    #     print('min up time check:',n, min(len_nonzeros), ',occurrence:',k)
    
    # runtime = (t_end-t_start)/ timedelta(minutes=1)
    
    
    
    ######### RESULTS TO COMPARE SCENARIOS #############
    
    #Check for non-empty rows in lists
    #i is day no., 1 is the 1st hr. If 1st hr non-empty, reoperation happened.
    # print('No. of days reoperated')     
    numdays_reop = []   
    for name in res_sysparam.keys():
        numdays_reop.append(res_ops[name]['numdays_reop'])
        # reop=0
        # for i in range(1,day+1):
        #     if res_ops[name]['r_search'][i][1]:
        #         reops = reop + 1
        # print(name,':', reop)
        # reops.append(reop)
    #############################
        
    res = ['Run'+str(run_no),ann_hydro,ann_solar,ann_coal,ann_imprt_Laos,ann_imprt_Thai,ann_imprt_Viet,ann_oil,ann_slack,
           sys_cost,ann_co2_coal,ann_co2_oil, 
           availhydro_all.sum().sum()/1000,hy_d.sum().sum()/1000] #, runtime
    res.extend(numdays_reop)
    # if run_no==1: 
    results = pd.DataFrame(res).T
    results.columns = ['Run','hydro','solar','coal','imp_Laos','imp_Thai','imp_Viet',
                       'oil','slack','cost','co2_coal','co2_oil',
                       'hyd_avail','hyd_dispatched',
                        'reopKAMCHAY','reopATAY','reopLRCHRUM','reopTATAY','reopKIRIROM1','reopKIRIROM3'] #'runtime',
    # else:
    #     results = pd.read_csv('results_all.csv')
    #     results.loc[len(results)] = res
    results.to_csv(path +'results_all.csv', index=False)

    ### GENMIX ###
    # if run_no == 1:
    # genmix_daily.to_excel('results/R'+str(run_no)+'/genmix_all.xlsx',sheet_name="run"+str(run_no), index=False)
    # else:
    #     with pd.ExcelWriter('genmix_all.xlsx', engine='openpyxl', mode='a') as writer: 
    #         genmix_daily.to_excel(writer,sheet_name="run"+str(run_no), index=False)
    
    ### HYDROPOWER USAGE ###
    # if run_no == 1:
    # availhydro_all.to_excel('results/R'+str(run_no)+'/HP_all.xlsx',sheet_name="run"+str(run_no)+"avail", index=False)
    # with pd.ExcelWriter('results/R'+str(run_no)+'/HP_all.xlsx', engine='openpyxl', mode='a') as writer: 
    #     hy_d.to_excel(writer,sheet_name="run"+str(run_no)+"used", index=False)
    # else:
    #     with pd.ExcelWriter('HP_all.xlsx', engine='openpyxl', mode='a') as writer: 
    #         availhydro_all.to_excel(writer,sheet_name="run"+str(run_no)+"avail", index=False)
    #         hy_d.to_excel(writer,sheet_name="run"+str(run_no)+"used", index=False)
            


postprocess()

# del day,end,i,index,instance,K,model,num_yrs,opt,start,varobject,yr,yr_start,yr_end,z,t_start,t_end,result,turbine_factor,tol
