def get_var(day, name):
    
    hydhead = res_sysparam[name]['hydhead']
    Hmin = res_sysparam[name]['Hmin']
    Hmax = res_sysparam[name]['Hmax']
    V = res_sysparam[name]['V']
    Qmax = res_sysparam[name]['Release']
    Q_in = res_ops[name]['Q_in']       
    # Q = res_ops[name]['Q']

    s_ini = res_ops[name]['s'][day-1]
    r_ini = res_ops[name]['r'][day-1]
    
    
    Q_MEF = res_ops[name]['Q_MEF']
    # if (end-day)<p1:
    #     Q_MEF = Q_MEF[day:] #
    # else:
    Q_MEF = Q_MEF[day:day+p1]
    
    if FC_typ == 'perfect':
        # Q_forecast = Q_in.iloc[day-1,:]
        Q_forecast = Q_in[day:day+p1]
    
    # if FC_typ == 'perfect':
    #     Q_perfectFC  = res_ops[name]['Q_in']    
    #     if (end-day)<p1:
    #         Q_forecast = Q_perfectFC[day:] #
    #     else:
    #         Q_forecast = Q_perfectFC[day:day+p1]
            
    # if FC_typ == 'persistent':
    #     Q_perfectFC  = res_ops[name]['Q_in']   
    #     # Q_persistFC = pd.concat( [Q_perfectFC[-1:], Q_perfectFC[:-1] ]).reset_index(drop=True)
    #     Q_persistFC = np.append(Q_perfectFC[0],Q_perfectFC[-1])
    #     Q_persistFC = np.append(Q_persistFC,Q_perfectFC[1:-1])
    #     if (end-day)<p1:
    #         Q_forecast = Q_persistFC[day:] #
    #     else:
    #         Q_forecast = Q_persistFC[day:day+p1]
    
    elif FC_typ == 'clim':
        Q_MAF = res_ops[name]['Q_MAF']
        # if (end-day)<p1:
        #     Q_forecast = Q_MAF[day:] #
        # else:
        Q_forecast = Q_MAF[day:day+p1]
        
    elif FC_typ == 'FC':
        Q_forecast = res_ops[name]['Q_FC'].iloc[day-1,:]
            
    return hydhead, Hmin, Hmax, V, Qmax, Q_in, s_ini,r_ini, Q_MEF, Q_forecast




def mpc_Test(day, s_ini, r_ini, hydhead, Hmin, Hmax, V, Qmax, Q_MEF):            
    # % define the length of the prediction horizon, over which the optimal
    # % decision vector x must be optimized
    
    H_ = len(Q_forecast)
    H  = int(H_/num_res) #Prediction horizon
    
    # % initialize vectors
    r = np.nan * np.ones(shape = (H_,1))
    s = np.nan * np.ones(shape = (H_,1))
    Qspil = np.nan * np.ones(shape = (H_,1))
    
    for n in range(num_res):      
        # % define initial conditions and the integration time-step
        s[n*H]= s_ini[n]
        r[n*H]= r_ini[n]
    
    # % minimum and maximum value of the control variables
    max_u = np.repeat(Qmax,H-1)
    # min_u = Q_MEF
    
    
        # print(name+' aft add u/s (pred): ',Q_forecast)
        
    # ===========================================================================
    
    def optfun(x):               
    # % Nested function that computes the objective function
        for n in range(num_res):
            res = cascade_list[n]
            # print('optfun res ',res)
            
            Q_forecast_ = Q_forecast[n*H:(n+1)*H]
            
            if res in list(cascade['Lower'].values()):
                # print('Q_forecast_ top ',Q_forecast_)
                # print('Lower:',res)
                # cascade_no = [k for k in cascade['Lower'].keys() if cascade['Lower'][k] == res][0]
                # Find inflow from upper reservoir(s)
                # Q_us = 0
                for res_up in cascade['Upper'][cascade_no]:
                    # print('Upper:',res_up)
                    Q_us = res_ops[res_up]['mpc_r'] + res_ops[res_up]['mpc_Qspil']
                # print(name+' u/s (pred): ',Q_us)
                Q_forecast_ = np.reshape(Q_forecast_,(-1,1)) + Q_us    
                # print('Q_forecast_ bottom ',Q_forecast_)
            
            # % simulation
            for t in range (H-1):
                # print('t',t)
                # print('n*H+t',n*H+t)
                s2_new = s[n*H+t] + Q_forecast_[t+1]*24*3600
                r[n*H+t+1] = x[n*(H-1)+t]
                if abs(r[n*H+t+1]-r[n*H+t]) > 0.05*Qmax[n]:#*24*3600:
                    if r[n*H+t+1] > r[n*H+t]:
                        r[n*H+t+1] = r[n*H+t] + 0.05*Qmax[n] #*24*3600
                    else:
                        r[n*H+t+1] = r[n*H+t] - 0.05*Qmax[n] #*24*3600
                    
                Qspil[n*H+t+1] = max((s2_new-r[n*H+t+1]*24*3600-V[n])/24/3600,0)
                s[n*H+t+1] = s2_new - r[n*H+t+1]*24*3600 - Qspil[n*H+t+1] *24*3600
            
        # if name in cascade['Upper'].values():
        #     for t in range (H-1):       
        #         s2_new = s[t] + Q_forecast[t+1]*24*3600
        #         r[t+1] = x[t]
        #         Qspil[t+1] = max((s2_new-r[t+1]*24*3600-V)/24/3600,0)
        #         s[t+1] = s2_new - r[t+1]*24*3600 - Qspil[t+1] *24*3600
          
            res_ops[res]['mpc_r'] = r[n*H:(n+1)*H]
            res_ops[res]['mpc_s'] = s[n*H:(n+1)*H]
            res_ops[res]['mpc_Qspil'] = Qspil[n*H:(n+1)*H]
                
        # % compute the step-costs over the simulation horizon and the aggregated
        # % cost, which correspond to the objective function
        if opt_mtd == 'mpc':
            g = step_cost_2(day, s.flatten(), r.flatten(), Qspil.flatten())#, min_u, max_u)        
            f = np.mean(g)
        
        elif opt_mtd == 'mpc_rc':
            f = step_cost_2_rc(day,H, s.flatten(), r.flatten(), Qspil.flatten())
        return f
   
# =============================================================================
    x0 = np.repeat(r_ini,H-1)
    # print('x0',x0)
    
    A = np.zeros(shape=((H-1)*num_res*2,(H-1)*num_res))
    # print('A before',A)
    b = np.zeros(shape=(2*num_res*(H-1)))
    for n in range(num_res):
        A[2*n*(H-1):2*(n+1)*(H-1),n*(H-1):(n+1)*(H-1)] = np.concatenate((np.identity(H-1), -np.identity(H-1)), axis=0)
        
        max_u = res_sysparam[cascade_list[n]]['Release']
        # min_u = res_ops[cascade_list[n]]['Q_MEF'][day:day+p1]
        min_u = 0
        # b[2*n*(H-1):2*(n+1)*(H-1)] = np.concatenate((max_u*np.ones(shape = (H-1,)), -np.reshape(min_u,(H-1,))), axis=0)
        b[2*n*(H-1):2*(n+1)*(H-1)] = np.concatenate((max_u*np.ones(shape = (H-1,)), min_u*np.ones(shape = (H-1,))), axis=0)
    
    # print('A ',A)
    # else:
    #     A = np.concatenate((np.identity(H-1), -np.identity(H-1)), axis=0)
    #     b = np.concatenate((max_u*np.ones(shape = (H-1,)), -np.reshape(min_u,(H-1,))), axis=0)
    #     x0 = np.ones(shape = (H-1,))*r_ini
# =============================================================================    
    
    from scipy.optimize import LinearConstraint
    linear_constraint = LinearConstraint(A, -np.inf*np.ones(shape = b.shape), b)
    # linear_constraint = LinearConstraint(A, 0*np.ones(shape = (H-1,)), max_u*np.ones(shape = (H-1,)))
    
    # from scipy.optimize import Bounds
    # lb = -V + s[:-1].flatten() + Q_forecast[1:] - Qspil[1:].flatten()
    # ub = s[:-1].flatten() + Q_forecast[1:] - Qspil[1:].flatten()
    
    # % Determine the initialization vector
    # x0 = np.ones(shape = (H-1,))*r_ini#*0.5*V#min_u

    from scipy.optimize import minimize
    x = minimize(optfun, x0, args=(), method='trust-constr', 
             hess=None, hessp=None, bounds=None, constraints=linear_constraint, tol=1e-03, callback=None, 
             options={'maxiter': 300}) #bounds=Bounds(lb,ub)
    
    return x#,s,r

def step_cost_2_rc(day,H, s, r, Qspil):#, min_u, max_u    
    goal_HP = []
    goal_finalV = []
    # H = len(s)-1
    
    # pre-allocate the memory
    Hd_m = s.copy()
    availhydro_ = np.zeros (shape = (H,1))
    
    for n in range(num_res):
        res = cascade_list[n]
        # print('res ',res)
        
        Hd_m_ = s[n*H:(n+1)*H]/V[n]*(Hmax[n]-Hmin[n]) + Hmin[n]
        Hd_m[n*H:(n+1)*H] = Hd_m_
        # print('Hd_m',Hd_m)
        r_ = r[n*H:(n+1)*H]
        
    # Hd_m = s/V*(Hmax-Hmin) + Hmin
    # availhydro = np.zeros (shape = (H,1))
    # # hydropeak = np.zeros (shape = (H,1))
    
        for i in range(H-1):     
          # Available HP
          energy = turbine_factor*(hydhead[n] -(Hmax[n]-(Hd_m_[i]+Hd_m_[i+1])/2))*r_[i+1]*9.81/1000
          # energy=energy[0]
          if energy <0:
            energy=0
          availhydro_[i] = energy*24
      # # Hydropeaking
      # hydropeak[i] = max( abs(r[i+1]-r[i])-0.05*max_u ,0)     
       
    # # For normalization
    # # Hrange = Hmax-Hmin
    # HPrange = turbine_factor*(hydhead -(Hmax-(Hmax+Hmin)/2))*Qmax*9.81/1000 *24
        
        pot_HP = 24*res_sysparam[res]['plant_cap']
        goal_HP_ = np.mean( availhydro_.flatten()/ pot_HP )
        goal_HP = np.append(goal_HP,goal_HP_)
        # print(res,' plant_cap:',res_sysparam[res]['plant_cap'],', availHP:',availhydro_,
              # 'goal_HP:',goal_HP)
    
        # if (end-day)<p1:
        #     Vdes = res_ops[res]['Vdesign'][-1]
        # else:
        #     Vdes = res_ops[res]['Vdesign'][day+p1-1]+1 #To counter the case when Vdes=0
        
        ## Vdesign to be looped round the year. 
        # e.g. opt cost to go for day 337 for 30-day horizon to be Vdes on day 1##
        Vdes_day = (day+p1-1)%365
        if Vdes_day == 0:
            Vdes_day = 365
        
        Vdes = res_ops[res]['Vdesign'][Vdes_day]+1 #To counter the case when Vdes=0
        
        goal_finalV_ = abs(Vdes-res_ops[res]['mpc_s'].flatten()[-1]) /V[n]
        goal_finalV = np.append(goal_finalV,goal_finalV_)
        
        # print(res,'day ',day,'Vdes',Vdes)
        # print('goal_finalV',goal_finalV)
    
    # print('Vdes:',Vdes)
    
    # goal_flood = np.max(Qspil.flatten()[-1],0)
    # alpha = 1
    # beta = 100
    if name == 'KIRIROM1' or name == 'KIRIROM3':
        alpha=1; beta = 100
    elif name  == 'ATAY': #name == 'LRCHRUM' or 
        alpha=5; beta = 30 #alpha=1; beta = 30
        if res == 'LRCHRUM':
            alpha=20; beta = 5 #alpha=1; beta = 30
    elif name  == 'TATAY':
        alpha=1; beta = 50
    elif name  == 'KAMCHAY':
        alpha=1; beta=50
    # % aggregate step-costs:
    G = -alpha * np.mean(goal_HP) + beta * np.mean(goal_finalV) \
            
    # print(name,'goal_HP:',goal_HP,', goal_V:',goal_finalV)
    
    # res_ops[name]['goal_HP'][day] = goal_HP
    # res_ops[name]['goal_finalV'][day] = goal_finalV
    
    return G

# g=res_ops['KAMCHAY']['mpc_HP']+np.maximum( -test_s[1:] ,0)*1000
# =============================================================================
# # Forecast operations
# =============================================================================
# Assume forecast is mean annual flow
# Q_MAF = pd.read_csv('../Data/data_camb_Q_MAF.csv')
# Q_forecast_ = Q_MAF[day-1:day+p1-1]
# if (end-day)<p1:
#     Q_forecast_ = Q_MAF[day-1:] #

# Assume perfect forecast 
# Q_perfectFC  = pd.read_csv('../Data/Q_1000/data_camb_Q_syn_run'+str(run_no)+'.csv')
# Q_forecast_ = Q_perfectFC[day-1:day+p1-1]

# =============================================================================
cascade_list = []
for name in res_sysparam.keys():
    if name in cascade_list:
        continue
    # if name == 'LRCHRUM':
    #     haha = haha
    # elif name == 'TATAY':
    #     haha = hoho
    
        
    for k in list(cascade['Upper'].keys()):
        if name in list(cascade['Upper'][k]):
            # print(k)
            cascade_no = k
            cascade_list = np.append (list(cascade['Upper'][cascade_no]),cascade['Lower'][cascade_no] ) 
            num_res = len(cascade_list) #Find no. of reservoirs to build matrix later
            
    
    # if name in list(cascade['Upper'].values())[0]:
    #     cascade_no = [k for k in cascade['Upper'].keys() if cascade['Upper'][k] == name][0]
    #     cascade_list = cascade['Upper'][cascade_no],cascade['Lower'][cascade_no]    
    #     num_res = len(cascade_list) #Find no. of reservoirs to build matrix later
        else: 
            cascade_list = [name]
            num_res = 1
    
    # initialize values and TS
    Hmin=[]; Hmax=[]; hydhead=[]; V=[]; Qmax=[]; s_ini=[];r_ini=[];
    Q_forecast=[]; Q_actual=[]; Q_MEF=[];
    
    
    for res in cascade_list:
        hydhead_, Hmin_, Hmax_, V_, Qmax_, Q_in_, s_ini_,r_ini_, Q_MEF_, Q_forecast_ = get_var(day,res)
        hydhead= np.append(hydhead, hydhead_)
        Hmin= np.append(Hmin, Hmin_)
        Hmax= np.append(Hmax, Hmax_)
        V= np.append(V, V_)
        Qmax= np.append(Qmax, Qmax_)
        # Q_in= np.append(Q_in, Q_in_) ### To deal with during actual operations ###
        s_ini= np.append(s_ini, s_ini_)
        r_ini= np.append(r_ini, r_ini_)
        Q_MEF= np.append(Q_MEF, Q_MEF_)
        Q_forecast= np.append(Q_forecast, np.append(np.nan,Q_forecast_))
        
    ################################## HERE ################################################
    # % Determine the trajectory of the optimal controls
    x  = mpc_Test(day, s_ini, r_ini, hydhead, Hmin, Hmax, V, Qmax, Q_MEF).x
    # res_ops[name]['x'] = x
    #### Prediction horizon: Account for extra # days at the end ###
    if (end-day) < (2*p2):
        H = end-day+1
    else: 
        H = p2
    
# =============================================================================
#     # Actual operations
# =============================================================================

    for n in range(num_res):
        name = cascade_list[n]
        hydhead, Hmin, Hmax, V, Qmax, Q_in, s_ini,r_ini, Q_MEF, Q_forecast = get_var(day,name)
        Q_actual = Q_in[day:day+H] #One week actual Q
        
        if name in list(cascade['Lower'].values()):
        #     for res_up in cascade['Upper'][cascade_no]:
        #         print(res_up)
        
        # if name in cascade['Lower'].values():
        #     cascade_no = [k for k in cascade['Lower'].keys() if cascade['Lower'][k] == name][0]
            # Find inflow from upper reservoir(s)
            for res_up in cascade['Upper'][cascade_no]:
                # print (res_up)
                # print (str(res_ops[res_up]['r'][day]) +' '+ str(res_ops[res_up]['Qspil'][day]))
                Q_us = res_ops[res_up]['r'][day:day+H] + res_ops[res_up]['Qspil'][day:day+H]
            # print(name+' before add u/s: ',str(Q_in[day]))
            Q_actual = Q_actual + Q_us.flatten()    
        
        # r = x[:p2]
        Hd_m_ini = s_ini/V*(Hmax-Hmin) + Hmin    
        
        # Initialize
        # r_ = np.nan * np.ones(shape = (p2+1,1))
        s_ = np.nan * np.ones(shape = (H+1,1))
        r_ = np.nan * np.ones(shape = (H+1,1))
        Hd_m_ = np.nan * np.ones(shape = (H+1,1))
        Qspil_ = np.nan * np.ones(shape = (H+1,1))
        availhydro_ = np.nan * np.ones(shape = (H,1))
            
        # r_ = np.append(r_ini,x[:H]) # Release results from forecast
        # r_ = np.reshape(r_,(-1,1))
        s_[0] = s_ini
        r_[0] = r_ini
        Hd_m_[0] = Hd_m_ini
        Qspil_[0] = np.nan
    
        # if opt_mtd == 'mpc':
        #     Vmpc = res_ops[name]['mpc_s']
        #     Vdesign=res_ops[name]['Vdesign'];            Vdesign[day:day+H]=Vmpc[1:H+1]
        
        for t in range(H):
            s2_new = s_[t] + Q_actual[t]*24*3600
            # r[t+1] = x[t]
            
            if opt_mtd == 'mpc':
                if s2_new>Vmpc[t+1]:               
                    if (s2_new - Vmpc[t+1] <= Qmax*24*3600):
                        r_new = s2_new - Vmpc[t+1]
                    else:
                        r_new = Qmax*24*3600
                else:                        
                    r_new = 0
                    
            elif opt_mtd == 'mpc_rc': 
                # r_new = x[t]*24*3600
                r_new = res_ops[name]['mpc_r'][t+1]*24*3600
                # print(name,' day',day,'x=',x[t])
                # break
            Qdischarge_new = r_new/24/3600    
            # Account for hydropeaking
            if abs(Qdischarge_new-r_[t]) > 0.1*Qmax:#*24*3600:
                if Qdischarge_new > r_[t]:
                    Qdischarge_new = r_[t] + 0.1*Qmax#*24*3600
                else:
                    Qdischarge_new = r_[t] - 0.1*Qmax#*24*3600
            
            Qdischarge_new = max(Qdischarge_new, Q_MEF[t])
            Qdischarge_new = min(Qdischarge_new, Qmax)
            
            r_[t+1] = Qdischarge_new
            Qspil_[t+1] = max((s2_new-Qdischarge_new*24*3600-V)/24/3600,0)
            s_[t+1] = s2_new - Qdischarge_new*24*3600 - Qspil_[t+1] *24*3600
            if s_[t+1]<0: #Negative storage not allowed
                s_[t+1] = 0
                r_[t+1] = s2_new /24/3600
            Hd_m_[t+1] = s_[t+1]/V*(Hmax-Hmin) + Hmin        
                    
            energy = turbine_factor*(hydhead -(Hmax-(Hd_m_[t]+Hd_m_[t+1])/2))*r_[t+1]*9.81/1000
            availhydro_[t] = energy*24
                
    
        
        availhydro=res_ops[name]['availhydro'];  availhydro[day:day+H]=availhydro_
        s=res_ops[name]['s'];                    s[day:day+H]=s_[1:]
        r=res_ops[name]['r'];                    r[day:day+H]=r_[1:]
        Hd_m=res_ops[name]['Hd_m'];              Hd_m[day:day+H]=Hd_m_[1:]
        Qspil=res_ops[name]['Qspil'];            Qspil[day:day+H]=Qspil_[1:]    
                             
        res_ops[name].update(energy=energy, availhydro=availhydro, 
                             s=s, r=r, Hd_m=Hd_m, Qspil=Qspil)

    # print('finish res ',name)
# fig,ax = plt.subplots()
# ax.plot(res_ops[name]['Q_in'],label="Q_in")
# ax.plot(res_ops[name]['Q_MAF'],label="Q_MAF")
# ax.plot(res_ops[name]['Q_MEF'],label="Q_MEF")
# ax.plot(res_ops[name]['r'],label="r")
# # ax.legend(loc="upper left")
# ax.set_ylabel('[m3/s]')
# ax.legend(loc='upper left', bbox_to_anchor=(-0.3, 1))

# ax2 = ax.twinx()
# ax2.plot(res_ops[name]['Vdesign'], color='black',label="Vdesign")
# ax2.plot(res_ops[name]['s'], color='grey',Linestyle=':',label="S")
# ax2.plot(np.zeros(len(res_ops[name]['s'])), color='brown',Linestyle=':',label="s=0")
# ax2.legend(loc='upper left', bbox_to_anchor=(1.05, 1))
# ax2.set_ylabel('[m3]')

# ax.set_title(name)