def ResModel2(day,res_sysparam,res_ops,hydro_,cascade):
    
    tol = 1e-4
    should_restart = True
    
    while should_restart:
        should_restart = False
        lower_reop = []
        
        for name in res_reop: 
            Qmax = res_sysparam[name]['Release']
        # =============================================================================
        # ### For each reservoir: Compare avail & dispatched HP ###  
        # =============================================================================
            dispatched = hydro_.copy()
            
            dispatched_pd=pd.DataFrame(dispatched,columns=('Node','Time','Value'))
            
            hydro_dispatched = dispatched_pd[dispatched_pd.Node==name][-24:]
            hydro_tot_dispatched = sum(hydro_dispatched.Value)
            
            if name in lower_reop:  
                from ResModel1b import ResModel1b
                res_ops = ResModel1b(name,day,res_sysparam,res_ops,turbine_factor,cascade)
            
            hydro_tot_avail = res_ops[name]['availhydro'][day]
            
            # =============================================================================
            # ### SCENARIOS ###  
            # =============================================================================
            
            ### Case 1: Fully dispatched ###
            if (abs(hydro_tot_avail - hydro_tot_dispatched) < tol ):
                print(name+' fully dispatched')
            
            ### Case 3: Dispatched > avail --> Under-supply in lower res ###
            elif (hydro_tot_dispatched - hydro_tot_avail >= tol):#sum(hydro_tot_avail) >= tol):#round( hydro_tot_avail - hydro_dispatched ).any() < 0 : #
                if name in cascade['Lower'].values():
                    print('Repeat day ',day)
                    
                    from Solver import Solver
                    instance.HorizonHydro[name] = hydro_tot_avail
                    hydro_,_,_,_,_,_,_,_ = Solver(instance,day)
                    
                    should_restart = True
                    break
                else:
                    raise Exception('ERROR -- Avail hydro < dispatched!!')
                    print(name,'avail:',hydro_tot_avail,'dispatched:',hydro_tot_dispatched)
            
            ### Case 2: Avail > dispatched --> Over-supply ###
            elif (hydro_tot_avail - hydro_tot_dispatched >= tol ):
                    print('Reoperate ',name)
                    
                    V = res_sysparam[name]['V']; hydhead = res_sysparam[name]['hydhead']
                    Hmax = res_sysparam[name]['Hmax']; Hmin = res_sysparam[name]['Hmin']
                    Q_in = res_ops[name]['Q_in'].copy()
                    Q_MEF = res_ops[name]['Q_MEF']
                
                    s = res_ops[name]['s']; 
                    r = res_ops[name]['r']
                    Hd_m = res_ops[name]['Hd_m']
                    Qspil = res_ops[name]['Qspil']
                    availhydro = res_ops[name]['availhydro']
                    numdays_reop = res_ops[name]['numdays_reop']     
                    day_reop = res_ops[name]['day_reop']
                    
                    
                    if name in cascade['Lower'].values():
                        cascade_no = [k for k in cascade['Lower'].keys() if cascade['Lower'][k] == name][0]
                        # Find inflow from upper reservoir(s)
                        Q_us = 0
                        for res_up in cascade['Upper'][cascade_no]:
                            print ('Re-op: upper res ',res_up)
                            print (str(res_ops[res_up]['r'][day]) +' '+ str(res_ops[res_up]['Qspil'][day]))
                            Q_us = Q_us + res_ops[res_up]['r'][day] + res_ops[res_up]['Qspil'][day]
                        
                        Q_in[day] = Q_in[day] + Q_us
                    
                    
                    max_it = 10
                    hydro_dispatched = hydro_dispatched.reset_index().drop(columns = ['index'])
                    
                    
                    goal = hydro_tot_dispatched/24                    
                    k = goal/9.81/turbine_factor*1000
                    
                    diff_r = np.inf
                    diff_E = -1
                    count = 0
                    
                    # First assume that WL remains the same
                    Hd_m_new = Hd_m[day-1].copy()
                    
                    while (diff_r >= tol) or (diff_E<0) and (count<=max_it):
                        count = count+1
                        
                        r_new = k/(hydhead-(Hmax-(Hd_m[day-1]+Hd_m_new)/2))
                        r_new = max( Q_MEF[day], r_new )
                        
                        Q_spil_new = max( (s[day-1] + (Q_in[day]-r_new)*24*3600 - V)/24/3600, 0 )
                        
                        s_new = s[day-1] + (Q_in[day] - Q_spil_new - r_new)*24*3600
                        Hd_m_new = s_new/V*(Hmax-Hmin) + Hmin
                                                
                        if count == 1:
                            diff_r = np.inf
                        
                        # 
                        else:
                            diff_r = abs(r_new - r_it)
                            availHP_new = turbine_factor*(hydhead - (Hmax-(Hd_m[day-1]+Hd_m_new)/2))*r_new*9.81/1000
                            diff_E = availHP_new[0] - goal
                    
                        r_it = r_new
                    
                    # Account for hydropeaking
                    if abs(r_new-r[day-1]) > 0.1*Qmax:
                        if r_new > r[day-1]:
                            r_new = r[day-1] + 0.1*Qmax
                        else:
                            r_new = r[day-1] - 0.1*Qmax
                    
                    
                    r[day] = r_new 
                    s[day] = s_new 
                    Hd_m[day] = Hd_m_new 
                    Qspil[day] = Q_spil_new
                    availhydro[day] = availHP_new[0]*24
                    numdays_reop = numdays_reop+1
                    day_reop[day] = 1
                                            
                    res_ops[name].update(s=s, r=r, Hd_m=Hd_m, Qspil=Qspil, availhydro=availhydro,
                                         day_reop=day_reop,numdays_reop=numdays_reop)
    return res_ops
res_ops = ResModel2(day,res_sysparam,res_ops,hydro_,cascade)
