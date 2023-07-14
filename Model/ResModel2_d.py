# dispatched = hydro_.copy()
# dispatched_pd=pd.DataFrame(dispatched,columns=('Node','Time','Value'))
# tol = 1e-4
# turbine_factor = turbine_factor

def ResModel2(day,res_sysparam,res_ops,hydro_,cascade):
    # cascade = res_ops['cascade']
    
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
            # print('reoperating ',name)
            
            hydro_dispatched = dispatched_pd[dispatched_pd.Node==name][-24:]
            hydro_tot_dispatched = sum(hydro_dispatched.Value)
            # hydro_tot_avail = res_ops[name]['energy']*24
                    
            # print(name, '  Avail HP:',hydro_tot_avail,'; Dispatched HP:',hydro_tot_dispatched)
            # print('diff', hydro_tot_avail - hydro_dispatched)
            
            # if (hydro_tot_dispatched == hydro_tot_avail): #sum(hydro_tot_avail)):#( hydro_tot_avail - hydro_dispatched ).all() == 0 : #
            #     print('All avail hydro dispatched')
            
            ### Update HP availability of downstream reservoirs in case upstream re-operated ###
            # if name in cascade['Lower'].values():  
            if name in lower_reop:  
                # print(name+' before (ResModel2):',res_ops[name]['energy']*24)
                from ResModel1b import ResModel1b
                res_ops = ResModel1b(name,day,res_sysparam,res_ops,turbine_factor,cascade)
            
            hydro_tot_avail = res_ops[name]['availhydro'][day]
            # print(name+': avail ',hydro_tot_avail,'dispatched ', hydro_tot_dispatched)
            
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
                    # print(name+': input to solver ',instance.HorizonHydro[name])
                    hydro_,_,_,_,_,_,_,_ = Solver(instance,day)
                    # print('hydro_',hydro_)
                    
                    should_restart = True
                    # day=day-1
                    # repeat = True
                    break
                else:
                    raise Exception('ERROR -- Avail hydro < dispatched!!')
                    print(name,'avail:',hydro_tot_avail,'dispatched:',hydro_tot_dispatched)
            
            ### Case 2: Avail > dispatched --> Over-supply ###
            elif (hydro_tot_avail - hydro_tot_dispatched >= tol ):#( hydro_tot_avail - hydro_dispatched ).any() > tol : #
                    print('Reoperate ',name)
                    
                    V = res_sysparam[name]['V']; hydhead = res_sysparam[name]['hydhead']
                    Hmax = res_sysparam[name]['Hmax']; Hmin = res_sysparam[name]['Hmin']
                    Q_in = res_ops[name]['Q_in'].copy()
                    Q_MEF = res_ops[name]['Q_MEF']
                
                    s = res_ops[name]['s']; #s_rc = res_ops[name]['s_rc']
                    r = res_ops[name]['r']
                    Hd_m = res_ops[name]['Hd_m']
                    Qspil = res_ops[name]['Qspil']
                    availhydro = res_ops[name]['availhydro']
                    # dispatchedhydro = res_ops[name]['dispatchedhydro']
                    # dispatchedhydro_final = res_ops[name]['dispatchedhydro_final']
                    # r_search = res_ops[name]['r_search']; #s_search = res_ops[name]['s_search'];Hd_m_search = res_ops[name]['Hd_m_search']; 
                    # Qspil_hr = res_ops[name]['Qspil_hr']
                    numdays_reop = res_ops[name]['numdays_reop']     
                    day_reop = res_ops[name]['day_reop']
                    
                    # dispatchedhydro[day-1] = hydro_tot_dispatched
                    # dispatchedhydro_final[day-1] = hydro_tot_dispatched
                
                    # s_hr = np.nan*np.ones(shape=(25,1))
                    # s_hr[0] = s[day-1]
                    # r_hr_final = [] #optimized hourly release
                    # Hd_m_hr_final = [] #optimized hourly WL
                    # # dispatchedHP=[]
                    # Q_spil_hr = [] #hourly spills
                    # Q_in_hr = Q_in[day].copy()
                    
                    
                    # if name == 'LRCHRUM':
                    #     Q_in_hr = Q_in_hr + res_ops['ATAY']['r'][day] + res_ops['ATAY']['Qspil'][day]
                    
                    if name in cascade['Lower'].values():
                        cascade_no = [k for k in cascade['Lower'].keys() if cascade['Lower'][k] == name][0]
                        # Find inflow from upper reservoir(s)
                        Q_us = 0
                        for res_up in cascade['Upper'][cascade_no]:
                            print ('Re-op: upper res ',res_up)
                            print (str(res_ops[res_up]['r'][day]) +' '+ str(res_ops[res_up]['Qspil'][day]))
                            Q_us = Q_us + res_ops[res_up]['r'][day] + res_ops[res_up]['Qspil'][day]
                        # print(name+' before: ',str(Q_in[day]))
                        Q_in[day] = Q_in[day] + Q_us
                        # print('after: ',Q_in[day])
                    
                    
                    max_it = 10
                    hydro_dispatched = hydro_dispatched.reset_index().drop(columns = ['index'])
                    
                    # s2_new = s[day-1] + Q_in[day]*24*3600
                    # Q_spil_new = max( (s2_new-V)/24/3600, 0 )
                    
                    goal = hydro_tot_dispatched/24                    
                    k = goal/9.81/turbine_factor*1000
                    
                    diff_r = np.inf
                    diff_E = -1
                    count = 0
                    # r_it = []
                    
                    # First assume that WL remains the same
                    Hd_m_new = Hd_m[day-1].copy()
                    
                    while (diff_r >= tol) or (diff_E<0) and (count<=max_it):
                        count = count+1
                        
                        r_new = k/(hydhead-(Hmax-(Hd_m[day-1]+Hd_m_new)/2))
                        r_new = max( Q_MEF[day], r_new )
                        
                        Q_spil_new = max( (s[day-1] + (Q_in[day]-r_new)*24*3600 - V)/24/3600, 0 )
                        
                        s_new = s[day-1] + (Q_in[day] - Q_spil_new - r_new)*24*3600
                        Hd_m_new = s_new/V*(Hmax-Hmin) + Hmin
                        
                        # print('r:',r_new,',s:',s_new)
                        # r_it.append(r_new)
                        
                        if count == 1:
                            diff_r = np.inf
                        
                        # 
                        else:
                            diff_r = abs(r_new - r_it)
                            availHP_new = turbine_factor*(hydhead - (Hmax-(Hd_m[day-1]+Hd_m_new)/2))*r_new*9.81/1000
                            diff_E = availHP_new[0] - goal
                            # print('goal:',goal, ',HPnew:',availHP_new[0])
                    
                        r_it = r_new
                    
                    # print(name,' r before:',r[day],', r aft:',r_new) #average over 24 hours or MEF
                    # Account for hydropeaking
                    if abs(r_new-r[day-1]) > 0.1*Qmax:#*24*3600:
                        if r_new > r[day-1]:
                            r_new = r[day-1] + 0.1*Qmax#*24*3600
                        else:
                            r_new = r[day-1] - 0.1*Qmax#*24*3600
                    
                    
                    r[day] = r_new #max(Q_MEF[day],
                    s[day] = s_new #s[day-1] + (Q_in[day]*24 - Q_spil_new - r[day]*24)*3600#s_hr[-1][0]   #s_hr[0] + (Q_in[day-1]*24 - sum(Q_spil_hr) - sum(r_hr_final))*3600
                    Hd_m[day] = Hd_m_new #s[day]/V*(Hmax-Hmin) + Hmin #Hd_m_new_hr
                    Qspil[day] = Q_spil_new
                    availhydro[day] = availHP_new[0]*24
                    # Qspil[day] = sum(Q_spil_hr)/24 #average over 24 hours
                    # dispatchedhydro_final[day] = sum(dispatchedHP)/24 #average over 24 hours
                    numdays_reop = numdays_reop+1
                    day_reop[day] = 1
                    # print('Storage by rule curve', s_rc[day], ', aft reop:', s[day])   
                        
                    res_ops[name].update(s=s, r=r, Hd_m=Hd_m, Qspil=Qspil, availhydro=availhydro,
                                         day_reop=day_reop,numdays_reop=numdays_reop)
                                 # dispatchedhydro=dispatchedhydro,r_search=r_search, s_search=s_search
                                 # dispatchedhydro_final=dispatchedhydro_final,Hd_m_search=Hd_m_search,
# =============================================================================
### Removed for forecast run
#                     if name in sum(list(cascade['Upper'].values()),[]):
#                         lower_reop_ = [k for k in cascade['Upper'].keys() if name in cascade['Upper'][k]][0]
#                         lower_reop.append(cascade['Lower'][lower_reop_])
# =============================================================================
                        # print('Upper res:',name,'  To reop:',lower_reop)
                        # print('Lower reop:',lower_reop)
                    #     # print('ATAY: avail ',hydro_tot_avail,'dispatched ', hydro_tot_dispatched)
                    #     from ResModel1b import ResModel1b
                    #     res_ops = ResModel1b('LRCHRUM',day,res_sysparam,res_ops,turbine_factor)
                        # print('LRCHRUM avail aft op ATAY:',res_ops['LRCHRUM']['energy']*24)
                    
    return res_ops
res_ops = ResModel2(day,res_sysparam,res_ops,hydro_,cascade)
