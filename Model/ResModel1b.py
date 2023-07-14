# import numpy as np
# =============================================================================
def ResModel1b(name,day,res_sysparam,res_ops,turbine_factor):
    # for name in res_sysparam.keys(): 
        
    V = res_sysparam[name]['V']
    Vdesign = res_ops[name]['Vdesign']
    hydhead = res_sysparam[name]['hydhead']
    Hmax = res_sysparam[name]['Hmax']
    Hmin = res_sysparam[name]['Hmin']
    Qmax = res_sysparam[name]['Release']
    Q_in = res_ops[name]['Q_in'].copy()
    Q_MEF = res_ops[name]['Q_MEF']
    
    s = res_ops[name]['s']; #s_rc = res_ops[name]['s_rc']
    r = res_ops[name]['r']; #r_rc = res_ops[name]['r_rc']
    availhydro = res_ops[name]['availhydro']; #availhydro_rc = res_ops[name]['availhydro_rc']
    Hd_m = res_ops[name]['Hd_m']; #Hd_m_rc = res_ops[name]['Hd_m_rc']
    Qspil = res_ops[name]['Qspil']; #Qspil_rc = res_ops[name]['Qspil_rc']
    
    cascade = res_ops['cascade']
     
    # ==========Actual operations=========== #
    # for t in range(day,day+1): 
    # if name == 'LRCHRUM':    
    if name in cascade['Lower'].values():
        cascade_no = [k for k in cascade['Lower'].keys() if cascade['Lower'][k] == name][0]
        # Find inflow from upper reservoir(s)
        Q_us = 0
        for res_up in cascade['Upper'][cascade_no]:
            # print (res_up)
            # print (str(res_ops[res_up]['r'][day]) +' '+ str(res_ops[res_up]['Qspil'][day]))
            Q_us = Q_us + res_ops[res_up]['r'][day] + res_ops[res_up]['Qspil'][day]
        # print(name+' before add u/s: ',str(Q_in[day]))
        Q_in[day] = Q_in[day] + Q_us
        # print('after: ',str(Q_in[day]))
    
    #convert storage to water level
    # Hd_m_new = s[t]/V*(Hmax-Hmin) + Hmin
    #storage+inflow
    s2_new = s[day-1] + Q_in[day]*24*3600
    #convert storage to water level
    # H2_new = s2_new/V*(Hmax-Hmin) + Hmin
    #Calculate release
    if s2_new>Vdesign[day]:               
        if (s2_new - Vdesign[day] <= Qmax*24*3600):
            r_new = s2_new - Vdesign[day]
        else:
            r_new = Qmax*24*3600
    else:                        
        r_new = 0
    
    
    
    # if r_new < Q_MEF_:
    #     r_new = Q_MEF_
    
    #discharge flow rate
    # print(name,' r:',r_new/24/3600, ',Q_MEF:',Q_MEF[day])
    
    
    Qdischarge_new = r_new/24/3600
    # Account for hydropeaking
    if abs(Qdischarge_new-r[day-1]) > 0.1*Qmax:#*24*3600:
        if Qdischarge_new > r[day-1]:
            Qdischarge_new = r[day-1] + 0.1*Qmax#*24*3600
        else:
            Qdischarge_new = r[day-1] - 0.1*Qmax#*24*3600
    # Account for MEF
    Qdischarge_new = max(Qdischarge_new, Q_MEF[day])
    Qdischarge_new = min(Qdischarge_new, Qmax)
    # print('Qdischarge:',Qdischarge_new)
    #spillage when volume > reservoir capacity
    Qspil_new = max((s2_new-Qdischarge_new*24*3600-V)/24/3600,0)
    #mass balance : storage for next time step
    s_new = s[day-1] + (Q_in[day] - Qdischarge_new - Qspil_new )*24*3600
    
    if s_new<0: #Negative storage not allowed
        s_new = 0
        Qdischarge_new = s2_new /24/3600
    
    s[day] = s_new
    r[day] = Qdischarge_new
    Hd_m[day] = s_new/V*(Hmax-Hmin) + Hmin
    
    Qspil[day] = Qspil_new

# for t in range(day,day+1): 
    energy = turbine_factor*(hydhead -(Hmax-(Hd_m[day]+Hd_m[day-1])/2))*r[day]*9.81/1000
    energy=energy[0]
    if energy <0:
        energy=0
    availhydro[day] = energy *24
    
    # if name == 'LRCHRUM':
    #     print('Day',day,', LRCHRUM avail energy:', energy*24, ',Storage=',s_new)
    
    # =============================================================================
     
    # =============================================================================
# =============================================================================
#         # # ==========Rule curve operations=========== #
# =============================================================================
#     if name == 'LRCHRUM':
#         Q_in[day] = Q_in[day] + res_ops['ATAY']['r_rc'][day] + res_ops['ATAY']['Qspil_rc'][day]
#         
#     
#     #convert storage to water level
#     # Hd_m_new = s[t]/V*(Hmax-Hmin) + Hmin
#     #storage+inflow
#     s2_new = s[day-1] + Q_in[day]*24*3600
#     #convert storage to water level
#     # H2_new = s2_new/V*(Hmax-Hmin) + Hmin
#     #Calculate release
#     if s2_new>Vdesign[day]:               
#         if (s2_new - Vdesign[day] <= Qmax*24*3600):
#             r_new = s2_new - Vdesign[day]
#         else:
#             r_new = Qmax*24*3600
#     else:                        
#         r_new = 0
#     #discharge flow rate
#     Qdischarge_new = r_new/24/3600
#     #spillage when volume > reservoir capacity
#     Qspil_new = max((s2_new-V)/24/3600,0)
#     #mass balance : storage for next time step
#     s_new = s[day-1] - (Qdischarge_new + Qspil_new - Q_in[day])*3600*24
#     s[day] = s_new
#     r[day] = Qdischarge_new
#     Hd_m[day] = s_new/V*(Hmax-Hmin) + Hmin
#     
#     Qspil[day] = Qspil_new
# 
# # for t in range(day,day+1): 
#     energy = turbine_factor*(hydhead -(Hmax-(Hd_m[day]+Hd_m[day-1])/2))*r[day]*9.81/1000
#     energy=energy[0]
#     if energy <0:
#         energy=0
#     availhydro[day] = energy
# =============================================================================
#         Q_in = list(res_ops[name]['Q_in'])
#         
#         for t in range(day-1,day+1): 
#             if name == 'LRCHRUM':
#                 Q_in[t] = Q_in[t] + res_ops['ATAY']['r_rc'][t+1]
#                 Q_in_1a = Q_in.copy()
#             #convert storage to water level
#             Hd_m_new_rc = s_rc[t]/V*(Hmax-Hmin) + Hmin
#             #storage+inflow
#             s2_new_rc = s_rc[t] + Q_in[t]*24*3600
#             #convert storage to water level
#             H2_new_rc = s2_new_rc/V*(Hmax-Hmin) + Hmin
#             #Calculate release
#             if s2_new_rc>Vdesign[t]:               
#                 if (s2_new_rc - Vdesign[t] <= Qmax*24*3600):
#                     r_new_rc = s2_new_rc - Vdesign[t]
#                 else:
#                     r_new_rc = Qmax*24*3600
#             else:                        
#                 r_new_rc = 0
#             #discharge flow rate
#             Qdischarge_new_rc = r_new_rc/24/3600
#             #spillage when volume > reservoir capacity
#             Qspil_new_rc = max((s2_new_rc-V)/24/3600,0)
#             #mass balance : storage for next time step
#             s_new_rc = s_rc[t] - (Qdischarge_new_rc + Qspil_new_rc - Q_in[t])*3600*24
#                 
#             Qspil_rc[t] = Qspil_new_rc
#             s_rc[t+1] = s_new_rc
#             r_rc[t+1] = Qdischarge_new_rc
#             Hd_m_rc[t+1] = Hd_m_new_rc
#         
#         for t in range(day,day+1): 
#             energy_rc = turbine_factor*(hydhead-(Hmax-(Hd_m_rc[t]+Hd_m_rc[t+1])/2))*r_rc[t]*9.81/1000  # 9.81 is the gravitational accelerator
#             availhydro_rc[t] = energy_rc
# =============================================================================
    # =============================================================================
        
    
    res_ops[name].update(energy=energy, availhydro=availhydro, #availhydro_rc=availhydro_rc,
                         s=s, r=r, Hd_m=Hd_m, Qspil=Qspil)
                         #s_rc=s_rc, r_rc=r_rc, Hd_m_rc=Hd_m_rc, Qspil_rc=Qspil_rc)
    
    return res_ops


# del energy,energy_rc,availhydro,availhydro_rc,s,s_rc,r,r_rc,Hd_m,Hd_m_rc,Qspil,Qspil_rc,
# del Hmin,Hmax,V,Vdesign,Qmax,Q_in,hydhead