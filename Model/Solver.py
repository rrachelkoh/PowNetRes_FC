 def Solver(instance,day):
    from pyomo.opt import SolverFactory, SolverStatus
    from pyomo.core import Var
    from pyomo.core import Param
    from operator import itemgetter
    
    opt = SolverFactory("gurobi") ##SolverFactory("cplex")
    opt.options["threads"] = 1
    opt.options["TimeLimit"] = 30 ##in seconds

    result = opt.solve(instance) ##,tee=True to check number of variables
    if result.solver.status == SolverStatus.aborted: #max time limit reached #v1.3
        result.solver.status = SolverStatus.warning #change status so that results can be loaded
        
        
    instance.solutions.load_from(result) 
     
    #  #The following section is for storing and sorting results
    for v in instance.component_objects(Var, active=True):
        varobject = getattr(instance, str(v))
        a=str(v)
        if a=='hydro':   
             hydro_ = []
             for index in varobject:
                 if int(index[1]>0 and index[1]<25):
                    if index[0] in instance.h_nodes:
                        hydro_.append((index[0],index[1]+((day-1)*24),varobject[index].value)) 
    
        # if a=='solar':
        #       solar_ = []          
        #       for index in varobject:
        #           if int(index[1]>0 and index[1]<25):
        #             if index[0] in instance.s_nodes:
        #                 solar_.append((index[0],index[1]+((day-1)*24),varobject[index].value)) 
    
        #         if a=='wind':      
        #              for index in varobject:
        #                  if int(index[1]>0 and index[1]<25):
        #                     if index[0] in instance.w_nodes:
        #                         wind.append((index[0],index[1]+((day-1)*24),varobject[index].value))   
                            
    
        if a=='hydro_import':     
              hydro_import_ = []
              for index in varobject:
                  if int(index[1]>0 and index[1]<25):
                    hydro_import_.append((index[0],index[1]+((day-1)*24),varobject[index].value)) 
    
        if a=='vlt_angle':
             vlt_angle_ = []
             for index in varobject:
                 if int(index[1]>0 and index[1]<25):
                    if index[0] in instance.nodes:
                        vlt_angle_.append((index[0],index[1]+((day-1)*24),varobject[index].value))
    
        if a=='mwh':  
            mwh_ = []
            ini_mwh_ = {} #v1.3
            for index in varobject:
                if int(index[1]>0 and index[1]<25): 
                    mwh_.append((index[0],index[1]+((day-1)*24),varobject[index].value))  
                if int(index[1])==24:
                    ini_mwh_[index[0]] = varobject[index].value                            
    
                        
        if a=='on':   
            on_ = []
            ini_on_ = {}  
            for index in varobject:
                if int(index[1]>0 and index[1]<25):
                    # on.append((index[0],index[1]+((day-1)*24),varobject[index].value))
                    on_.append((index[0],index[1]+((day-1)*24),varobject[index].value))
                if int(index[1])==24:
                    ini_on_[index[0]] = varobject[index].value    
    
        if a=='switch':  
            switch_ = []
            for index in varobject:
                if int(index[1]>0 and index[1]<25):
                    switch_.append((index[0],index[1]+((day-1)*24),varobject[index].value))
    
        if a=='srsv':   
            srsv_ = []
            for index in varobject:
                if int(index[1]>0 and index[1]<25):
                    srsv_.append((index[0],index[1]+((day-1)*24),varobject[index].value))
                        
        if a=='nrsv':   
            nrsv_ = []
            for index in varobject:
                if int(index[1]>0 and index[1]<25):
                    nrsv_.append((index[0],index[1]+((day-1)*24),varobject[index].value))           
    
    # Update initialization values for "on" 
    for z in instance.Generators:
        instance.ini_on[z] = round(ini_on_[z])
        instance.ini_mwh[z] = max(ini_mwh_[z],0) #v1.3
        
    return hydro_,hydro_import_,vlt_angle_,switch_,on_,mwh_,nrsv_,srsv_  #,solar_

