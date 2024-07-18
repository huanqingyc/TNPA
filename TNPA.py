from TNPA_Epidemics import *
import argparse

'''
对代码进行了初步调整，将各算法分离出来，增强可读性
'''

def Interaction_module(): # unfinished
    print('Graph:')
    gtype = input("graph_type: ")
    if gtype in ['club','macrotree']:
        G,gname = graph(gtype)
    else:
        gn = eval(input("graph_size: "))
        if gtype == 'qubic':
            num = gn//2-1
            G,gname = graph(gtype,[num])
        else:
            seed = eval(input("seed: "))
            if gtype == 'tree':
                G,gname = graph(gtype,[gn,seed])
            else:
                c = eval(input("average_degree: "))
                if gtype == 'smallworld':
                    p = eval(input("reconnect: "))
                    G,gname = graph(gtype,[gn,c,p,seed])
                else:
                    if gtype == 'erg':
                        G,gname = graph(gtype,[gn,c/gn,seed])
                    else:
                        G,gname = graph(gtype,[gn,c,seed])
    clique = input('Add clique? y or n ')
    if clique == 'y':
        n = eval(input('number of clique: '))
        stepped = input("stepped size or constant size? s or c ")
        if stepped == 'c':
            stepped = False
        elif stepped == 's':
            stepped = True
        size = eval(input(int(stepped) * 'max' + 'size of clique(>=3): '))
        G = add_clique(G,n,size,seed,stepped)
        gname = str(n) + '*' + int(stepped) * '(3,' + str(size) + stepped*')' + 'cliqued_' + gname
    
    propose = input("target(r for region,e for evolution and epidemic): ")

    if propose == 'r': # 测试Region_generator
        method = input("method(s for one group of set, m for multi): ")
        if method == 's':
            R = eval(input("R= "))
            N = eval(input("N= "))
            print_region(G,R,N)
        elif method == 'm':
            R = eval(input("R_max= "))
            for r in range(3,R+1):
                print_region(G,r,r)
    elif propose == 'e':
        print("Evolution parameter")
        print("About time:")
        default = input("use default setting (200/0.1)?:y or n ")
        if default == 'y':
            t_max = 200
            tau = 0.1
        else:
            t_max = eval(input('t_max:'))
            tau = eval(input('tau:'))
        print("About Epidemic")
        default = input("use default setting (SIR:0.1,0.05)?:y or n ")
        if default == 'y':
            etype = 'SIR'
            l = 0.1
            r = 0.05
        else:
            etype = input("etype:")
            l = eval(input('lambda:'))
            r = eval(input('rho:'))
        epar = [l,r]
        init = None
        inits = []
        while init != '':
            init = input('init(press enter to continue): ')
            if init != '':
                inits.append(eval(init))
        
        print("About method")
        target = input("target(e for data and error up to MCMC, d for data only): ")
        if target == 'e':
            print('set MCMC:')
            default = input("use default setting (1000000,40)?:y or n ")
            if default == 'y':
                repeats = 1000000
                mp_num = 40
            elif default == 'n':
                repeats = eval(input('MCMC_repeats:'))
                mp_num = eval(input('MCMC_multiprocess_number:'))
        print("input methods("+(target == 'd')*"MCMC,"+"DMP,PA,TNPA), press enter to continue")
        method = None
        methods = []
        while True:
            method = input("method("+(target == 'd')*"MCMC,"+"DMP,PA,TNPA):")
            if method in methods:
                print("Repeated inputs. Try another one?")
            elif target == 'e' and method == 'MCMC':
                print("Obviously it's 0. What can I say? Try another one?")
            elif method in ['MCMC','DMP','PA','TNPA']:
                if method == 'TNPA':
                    RNs = []
                    while True:
                        print("print all R&N pairs you want:(press enter at any step to finish)")
                        r = input('R:')
                        if r != '':
                            n = input('N:')
                            if n != '':
                                RNs.append([eval(r),eval(n)])
                            else:
                                break
                        else:
                            break
                if method == 'MCMC':
                    repeats = input('repeats:')
                    mp_num = input('multiprocess_number:')
                methods.append(method)
            elif method == '':
                break
            else:
                print("Sorry, this method is unavailable now. Try another one?")

        if target == 'e' or 'MCMC' in methods:
            s_mc = Epidemic(G,gname,etype,epar,tau)
            s_mc.sys_init(inits)
            repeats = input('')
            s_mc.MCMC_init(repeats= repeats,mp_num=mp_num)
            s_mc.update_to(t_max)
            s_mc.save_data(precision)
        methods.remove('MCMC')
        simulations = []
        for method in methods:
            s = Epidemic(G,gname,etype,epar,tau)
            s.sys_init(inits)
            if method == 'TNPA':
                while len(RNs)>0:
                    rn = RNs.pop()
                    [r,n] = rn
                    s.TNPA_init(r,n)
                    if len(RNs)>0:
                        s.update_to(t_max)
                        s.save_data(precision)
                        simulations.append(s)
                        s = Epidemic(G,gname,etype,epar,tau)
                        s.sys_init(inits)
            else:
                if method == 'DMP':
                    s.DMP_init()
                elif method == 'PA':
                    s.PA_init()
            s.update_to(t_max)
            s.save_data(precision)

        if target == 'e':
            for s in simulations:
                s.save_error(s_mc.marginal_all,precision)

                
    else:
        print('Wrong command')

if __name__ == '__main__':

    eps = 1e-7
    precision = 5
    np.set_printoptions(precision=3,floatmode='maxprec',suppress=True)

    if False:
        Interaction_module()

    # n = 0
    # G,gname= graph('club')
    G,gname= graph('USA')
    # t_max = 50

    # G,gname= graph('rrg',[200,3,1])
    # G,gname= graph('rrg',[1000,3,1])
    # n = 20
    # k = 3
    # seed = 1
    # t_max = 200

    # G,gname = graph('triangle_tree',[6,1,50])
    # G,gname = graph('squares',[20])
    # G,gname = graph('lattice',[10])
    
    # G,gname = graph('random_tree',[150,35])
    t_max = 200

    # n = 12
    # k = 3
    # seed_clique = 48

    # G = add_clique(G,n,k,seed_clique)
    # gname = str(n) + '*' + str(k) + 'cliqued_seed=' + str(seed_clique) + gname

    tau = 0.1

    # t_max = 2
    # tau = 1.
    
    # rho = 0.1
    rho = 0.05
    lamda = 0.1

    etype = 'SIR'
    epar = [lamda,rho]

    ini = [0]
    
    # elist = [(0,1),(0,2),(0,3),(4,1),(4,2),(5,1),(5,3),(6,2),(6,3)]
    # elist = [(0,1),(0,2),(1,2),(0,3),(4,0),(3,4)]
    # G,gname = graph('elist',['test',elist])

    # elist = [(0,1),(0,2),(1,3),(2,3),(4,2),(5,3),(5,4),(5,6),(6,0)]
    # G,gname = graph('elist',['test2',elist])

    # 模型参数
    if True:
        s = MCMC_mp(G,gname,etype,epar,tau,ini)
        s.evolution(t_max,repeats= 1000000,mp_num=20)
        s.save_data(precision)

        s = PA(G,gname,etype,epar,tau,ini)
        s.evolution(t_max)
        s.save_data(precision)

        s = DMP(G,gname,etype,epar,tau,ini)
        s.evolution(t_max)
        s.save_data(precision)

        region_dict = get_Regions_diction(G,3,7)
        label = '_R3N7'
        s = TNPA(G,gname,etype,epar,tau,ini,region_dict[3],label)
        s.evolution(t_max)
        s.save_data(precision)

        # R = 4
        # N = 4
        # region_dict = get_Regions_diction(G,R,N)
        # for r in range(3,R+1):
        #     print(r)
        #     label = '_R' + str(r)  + 'N' + str(N)
        #     s = TNPA(G,gname,etype,epar,tau,ini,region_dict[r],label)
        #     s.evolution(t_max)
        #     s.save_data(precision)

        for N in range(4,12,2):
            R = [4]
            region_dict = get_Regions_diction(G,max(R),N)
            for r in R:
                # print(r)
                label = '_R' + str(r)  + 'N' + str(N)
                s = TNPA(G,gname,etype,epar,tau,ini,region_dict[r],label)
                s.evolution(t_max)
                s.save_data(precision)
        
    else:
        print_region_diction(G,4,4,False,True)

# nohup python -u TNPA.py >/dev/null 2>&1 &