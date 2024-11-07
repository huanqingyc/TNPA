import numpy as np
import networkx as nx
from TNPA import *
from scipy.stats import norm

def plot_temporal_CI(name1,name2,labels_TNPA,t_max = None):
    # mDMP = read_data(name1 + 'mDMP' + name2)[1][m]
    DMP = 1.-read_data(name1 + 'DMP' + name2)[1][0]
    MCMC = 1.-read_data(name1 + 'MCMC' + name2)[1][0]
    PA = 1.-read_data(name1 + 'PA' + name2)[1][0]

    TNPA = []
    for label in labels_TNPA:
        TNPA.append(1.-read_data(name1 + 'TNPA_' + label + name2)[1][0])
    
    if t_max == None:
        t_max = list(DMP.shape)[0]-1 # 包含零时刻
    # print(t_max)
    plt.figure(1)
    plt.scatter(range(0,t_max,10), MCMC[:t_max:10], label = 'MC', marker = 'o', edgecolors='black' ,c = 'None')

    # plt.plot(range(t_max+1),  mDMP[:t_max+1], label = 'mDMP')
    plt.plot(range(t_max+1),  PA[:t_max+1] , label = 'PA')
    plt.plot(range(t_max+1),  DMP[:t_max+1], label = 'DMP',linestyle='--',dashes = (8,8))
    if len(TNPA) > 1:
        for i in range(len(TNPA)):
            data = TNPA[i]
            label = labels_TNPA[i]
            plt.plot(range(t_max+1), data[:t_max+1], label = 'TNPA_'+label)
    else:
        data = TNPA[0]
        plt.plot(range(t_max+1), data[:t_max+1],'r', label = 'TNPA')

    plt.ylabel('Cumulative Infections')
    plt.xlabel('Time t')

    y_max = np.ceil(100*np.max(PA))/100.
    plt.ylim(0,y_max)
    # plt.ylim(0.1,1)
    plt.xlim(0,t_max)
    plt.legend()

def plot_error_CI(name1,name2,labels_TNPA,t_max = None):
    MCMC = read_data(name1 + 'MCMC' + name2)[0][:,0,:]
    # mDMP = np.average(abs(read_data(name1 + 'mDMP' + name2)[0][:,m,:]-MCMC),axis=0)
    DMP = np.average(abs(read_data(name1 + 'DMP' + name2)[0][:,0,:]-MCMC),axis=0)
    PA = np.average(abs(read_data(name1 + 'PA' + name2)[0][:,0,:]-MCMC),axis=0)
    
    TNPA = []
    for label in labels_TNPA:
        TNPA.append(np.average(abs(read_data(name1 + 'TNPA_' + label + name2)[0][:,0,:]-MCMC),axis=0))

    if t_max == None:
        t_max = list(DMP.shape)[0]-1 # 包含零时刻
        
    plt.plot(range(t_max+1),  PA[:t_max+1] , label = 'PA')
    plt.plot(range(t_max+1),  DMP[:t_max+1], label = 'DMP',linestyle='--',dashes = (8,8))
    if len(TNPA) > 1:
        for i in range(len(TNPA)):
            data = TNPA[i]
            label = labels_TNPA[i]
            plt.plot(range(t_max+1), data[:t_max+1], label = 'TNPA_'+label)
    else:
        data = TNPA[0]
        plt.plot(range(t_max+1), data[:t_max+1],'r', label = 'TNPA')


    plt.ylabel('$L_1$ Error')

    plt.xlabel('Time t')
    y_max = np.ceil(100*np.max(np.max(PA)))/100.
    plt.ylim(0,y_max)
    # plt.ylim(0.1,1)
    plt.xlim(0,t_max)
    plt.legend()

def plot_snap_CI(name1,name2,label_TNPA,t,bins = 10,other_label = False):
    # mDMP = 1-read_data(name1 + 'mDMP' + name2)[0][:,0,t]
    DMP = 1-read_data(name1 + 'DMP' + name2)[0][:,0,t]
    MCMC = 1-read_data(name1 + 'MCMC' + name2)[0][:,0,t]
    PA = 1-read_data(name1 + 'PA' + name2)[0][:,0,t]

    TNPA = 1-read_data(name1 + 'TNPA_' + label_TNPA + name2)[0][:,0,t]
    if other_label:
        other = 1-read_data(name1 + 'TNPA_' + other_label + name2)[0][:,0,t]

    m_min = (10*np.min(MCMC)//1)*0.1
    m_max = np.ceil(10*np.max(DMP))/10.
    bins = np.linspace(m_min,m_max,bins+1)

    # mDMP = norm.pdf(bins, np.mean(mDMP), np.std(mDMP))
    DMP = norm.pdf(bins, np.mean(DMP), np.std(DMP))
    MCMC = norm.pdf(bins, np.mean(MCMC), np.std(MCMC))
    PA = norm.pdf(bins, np.mean(PA), np.std(PA))
    TNPA = norm.pdf(bins, np.mean(TNPA), np.std(TNPA))
    if other_label:
        other = norm.pdf(bins, np.mean(other), np.std(other))

    plt.figure(1,figsize = (8,5))
    
    # print(bins,counts)
    # mid = (bins[:-1]+bins[1:])/2.

    plt.plot(bins,  MCMC,'black', label = 'MC')
    # plt.plot(bins,  mDMP,'g', label = 'mDMP')
    plt.plot(bins,  PA , label = 'PA')
    plt.plot(bins,  DMP, label = 'DMP',linestyle='--',dashes = (8,8))
    if other_label:
        plt.plot(bins, TNPA ,'r', label = 'TNPA_'+label_TNPA)
        plt.plot(bins, other ,'g', label = 'TNPA_'+other_label)
    else:
        plt.plot(bins, TNPA ,'r', label = 'TNPA')
        

    # plt.xlim(0,I_max)

    # 添加标题和坐标轴标签
    # plt.title('Comparison of Data Distributions')
    plt.xlim(m_min,m_max)
    # y_max = np.ceil(10*np.max(MCMC))/10.
    plt.ylim(0,)
    plt.xlabel('Cumulative Infection')
    plt.ylabel('Density Function')

    plt.legend()
    
def print_late_time_CI(name1,name2,labels_TNPA,precision=5):
    print("Cumulative Infection:")
    MCMC = 1-read_data(name1 + 'MCMC' + name2)[0][:,0,-1]
    print('MCMC:',np.round(np.average(MCMC),precision))
    # mDMP = read_data(name1 + 'mDMP' + name2)[0][:,0,-1]
    # print('mDMP:',np.round(np.average(mDMP),precision))
    # print('error_mDMP:',np.round(np.average(abs(mDMP-MCMC)),precision))
    DMP = 1-read_data(name1 + 'DMP' + name2)[0][:,0,-1]
    print('DMP:',np.round(np.average(DMP),precision))
    print('error_DMP:',np.round(np.average(abs(DMP-MCMC)),precision))
    PA = 1-read_data(name1 + 'PA' + name2)[0][:,0,-1]
    print('PA:',np.round(np.average(PA),precision))
    print('error_PA:',np.round(np.average(abs(PA-MCMC)),precision))
    for label in labels_TNPA:
        TNPA = 1-read_data(name1 + 'TNPA_' + label + name2)[0][:,0,-1]
        print('TNPA'+label+':',np.round(np.average(TNPA),precision))
        print('error_TNPA'+label+':',np.round(np.average(abs(TNPA-MCMC)),precision))

# I
name1 = '[0.1, 0.05]_'

# # RRG
# name2 = '_T=200_tau=0.1'
# gname = 'rrg_c=3_seed=1_n=200_'
# figname = 'rrg'
# # gname = '10*3cliqued_rrg_c=3_seed=1_n=200_'
# # figname = 'cliqued_rrg'
# labels_TNPA = ['R'+str(r)+'N8' for r in range(3,9)]
# label = 'R6N8'

# # gname = 'rrg_c=3_seed=1_n=1000_'
# # gname = '20*3cliqued_rrg_c=3_seed=1_n=1000_'
# labels_TNPA = ['R'+str(r)+'N9' for r in range(3,9)]
# label = 'R6N9'
# figname = 'cliqued_rrg'

# # RRG
# name2 = '_T=300_tau=0.1'
# gname = 'rrg_c=3_seed=1_n=500_'
# figname = 'rrg'
# # gname = '30*3cliqued_rrg_c=3_seed=1_n=500_'
# # figname = 'cliqued_rrg'
# labels_TNPA = ['R'+str(r)+'N7' for r in range(3,8,2)]
# label = 'R7N7'

# # Cliqued_tree
# name2 = '_T=500_tau=0.1'
# gname = '12*3cliqued_random_tree_seed=35_n=150_'
# labels_TNPA = ['R'+str(r)+'N8' for r in range(3,9) if r != 4]
# label = 'R8N8'
# other_label = 'R3N8'
# figname = 'cliqued_rand_tree'

# # Cliqued_tree
# name2 = '_T=500_tau=0.1'
# gname = '15*3cliqued_random_tree_seed=29_n=150_'
# labels_TNPA = ['R'+str(r)+'N9' for r in range(3,8)]
# label = 'R9N9'
# other_label = 'R3N9'
# figname = 'cliqued_rand_tree'

# Macrotree
# name2 = '_T=100_tau=0.1'
# gname = 'macrotree_n=24_'
# labels_TNPA = ['R3N7','R5N7']

# # Club
# name2 = '_T=200_tau=0.1'
# gname = 'karate_club_n=34_'
# # labels_TNPA = ['R'+str(r)+'N'+str(n) for r in [4] for n in range(4,7)]
# labels_TNPA = ['R'+str(r)+'N'+str(n) for r in [3,4] for n in range(4,7)] # 末态好，但过程差
# label = 'R4N6'
# figname = 'club'

# USA
name2 = '_T=200_tau=0.1'
gname = 'contiguous_usa_n=49_'
labels_TNPA = ['R'+str(r)+'N'+str(n) for r in [3] for n in range(3,7)]
label = 'R3N7'
figname = 'usa'

# # Lattice
# name2 = '_T=300_tau=0.1'
# gname = 'lattice_n=100_'
# labels_TNPA = ['R4N'+str(n) for n in range(4,12,2)]
# figname = 'lattice'

# # Squares_chain
# name2 = '_T=400_tau=0.1'
# gname = 'squares_n=42_'
# labels_TNPA = ['R4N4']
# label = 'R4N4'
# figname = 'square_chain'

# # Cross_square_chain
# name2 = '_T=400_tau=0.1'
# gname = 'cross_square_n=31_'
# labels_TNPA = ['R4N4']
# label = 'R4N4'
# figname = 'cross_square'

# # Triangle_tree
# name2 = '_T=400_tau=0.1'
# gname = 'triangle_tree_n=290_'
# labels_TNPA = ['R3N3']
# label = 'R3N3'
# figname = 'triangle_tree'

name1 += gname

# plot_temporal_CI(name1,name2,labels_TNPA,200)
# # plt.savefig(figname+'_CI.pdf')

plot_error_CI(name1,name2,labels_TNPA,100)
# plt.savefig(figname+'_CI_error.pdf')

# t = 150
# plot_snap_CI(name1,name2,label,t,bins = 30)
# plt.savefig(figname+'_snap'+'_t='+str(t)+'.pdf')

# print(figname)
# print_late_time_CI(name1,name2,labels_TNPA,5)
