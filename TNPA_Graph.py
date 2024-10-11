from copy import deepcopy
import networkx as nx
import numpy as np
import time

'''
现行版本，intersection就当他不存在吧
调整subg的分割方案，以保证R增大结果一定变好的一致性。
具体实现为：先挑出R的Region后，在这个Region中搜索R-1的Region，先生成内部的分割方案，再在其基础上向外拓展，生成剩余部分的subg分割方案。
'''

def print_region(G,R,N):
    partition = get_partition(G,R,N) 
    nodes = set()
    for region in partition[R]:
        nodes.update(list(region))
        print(list(region.graph),list(region.graph.edges()))
    print(len(nodes))

def print_region_diction(G,R,N,edges:bool,macro_G:bool = False):
    regions_dict = get_partition(G,R,N) 
    # print(len(regions))
    for r in range(3,R+1):
        all_nodes = set()
        print('R='+str(r))
        for region in regions_dict[r]:
            all_nodes.update(list(region.graph))
            print(list(region.graph)) #,list(region.graph.edges())
            if edges:
                print(list(region.graph.edges()))
            if region.subg:
                print('subregions:')
                for part in region.subregions:
                    print(list(part))
                if macro_G:
                    print(list(region.macro_g.edges()))
            print()
        print(len(all_nodes))

def cut(G,nodes = False):
    G_loop = remove_empty_nodes(G,nodes)
    if nodes != False:
        degree = dict(nx.degree(G_loop,nodes)) # 很好的鲁棒函数，不存在的点直接略过,这样即使前一步已经把无连边的点去掉了，但是还是能给出想要的结果 
    else:
        nodes = list(G_loop.nodes())
        degree = dict(nx.degree(G_loop))

    if len(degree)>0: # 避免因为前面已经把所有点都删掉了，或者给定了参考范围nodes的点都已经被删了
        min_degree_nodes = [node for node, deg in degree.items() if deg == min(degree.values())]
        if degree[min_degree_nodes[0]] == 1:
            for i in min_degree_nodes:
                j = list(G_loop[i])[0]
                G_loop.remove_node(i)
                while len(G_loop[j]) == 1:
                    k = list(G_loop[j])[0]
                    G_loop.remove_node(j)
                    j = k
                if len(G_loop[j]) == 0:
                    G_loop.remove_node(j)
                    break

    return G_loop

def remove_empty_nodes(G,nodes = False):
    """
    
    """
    if nodes != False:
        degree = dict(nx.degree(G,nodes))
    else:
        degree = dict(nx.degree(G))
    if len(degree)>0:
        G_n = deepcopy(G)
        min_degree_nodes = [node for node, deg in degree.items() if deg == min(degree.values())]

        if degree[min_degree_nodes[0]] == 0:
            for i in min_degree_nodes:
                G_n.remove_node(i)
        return G_n
    else:
        return G

def Region_generator(G,e,R:int):
    """Generate a subgraph Region(R) from the given edge e on the corresponding graph G, 
    which satisfies that the shortest loop cross the Region and G\Region is longer than R (cross means share more than 1 edge here).
    
    Parameters
    ----------
    G : nx.Graph
        The complete graph.
    e : tuple of int
        The beginning edge, it should be an element of list(G.edges()).
    R : int

    Returns
    -------
    Region : nx.Graph
        The subgraph.
    """
    Region = nx.Graph()
    G_environment = deepcopy(G)
    G_environment.remove_edge(e[0],e[1])
    Region.add_edge(e[0],e[1])
    boundary = list(e)
    edge_new = True
    l = 2

    while edge_new:
        edge_new = False

        new_path = []

        for i in range(1,l):
            for j in range(i):
                n1 = boundary[i]
                n2 = boundary[j]
                if nx.has_path(G_environment,n1,n2):
                    shortest_paths_out = list(nx.all_shortest_paths(G_environment, n1, n2))
                    l_shortest_paths_in = nx.shortest_path_length(Region, n1, n2)
                    if len(shortest_paths_out[0]) + l_shortest_paths_in <= R+1:
                        edge_new = True
                        new_path += shortest_paths_out
                                    
        outnode = set()
        for path in new_path:
            outnode.update(path)
            for node_id in range(len(path)-1):
                if not Region.has_edge(path[node_id],path[node_id+1]):
                    Region.add_edge(path[node_id],path[node_id+1])
                    G_environment.remove_edge(path[node_id],path[node_id+1])

        if edge_new:
            boundary = []
            for node in outnode: # 由于我们考虑的是圈的长度，因此如果一个点在一轮中没有进入outnode，也就意味着这个点不属于任何一个和上一轮区域相交的小环，那么新添加了点之后也不会属于了（因为该点到新的点必然经过上一轮判定过的点，所属的圈也必然包含途径的点）
                if len(Region[node])<len(G[node]):
                    boundary.append(node)
        l = len(boundary)
        
        if l < 2:
            break
        
    return Region

def get_partition(G,R:int,N:int):
    G = cut(G)
    G_remain = deepcopy(G)
    Regions_dict = dict() # 字典的key为3到R,对应相应R取值的Region变量
    for r in range(3,R+1):
        Regions_dict[r] = []
    edges = list(G_remain.edges())
    for e in edges:
        # print(e)
        if G_remain.has_edge(e[0],e[1]):
            g = Region_generator(G_remain,e,R)
            G_remain.remove_edges_from(list(g.edges())) # 不管怎样，g的所有边都可以从G_remain上移除了
            if len(g) > 2:# 在一条边的基础上找到了其他区域
                local_partition = get_local_partition(g,R,N)
                Regions_dict = (Regions_dict,local_partition,R+1)
            if len(G_remain)>0:
                G_remain = cut(G_remain,list(g))
    return Regions_dict

def subregion_generator(g,g_full,e,R:int,N:int):
    # g是当前还没被划分的，g_full是完整的Region
    g_temp = deepcopy(g_full) #是临时用来找subg以外的边的
    g_left = deepcopy(g) # 
    subg = nx.Graph()
    # 关系如下：subg+g_left=g，subg+g_temp = g_full
    subg.add_edge(e[0],e[1])
    g_left.remove_edge(e[0],e[1])
    g_temp.remove_edge(e[0],e[1])
    edge_new = True
    set_g = set(g)
    # 先找到新的subregion
    
    while edge_new:
        boundary = [node for node in set(subg).intersection(set_g) if len(subg[node])<len(g[node]) or len(subg[node]) ==1]
        # for node in set(subg).intersection(set_g):
        #     if len(subg[node])<len(g[node]) or len(subg[node]) ==1: # 还有可用的连接,或是在边界处
        #         boundary.append(node)
        # if subg.has_edge(12,10) or subg.has_edge(42,36):
        edge_new = False
        full = False

        for r in range(3,R+1):
            l = len(boundary)
            if l < 2:
                break
            for i in range(1,l):
                for j in range(i):
                    n1 = boundary[i]
                    n2 = boundary[j]
                    if nx.has_path(g_temp,n1,n2):
                        paths_out = list(nx.all_shortest_paths(g_temp,n1,n2))
                        l_path_in = nx.shortest_path_length(subg,n1,n2)
                        if len(paths_out[0]) + l_path_in <= r+1:                           
                            flag = False
                            # 找一个有意义的路径
                            while len(paths_out)>0:
                                path_out = paths_out.pop()
                                if len(paths_out)>0 and len(set(path_out).intersection(set_g)) == 2:#就只有开头结尾两个点，没意义
                                    continue
                                else:
                                    flag = True
                                    break

                            if flag: # 找到的路径有意义的情况下才添加 
                                edge_new = True
                                # print(path_out,list(subg))
                                if len(path_out)+len(set(subg).intersection(set_g)) > N+2: # 添加新边会超范围，就跳过
                                    full = True
                                else:
                                    full = False
                                    for node_id in range(len(path_out)-1):
                                        [n1,n2] = path_out[node_id:node_id+2]
                                        subg.add_edge(n1,n2) # 在这里会被添加上虚边
                                        g_temp.remove_edge(n1,n2)
                                        if g_left.has_edge(n1,n2):
                                            g_left.remove_edge(n1,n2)
                                # print('getmp',list(g_left))
                                    break
                if edge_new: # 一次只能添加一条新的边
                    break
            if edge_new: # 一次只能添加一条新的边
                break 
        if full: # 当添加任一新路径都会超过上限或恰好抵达上限时break
            break
        elif len(set(subg).intersection(set_g))==N:
            # 查漏补缺，此时可能出现现有的点之间有连边，但是没必要再过一遍search path了
            set_subg = set(subg)
            for node in list(subg):
                if node in list(g_left):
                    neighs_inside = set(g_left[node]).intersection(set_subg)
                    for neigh in neighs_inside:
                        subg.add_edge(node,neigh)
                        g_left.remove_edge(node,neigh)
            break


    # 不再保留subg中有意义的与g的公共边
    nodes = set()
    for e in list(subg.edges()):
        [n1,n2] = e
        if not g.has_edge(n1,n2):
            subg.remove_edge(n1,n2)
            nodes.update([n1,n2])
    subg = remove_empty_nodes(subg,list(nodes))

    if len(subg)>0: # 其实是废话，一定会>0的
        g_left = remove_empty_nodes(g_left,list(subg))

    return [subg,g_left]

def partition_dict_update(dict_R,dict_r,R):
    for r in range(3,R):
        dict_R[r] = dict_R[r] + dict_r[r]
    return dict_R

def split_region(g,N):
    regions = []
    g_left = deepcopy(g)
    edges = list(g.edges())
    for e in edges:
        n1,n2 = e
        if g_left.has_edge(n1,n2):
            region = nx.Graph()
            region.add_edge(n1,n2)
            g_left.remove_edge(n1,n2)
            paths = list(nx.all_shortest_paths(g_left, n1, n2))
            l = len(paths[0])-2
            for path in paths:
                if len(region)+l<=N:
                    for node_id in range(len(path)-1):
                        region.add_edge(path[node_id],path[node_id+1])
                        g_left.remove_edge(path[node_id],path[node_id+1])
                else:
                    break
            regions.append(region)
            g_left = cut(g_left,list(region))

    return regions

def search_connection(G,N,b1,n1,b2,n2):
    connections = []
    for a in b1:
        for b in b2:
            if a==b:
                connections.append([])
            elif nx.has_path(G,a,b):
                connections += list(nx.all_simple_paths(G, n1, n2, N-n1-n2+2))
    return connections

def add_path():
    None
def remove_path():
    None
    # claim:一个圈只要被断掉，那么不管分成几份都等同于用PA即'目->口+||+口 =口+|+|+口  or 口+凵+凵'其中后者反而更麻烦一些
def get_local_partition(g,R,N):
    local_dict = dict() # 字典的key为3到R,对应相应R取值的Region变量
    for r in range(3,R+1):
        local_dict[r] = []
    # 检索并删除R-1的region
    G_remain = deepcopy(g) # 用于搜索子图的剩余区域，找不到region的边会被删掉，去掉子图后没有用的摇摆边(必然属于更长的圈)会被剪掉
    g_left = deepcopy(g) # 只删除region的d真正的剩余区域
    if R>3:
        edges = list(G_remain.edges())
        for e in edges:
            if G_remain.has_edge(e[0],e[1]):
                g = Region_generator(G_remain,e,R-1)
                G_remain.remove_edges_from(list(g.edges))
                if len(g) > 2:
                    local_partition = get_local_partition(g,R-1,N)
                    local_dict = (local_dict,local_partition,R)
                    for region in local_partition[R-1]:
                        g_left.remove_edges_from(list(region.edges)) 
                        # 把小region里的边去掉,该操作会保留在小region中因为N的限制而本质上并没有被考虑的边
                if len(G_remain)>0:
                    G_remain = cut(G_remain,list(g))

    if len(list(g))<=N:
        local_dict[R].append(g)
    else:
        # 检索新添加区域和在R-1部分被舍去的边构成的region
        g_left = remove_empty_nodes(g_left)# R新增的区域以及R-1时没有被考虑的部分
        regions_in_g_left = []
        G_parts = list(nx.connected_components(g_left)) # 经常会是不止一个连通区域
        for G_remain in G_parts:
            if not nx.is_tree(G_remain): # 树图就略过吧
                G_remain = cut(G_remain) # 去掉无意义的摇摆边
                if len(G_remain)<=N: # 剩下的都是圈，足够小就直接吃掉
                    regions_in_g_left.append(G_remain)
                else:
                    regions_in_g_left += split_region(G_remain,N)

        for g in regions_in_g_left:
            g_left.remove_edges_from(list(g.edges()))
        g_left = remove_empty_nodes(g_left)# 最后剩下的连接边

        if len(g_left)>0:# 万一正好清空呢XP
            last_regions = []
            for region in (local_dict[R-1]+regions_in_g_left):
                if len(region)<N: # 还能继续添加点的regions
                    last_regions.append(region)
                else:
                    local_dict[R].append(region)
            m = len(last_regions)
            boundaries = [list(nx.intersection(region,g_left)) for region in last_regions] # 这些region的边界

            # 自连接的合并：
            for i in range(m):
                region = last_regions[i]
                n = len(region)
                b = boundaries[i]
                l = len(b)
                if l>1:# 在遗留的边中没有自圈
                    connections = []
                    for j in range(l):
                        for k in range(j+1,l):
                            connections += list(nx.all_simple_paths(G, b[j], b[k], N-n+2)) # 添加后必然超上限的直接去掉
                    connections = sorted(connections,key = len)
                    for path in connections:
                        g_test = deepcopy(region)
                        for node_id in range(len(path)-1):
                            g_test.add_edge(path[node_id],path[node_id+1])
                        if len(g_test)<=N:
                            region = g_test
                            for node_id in range(len(path)-1):
                                g_left.remove_edge(path[node_id],path[node_id+1])
                    if len(region) == N: # 因为是先从R-1的region开始延伸，应该是不可能出现恰好重合，能合并的情况即len(g1+g2)=len(g1)
                        local_dict[R].append(region)
                        del last_regions[i]
                        del boundaries[i]

            # 剩余的小region间的合并：     
            # 先搜索g_left构成的连接,并合并尽量邻居
            m = len(last_regions)
            region_graph = nx.Graph()
            while True: # 反复尝试合并邻居，直到遍历后也不能合并为止
                merged = False
                bad_region_id = [] # 用来存储过程中发现的和所有邻居直接合并都太大的region
                for i in range(m):
                    for j in range(i+1,m):
                        n1,n2 = len(last_regions[i]),len(last_regions[j])
                        if len(set(list(last_regions[i])).update(list(last_regions[j])))>N:# 直接拼起来都不行直接跳过
                            No_small_neighbor = True
                            continue
                        No_small_neighbor = False
                        connection = []
                        for a in boundaries[i]:
                            for b in boundaries[j]:
                                if a==b:
                                    connection.append([])
                                elif nx.has_path(G,a,b):
                                    connection += list(nx.all_simple_paths(G, n1, n2, N-n1-n2+2))
                        if len(connection)>1: # 有两条及以上的路径则试着合并
                            # 通过最短的两条路径合并
                            connection = sorted(connection,key = len)
                            g_test = nx.union(last_regions[i],last_regions[j])
                            path = connection[0]
                            for node_id in range(len(path)-1):
                                g_test.add_edge(path[node_id],path[node_id+1])
                            path = connection[1]
                            for node_id in range(len(path)-1):
                                g_test.add_edge(path[node_id],path[node_id+1])
                            if len(g_test)<=N: # 足够小
                                region = g_test
                                merged = True
                                for node_id in range(len(path)-1):
                                    g_left.remove_edge(path[node_id],path[node_id+1])
                                path = connection[0]   
                                for node_id in range(len(path)-1):
                                    g_left.remove_edge(path[node_id],path[node_id+1])
                                # 尝试合并其他的路径(相当于是新region的自环),哪怕恰好是N也要试一下，避免第三小长度为2（直连）,因此==N没有单独处理
                                connection = connection[2:]
                                for path in connection:
                                    g_test = deepcopy(region)
                                    for node_id in range(len(path)-1):
                                        g_test.add_edge(path[node_id],path[node_id+1])
                                    if len(g_test)<=N:
                                        region = g_test
                                        for node_id in range(len(path)-1):
                                            g_left.remove_edge(path[node_id],path[node_id+1])
                                    else: # 因为是从小到大排列的，所以短的不行长的更不行
                                        break
                        if merged:
                            break
                    if merged:# 能merge一定不No_small_neighbor
                        break
                    if No_small_neighbor:
                        bad_region_id.append(i)

                # 每次开始下一次循环前删掉没用的东西
                bad_region_id = sorted(bad_region_id,reverse=True)# 得从后往前删
                for k in bad_region_id: 
                    local_dict[R].append(last_regions[k])
                    del last_regions[k]
                    del boundaries[k]

                if merged: #合并相应变量
                    del last_regions[j]
                    del boundaries[j]
                    if len(region)==N:# 同上，不考虑极小概率的合并情况
                        local_dict[R].append(region)
                        del last_regions[i]
                        del boundaries[i]
                    else:
                        last_regions[i] = region
                        boundaries[i] = list(nx.intersection(region,g_left))
                else:# 遍历后也没有成功
                    break
            
            # 邻居之间尝试合并完成，搭建剩余的不太大的region构成的graph用来搜索可多区块合并的结构
            '''
            好麻烦啊...试试看常见的图到这一步复不复杂，实在不行就遍历
            '''
            last_regions = sorted(last_regions,key = len)# 从小到大排列
            m = len(last_regions)
            connections = [[False for _ in range(m)] for __ in range(m)]
            region_graph = nx.Graph()
            for i in range(m):
                for j in range(i+1,m):
                    n1,n2 = len(last_regions[i]),len(last_regions[j])
                    if len(set(list(last_regions[i])).update(list(last_regions[j])))>N:# 直接拼起来都不行直接跳过
                        continue
                    connection = []
                    for a in boundaries[i]:
                        for b in boundaries[j]:
                            if a==b:
                                connection.append([])
                            elif nx.has_path(G,a,b):
                                connection += list(nx.all_simple_paths(G, n1, n2, N-n1-n2+2))
                    if len(connection) > 0:
                        connections[i,j] = connection
                        region_graph.add_edge(i,j)


            # 新的拓展方案：
            # 构建一个新的图：每一个R-1的region是一个节点，一条或多条边存在在节点和节点之间以及节点自身，节点有权重n，边有长度l
            # (n=N的节点可以省略，n+n_neighbour>N or all{n+n_neighbour+l_min+l_min2}>N的同理,后面两个放在循环过程中吸收自环之后可能更好)
            # 因为每一个圈的长度都是R，造成的误差是同一数量级的，因此假定不分先后
            # 重复以下循环：
            # 首先试着吸收节点到自身的边，若n+l<=N则吸收之，吸收后n+l=N则将其去除。
            # 接着吸收两点间的边：若n1+n2+l_min+l_min2<=N则吸收,吸收后n=N则将其去除。
        else:# R-1的partition已经覆盖了这个Region
            local_dict[R] = local_dict[R-1] + regions_in_g_left
    return local_dict

# about graph

def graph(G_type,parameter = []):
    '''
    G_type
    '''
    gname = G_type
    
    if G_type == 'rrg':
        [n,c,seed] = parameter
        g = nx.random_regular_graph(c,n,seed=seed)
        gname += '_c=' + str(c) + '_seed=' + str(seed)
    elif G_type == 'erg': # 也即泊松度分布
        [n,c,seed] = parameter
        g = nx.erdos_renyi_graph(n,c/n,seed=seed)
        gname += '_c=' + str(c) + '_seed=' + str(seed)
    elif G_type == 'tree':
        [n,seed] = parameter
        g = nx.random_tree(n,seed = seed)
    elif G_type == 'powerlaw': # 无标度网络
        [n,m,seed] = parameter
        g = nx.barabasi_albert_graph(n,m,seed = seed)
        gname += '_m=' + str(m) + '_seed=' + str(seed)
    elif G_type == 'smallworld': # 小世界网络
        [n,c,p,seed] = parameter
        g = nx.watts_strogatz_graph(n,c,p,seed = seed)
        gname += '_c=' + str(c) + '_p=' + str(p) + '_seed=' + str(seed)
    elif G_type == 'macrotree':
        g = macrotree()
        n = len(list(g))
    elif G_type == 'squares':
        num = parameter[0]
        g = qubic_graph(num)
        n = len(list(g))
    elif G_type == 'elist':
        [gname,elist] = parameter
        g = nx.empty_graph()
        g.add_edges_from(elist)
        n = len(list(g))
    elif G_type == 'club':
        gname = 'karate_club'
        g = nx.karate_club_graph()
        n = len(list(g))
    elif G_type == 'bcspwr':
        [num] = parameter
        gname += str(num)
        g = nx.empty_graph()
        with open('bcspwr0'+str(num)+'.txt','r') as f:
            while True:
                e_str = f.readline()[:-1]
                if e_str :
                    n1,n2= e_str.split()
                    n1 = int(n1)-1 # 存储的是矩阵，从1开始，有error
                    n2 = int(n2)-1
                    if n1 != n2:
                        g.add_edge(n1,n2)
                else:
                    break
        n = len(g)
    elif G_type == 'random_tree':
        [n,seed] = parameter
        g = nx.random_tree(n,seed=seed)
        gname += '_seed=' + str(seed)
    elif G_type == 'USA':
        gname = 'contiguous_usa'
        g = nx.empty_graph()
        with open('contiguous_usa.txt','r') as f:
            while True:
                e_str = f.readline()[:-1]
                if e_str :
                    n1,n2= e_str.split(sep=',')
                    n1 = int(n1)-1 # 存储的是矩阵，从1开始，有error
                    n2 = int(n2)-1
                    if n1 != n2:
                        g.add_edge(n1,n2)
                else:
                    break
        n = len(g)

    elif G_type == 'triangle_tree':
        [layer,seed,tri_num] = parameter
        np.random.seed(seed)
        g = nx.empty_graph()
        g.add_edges_from([(0,1),(0,2),(0,3)])
        surface = [1,2,3]
        for _ in range(1,layer):
            present = len(g)
            for node in surface:
                l = len(g)    
                g.add_edge(node,l)
                g.add_edge(node,l+1)
            surface = list(range(present,len(g)))

        replace_list = list(np.random.choice(present, tri_num, replace = False))
        present = len(g)
        for node in replace_list:
            neighs = list(g[node])
            g.remove_edge(node,neighs[1])
            g.remove_edge(node,neighs[2])
            g.add_edges_from([(node,present),(node,present+1),(present,present+1),(present,neighs[1]),(present+1,neighs[2])])
            present += 2
        n = len(g)

    elif G_type == 'lattice':
        [n] = parameter
        g = lattice(n)
        n = len(list(g))

    gname += '_n=' + str(n)

    return g,gname

def qubic_graph(num):
    elist = [(0,1)]
    for i in range(num):
        n = i*2
        elist+=[(n,n+2),(n+1,n+3),(n+2,n+3)]
    G = nx.empty_graph()
    G.add_edges_from(elist)

    return G

def macrotree():
    elist=[(0,1),(0,2),(1,2), # 基础三角形
        (0,3),(1,4),(2,5), # 三条外延的边
        (3,6),(3,7),(4,8),(4,9),(5,10),(5,11),# 再向外延伸一层
        (6,12),(6,13),(12,13),(7,14),(7,15),(14,15),(8,16),(8,17),(16,17), # 六个三角形 
        (9,18),(9,19),(18,19),(10,20),(10,21),(20,21),(11,22),(11,23),(22,23),
        (13,14),(17,18),(21,22)]
    G = nx.empty_graph(0)
    G.add_edges_from(elist)
    return G

def lattice(n):
    g = nx.empty_graph()
    g.add_edges_from([(cols+n*row,cols+n*row+1) for cols in range(n-1) for row in range(n)])
    g.add_edges_from([(cols+n*row,cols+n*row+n) for cols in range(n) for row in range(n-1)])
    return g

def add_clique(G,n_clique,size,seed,stepped=False):
    """在一张图上生成局部的clique
    Parameters
    ----------
    n_clique: int
        clique的数目
    size: int
        clique的(最大)大小
    stepped: bool
        生成的clique的大小递减(均匀分布在3-size上,不能整除的余数都生成3个点的)还是固定尺寸
    """
    np.random.seed(seed)
    # 参数重调
    if stepped:
        n = n_clique//(size-2)
        n3 = n_clique - n*(size-3)
        size_min = 3
        n_node = n3*3 + n*(4+size)*(size-3)/2
    else:
        n = n_clique
        size_min = size
        n_node = n_clique*3
    #选点
    vertices_choice = list(np.random.choice(len(G), n_node, replace = False)) # 先一次性把点都选出来

    for s in range(size_min,size+1):
        if stepped and s==3:
            num = n3
        else:
            num = n
        for _ in range(num):
            nodes = []
            for _ in range(s):
                v = vertices_choice.pop()
                for node in nodes:
                    if not G.has_edge(v,node):
                        G.add_edge(v,node)
                nodes.append(v)
    return G

def generate_graph_clique_wang(G_type,paremater,n,seed,num3,maxk):
    """生成全局rrg、局部clique的graph
    
    Parameters
    ----------
    base_degree: int
        RRG的节点度数
    num3: int
        为保证每种大小的clique边数一致,设定大小为k的clique的数目=int(6*num3/(k*(k-1)))
    maxk: int 
        最大clique的节点数
    """
    G = graph(G_type, paremater)
    for k in range(2,maxk+1):
        numk = int(6*num3/(k*(k-1)))
        np.random.seed(seed+k)
        for i_clique in range(numk):
            #在所有节点中任选k个点建立k阶clique
            vertices_choice = np.random.randint(low=0,high=n,size=(k))
            for i1 in range(len(vertices_choice)):
                for i2 in range(i1+1,len(vertices_choice)):
                    v1 = vertices_choice[i1]
                    v2 = vertices_choice[i2]
                    if v2 not in list(G.neighbors(v1)):
                        G.add_edge(v1,v2)

    #如果生成的图有多个连通分支，则通过在两个分支之间增加两条连边来使它们连通
    while nx.number_connected_components(G) != 1:
        cc1 = list(list(nx.connected_components(G))[0])
        cc2 = list(list(nx.connected_components(G))[1])
        v1 = cc1[0]
        v2 = cc2[0]
        if v2 not in list(G.neighbors(v1)):
            G.add_edge(v1,v2)
        if len(cc1) > 1 and len(cc2) > 1:
            v1 = cc1[-1]
            v2 = cc2[-1]
            if v2 not in list(G.neighbors(v1)):
                G.add_edge(v1,v2)
    
    return G

def safe_inv(v):
    inv = np.zeros_like(v)
    for i in range(len(v)):
        if v[i]>0:
            inv[i] = 1./v[i]
        else:
            inv[i] = 0.
    return inv
   
if __name__ == '__main__':
    # G,gname= graph('rrg',[200,3,1])
    G,gname= graph('USA')
    # G,gname = graph('random_tree',[150,35])
    # n = 12
    # k = 3
    # seed_clique = 48
    # G = add_clique(G,n,k,seed_clique)
    # G = cut(G)
    # print(len(G))
    # # t = time.time()
    # regions = get_Regions(G,3,3)
    # # print(time.time()-t)

    # print(len(regions))
    # for region in regions:
    #     print(list(region.graph))
    #     if region.subg:
    #         print('subregions:')
    #         for part in region.subregions:
    #             print(list(part))
    #         # print(part.edges())
    #     print()

    R = 4
    N = 4
    # t = time.time()
    region_dict = get_partition(G,R,N)
    # print(time.time()-t)

    for r in range(3,R+1):
        print('R='+str(r))
        for region in region_dict[r]:
            print(list(region.graph))
            # print(list(region.graph.edges()))
            if region.subg:
                print('subregions:')
                for part in region.subregions:
                    print(list(part))
            print()
    

