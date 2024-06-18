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
    regions = get_Regions(G,R,N) 
    n = 0
    # print(len(regions))
    for region in regions:
        # print(len(region[0].edges()))
        n += len(list(region.graph))
        print(list(region.graph),list(region.graph.edges()))
        if region.subg:
            print('subregions:')
            for subg in region.subregions:
                # print(list(subg))
                print(list(subg),list(subg.edges()),len(subg))
            if region.intersec:
                print('intersections:')
                for intersec in region.intersections:
                    # print(list(intersec))
                    print(list(intersec),list(intersec.edges()),len(intersec))
        print()
    print(n)

def print_region_diction(G,R,N):
    regions_dict = get_Regions_diction(G,R,N) 
    # print(len(regions))
    for r in range(3,R+1):
        n = 0
        print('R='+str(r))
        for region in regions_dict[r]:
            n += len(list(region.graph))
            print(list(region.graph)) #,list(region.graph.edges())
            if region.subg:
                print('subregions:')
                for part in region.subregions:
                    print(list(part))
            print()
        print(n)

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

def get_Regions(G,R:int,N:int):
    G = cut(G)
    G_remain = deepcopy(G)
    Regions = []
    # Regions = dict() # 字典的key为3到R,对应相应R取值的Region变量
    # for r in range(3,R):
    #     Regions[r] = []
    edges = list(G_remain.edges())
    for e in edges:
        # print(e)
        if G_remain.has_edge(e[0],e[1]):
            # print(e)
            # t = time.time()
            g = Region_generator(G_remain,e,R)
            # print('Region',time.time()-t)
            # print(len(g))
            G_remain.remove_edges_from(list(g.edges())) # 不管怎样，g的所有边都可以从G_remain上移除了
            # print(list(g))
            if len(g) > 2:# 在一条边的基础上找到了其他区域
                # t = time.time()
                region = region_graph(g,R,N)
                Regions.append(region)
                # Regions = region.in_dict(Regions)
                # print('subRegion',time.time()-t)
            if len(G_remain)>0:
                G_remain = cut(G_remain,list(g))
    return Regions

def get_Regions_diction(G,R:int,N:int):
    G = cut(G)
    G_remain = deepcopy(G)
    Regions_dict = dict() # 字典的key为3到R,对应相应R取值的Region变量
    for r in range(3,R+1):
        Regions_dict[r] = []
    edges = list(G_remain.edges())
    for e in edges:
        # print(e)
        if G_remain.has_edge(e[0],e[1]):
            # print(e)
            # t = time.time()
            g = Region_generator(G_remain,e,R)
            # print('Region',time.time()-t)
            # print(len(g))
            G_remain.remove_edges_from(list(g.edges())) # 不管怎样，g的所有边都可以从G_remain上移除了
            # print(list(g))
            if len(g) > 2:# 在一条边的基础上找到了其他区域
                # t = time.time()
                region = region_graph(g,R,N)
                Regions_dict[R] = Regions_dict[R]+[region]
                if R>3:
                    Regions_dict = region.update_dict(Regions_dict)
                # print('subRegion',time.time()-t)
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
        boundary = []
        for node in set(subg).intersection(set_g):
            if len(subg[node])<len(g[node]) or len(subg[node]) ==1: # 还有可用的连接,或是在边界处
                boundary.append(node)
        # print(boundary)
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
                                if len(set(path_out).intersection(set_g)) == 2:#就只有开头结尾两个点，没意义
                                    path_out = paths_out.pop()
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

    if len(subg)>0:
        g_left = remove_empty_nodes(g_left,list(subg))

    return [subg,g_left]

class region_graph:
    def __init__(self,g,R,n): #pseudo 代表的是小于限定值R的时候的前置步骤，没有必要生成subg
        self.graph = g
        self.subg = (len(list(g))>n)
        self.R = R
        self.n = n
        if self.subg:
            self.divide_region()
            self.n_sg = len(self.subregions)
            self.get_macro_region()
        else:
            self.macro_g = nx.Graph()
            self.macro_g.add_node(0)
            regions_inside = []
            if self.R>3:
                subg_inside = []
                G_remain = deepcopy(self.graph)
                g_left = deepcopy(self.graph)
                edges = list(G_remain.edges())
                for e in edges:
                    if G_remain.has_edge(e[0],e[1]):
                        g = Region_generator(G_remain,e,self.R-1)
                        G_remain.remove_edges_from(list(g.edges()))
                        if len(g) > 2:
                            g_left.remove_edges_from(list(g.edges)) # 顺便把小region里的边去掉
                            region_inside = region_graph(g,self.R-1,self.n)
                            if region_inside.subg: # 把已经找好的subg选出来
                                for subg in region_inside.subregions:
                                    subg_inside.append(subg)
                            else:
                                subg_inside.append(region_inside.graph)

                            regions_inside.append(region_inside)
                        if len(G_remain)>0:
                            G_remain = cut(G_remain,list(g))
            self.regions_inside = regions_inside

    def divide_region(self): 
        # 先搜索有没有更小的Region：
        regions_inside = []
        subg_inside = []
        G_remain = deepcopy(self.graph)
        g_left = deepcopy(self.graph)
        if self.R>3:
            edges = list(G_remain.edges())
            for e in edges:
                if G_remain.has_edge(e[0],e[1]):
                    g = Region_generator(G_remain,e,self.R-1)
                    G_remain.remove_edges_from(list(g.edges()))
                    if len(g) > 2:
                        g_left.remove_edges_from(list(g.edges)) # 顺便把小region里的边去掉
                        region_inside = region_graph(g,self.R-1,self.n)
                        if region_inside.subg: # 把已经找好的subg选出来
                            for subg in region_inside.subregions:
                                subg_inside.append(subg)
                        else:
                            subg_inside.append(region_inside.graph)

                        regions_inside.append(region_inside)
                    if len(G_remain)>0:
                        G_remain = cut(G_remain,list(g))

        self.regions_inside = regions_inside
        # 相对于R-1新添加部分的subg
        subgs = []
        # 保留连通区域
        g_remain = remove_empty_nodes(g_left)
        while len(g_remain)>self.n:
            # 找初始边
            degree = dict(nx.degree(g_remain))
            min_degree_node = [node for node, deg in degree.items() if deg == min(degree.values())][0]
            if degree[min_degree_node] <= 2:
                neigh = list(g_remain[min_degree_node])
                neigh_degree = dict(nx.degree(g_remain,neigh))
                min_degree_neigh = [node for node, deg in neigh_degree.items() if deg == min(neigh_degree.values())][0]
                e = (min_degree_node,min_degree_neigh)
            else:
                min_weight = 2*len(list(self.graph))
                e = None
                for edge in list(g_remain.edges()):
                    [n1,n2] = edge
                    weight = degree[n1]+degree[n2]
                    if min_weight>weight:
                        e = edge
            # print(e,degree[e[0]]+degree[e[1]])
            # 从这条边开始生成subg
            temp = subregion_generator(g_remain,self.graph,e,self.R,self.n)
            subg,g_remain = temp

            if len(subg)>0:
                subgs.append(subg)
                if len(g_remain)<=self.n:
                    break
        if len(g_remain)>0:
            subgs.append(g_remain)

        #联通区域断开
        # g_remain = list(nx.connected_components(remove_empty_nodes(g_left)))
        # if len(g_remain)>1:
        #     for g in g_remain:
        #         if len(g)<=self.n:
        #             g_remain.remove(g)
        #             subgs.append(g)

        # while len(g_remain)>0:
        #     g = g_remain.pop()
        #     while len(g)>self.n:
        #         # 找初始边
        #         degree = dict(nx.degree(g_remain))
        #         min_degree_node = [node for node, deg in degree.items() if deg == min(degree.values())][0]
        #         if degree[min_degree_node] == 2:
        #             neigh = list(g_remain[min_degree_node])
        #             neigh_degree = dict(nx.degree(g_remain,neigh))
        #             min_degree_neigh = [node for node, deg in neigh_degree.items() if deg == min(neigh_degree.values())][0]
        #             e = (min_degree_node,min_degree_neigh)
        #         else:
        #             min_weight = 2*len(list(self.graph))
        #             e = None
        #             for edge in list(g_remain.edges()):
        #                 [n1,n2] = edge
        #                 weight = degree[n1]+degree[n2]
        #                 if min_weight>weight:
        #                     e = edge
        #         # print(e,degree[e[0]]+degree[e[1]])
        #         # 从这条边开始生成subg
        #         temp = subregion_generator(g,self.graph,e,self.R,self.n)
        #         subg,g = temp

        #         if len(subg)>0:
        #             subgs.append(subg)
        #             if len(g)<=self.n:
        #                 break
        #     if len(g)>0:
        #         subgs.append(g)

        # 重塑，去掉公共边(和intersection)之后的subregion尺寸变小，如果可能可以合在一起
        # for subg in subgs:
        #     print(list(subg),list(subg.edges()))
        # subgs.sort(key = lambda x:len(x))

        # reform_subg = []
        # unvisited_subg = list(range(len(subgs)))
        # while len(unvisited_subg)>0:
        #     subg_id = unvisited_subg[0]
        #     subg = subgs[subg_id]
        #     set_subg = set(subg)
        #     for j in unvisited_subg:
        #         if j != subg_id:
        #             g = subgs[j]
        #             if len(set(g).intersection(set_subg))>0:
        #                 new_g_temp = nx.compose(subg,g)
        #                 if len(new_g_temp)<=self.n:
        #                     # print(subg_id,j,list(new_g_temp))
        #                     subg = new_g_temp
        #                     set_subg = set(subg)
        #                     unvisited_subg.remove(j)
        #     reform_subg.append(subg)
        #     unvisited_subg.remove(subg_id)
        # subgs = reform_subg
        
        self.subregions = subgs + subg_inside

    def get_macro_region(self):
        self.macro_g = nx.empty_graph()
        l = self.n_sg
        self.macro_g.add_nodes_from([i for i in range(l)])
        self.boundaries = [[None for _ in range(l)] for _ in range(l)]
        for i in range(l):
            g1 = set(self.subregions[i])
            for j in range(i):
                g2 = set(self.subregions[j])
                boundary = g1.intersection(g2)
                if len(boundary)>0:
                    self.macro_g.add_edges_from([(i,j)])
                    self.boundaries[i][j] = list(boundary)
                    self.boundaries[j][i] = list(boundary)

    def update_dict(self,region_diction):
        # 还没改好，争取一劳永逸
        regions_inside = self.regions_inside
        while len(regions_inside)>0:
            region_inside = regions_inside.pop()
            region_diction[self.R-1] = region_diction[self.R-1] +[region_inside]
            if self.R>4:
                region_inside.update_dict(region_diction)
        return region_diction

# about graph

def graph(G_type,paremater = []):
    '''
    G_type
    '''
    gname = G_type
    
    if G_type == 'rrg':
        [n,c,seed] = paremater
        g = nx.random_regular_graph(c,n,seed=seed)
        gname += '_c=' + str(c) + '_seed=' + str(seed)
    elif G_type == 'erg': # 也即泊松度分布
        [n,c,seed] = paremater
        g = nx.erdos_renyi_graph(n,c/n,seed=seed)
        gname += '_c=' + str(c) + '_seed=' + str(seed)
    elif G_type == 'tree':
        [n,seed] = paremater
        g = nx.random_tree(n,seed = seed)
    elif G_type == 'powerlaw': # 无标度网络
        [n,m,seed] = paremater
        g = nx.barabasi_albert_graph(n,m,seed = seed)
        gname += '_m=' + str(m) + '_seed=' + str(seed)
    elif G_type == 'smallworld': # 小世界网络
        [n,c,p,seed] = paremater
        g = nx.watts_strogatz_graph(n,c,p,seed = seed)
        gname += '_c=' + str(c) + '_p=' + str(p) + '_seed=' + str(seed)
    elif G_type == 'macrotree':
        g = macrotree()
        n = len(list(g))
    elif G_type == 'qubics':
        num = paremater[0]
        g = qubic_graph(num)
        n = len(list(g))
    elif G_type == 'elist':
        [gname,elist] = paremater
        g = nx.empty_graph()
        g.add_edges_from(elist)
        n = len(list(g))
    elif G_type == 'club':
        gname = 'karate_club'
        g = nx.karate_club_graph()
        n = len(list(g))
    elif G_type == 'bcspwr':
        [num] = paremater
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
    G,gname= graph('rrg',[200,3,1])
    # G,gname= graph('macrotree')

    # # t = time.time()
    # regions = get_Regions(G,7,7)
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

    R = 7
    N = 7
    t = time.time()
    region_dict = get_Regions_diction(G,R,N)
    print(time.time()-t)

    for r in range(3,R+1):
        print('R='+str(r))
        for region in region_dict[r]:
            print(list(region.graph))
            if region.subg:
                print('subregions:')
                for part in region.subregions:
                    print(list(part))
            print()
