import networkx as nx
import numpy as np
import time

'''
重新理清了PA的本质，多个公共点的情况不能再用联合概率的方式了，反而是最早的message形式更合理,因此也没有必要划分subregion
'''

def print_region(G,R,N):
    partition = get_partition(G,R,N) 
    nodes = set()
    for region in partition[R]:
        nodes.update(list(region))
        print(list(region),list(region.edges()))
    print(len(nodes))

def print_region_diction(G,R,N,edges:bool):
    regions_dict = get_partition(G,R,N) 
    # print(len(regions))
    for r in range(3,R+1):
        all_nodes = []
        print('R='+str(r))
        for region in regions_dict[r]:
            all_nodes += list(region)
            print(list(region))
            if edges:
                print(list(region.edges()))
        print(len(set(all_nodes)))

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
        G_n = nx.Graph(G)
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
    G_environment = nx.Graph(G)
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
            add_path(Region,path)
            remove_path(G_environment,path)

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
    G_remain = nx.Graph(G)
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
                Regions_dict = partition_dict_update(Regions_dict,local_partition,R+1)
            if len(G_remain)>0:
                G_remain = cut(G_remain,list(g))
    return Regions_dict

def partition_dict_update(dict_R,dict_r,R):
    for r in range(3,R):
        dict_R[r] = dict_R[r] + dict_r[r]
    return dict_R

def split_region(g,N):
    regions = []
    g_left = nx.Graph(g)
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
                    add_path(region,path)
                    remove_path(g_left,path)
                else:
                    break
            regions.append(region)
            g_left = cut(g_left,list(region))

    return regions

def add_path(g,path):
    for node_id in range(len(path)-1):
        g.add_edge(path[node_id],path[node_id+1])
def remove_path(g,path):
    for node_id in range(len(path)-1):
        if g.has_edge(path[node_id],path[node_id+1]):
            g.remove_edge(path[node_id],path[node_id+1])

# claim:一个圈只要被断掉，那么不管分成几份都等同于用PA即'目->口+||+口 =口+|+|+口  or 口+凵+凵'其中后者反而更麻烦一些
def get_local_partition(g,R,N):
    local_dict = dict() # 字典的key为3到R,对应相应R取值的Region变量
    for r in range(3,R+1):
        local_dict[r] = []

    # 检索并删除R=R-1的region
    G_remain = nx.Graph(g) # 用于搜索子图的剩余区域，找不到region的边会被删掉，去掉子图后没有用的摇摆边(必然属于更长的圈)会被剪掉
    g_left = nx.Graph(g) # 只删除region的d真正的剩余区域
    if R>3:
        edges = list(G_remain.edges())
        for e in edges:
            if G_remain.has_edge(e[0],e[1]):
                region = Region_generator(G_remain,e,R-1)
                G_remain.remove_edges_from(list(region.edges))
                if len(region) > 2:
                    local_partition = get_local_partition(region,R-1,N)
                    local_dict = partition_dict_update(local_dict,local_partition,R)
                    for region in local_partition[R-1]:
                        g_left.remove_edges_from(list(region.edges)) 
                        # 把region小region里的边去掉,该操作会保留在小region中因为N的限制而本质上并没有被考虑的边
                if len(G_remain)>0:
                    G_remain = cut(G_remain,list(region))
    
    # R=R的部分
    if len(list(g))<=N:
        local_dict[R].append(g)
    else:
        # 检索新添加区域和在R-1部分被舍去的边构成的region
        g_left = remove_empty_nodes(g_left)# R新增的区域以及R-1时没有被考虑的部分
        print(list(g_left))
        regions_in_g_left = []
        G_parts = list(nx.connected_components(g_left)) # 经常会是不止一个连通区域
        for nodes in G_parts:
            G_remain = g_left.subgraph(nodes)
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
            # 搜索合并的方案：
            # 构建一个新的图：每一个R-1的region是一个节点，一条或多条边存在在节点和节点之间以及节点自身，节点有权重n，边有长度l
            # (n=N的节点可以省略，n+n_neighbour>N or all{n+n_neighbour+l_min+l_min2}>N的同理,后面两个放在循环过程中吸收自环之后可能更好)
            # 因为每一个圈的长度都是R，造成的误差是同一数量级的，因此假定不分先后
            # 首先试着吸收节点到自身的边，若n+l<=N则吸收之，吸收后若n+l=N则将其去除，否则继续尝试
            # 接着重复以下循环：
            # 接着吸收两点间的边：若n1+n2+l_min+l_min2<=N则吸收,吸收后n=N则将其去除。
            
            last_regions = []
            for region in (local_dict[R-1] + regions_in_g_left):
                if len(region)<N: # 还能继续添加点的regions
                    last_regions.append(region)
                else:
                    local_dict[R].append(region)
            m = len(last_regions)
            boundaries = [list(nx.intersection(region,g_left)) for region in last_regions] # 这些region对于g_left的边界

            # 自连接的合并：
            remove_id = []
            for i in range(m):
                region = last_regions[i]
                n = len(region)
                b = boundaries[i]
                l = len(b)
                if l>1:# 在遗留的边中没有自圈
                    connections = []
                    for j in range(l):
                        if g_left.has_node(b[j]):
                            for k in range(j+1,l):
                                connections += list(nx.all_simple_paths(g_left, b[j], b[k], N-n+2)) # 添加后必然超上限的直接去掉
                    connections = sorted(connections,key = len)
                    for path in connections:
                        g_test = nx.Graph(region)
                        add_path(g_test,path)
                        if len(g_test)<=N:
                            region = g_test
                            remove_path(g_left,path)
                    if len(region) == N: # 因为是先从R-1的region开始延伸，应该是不可能出现恰好重合，能合并的情况即len(g1+g2)=len(g1)
                        remove_id.append(i)
            
            remove_id = sorted(remove_id,reverse=True)# 得从后往前删
            for k in remove_id: 
                local_dict[R].append(last_regions[k])
                del last_regions[k]
                del boundaries[k]

            # 剩余的小region间的合并：     
            boundaries = [[node for node in list(region) if len(region[node])<len(g[node])] for region in last_regions] # 考虑邻居合并会出现公共点的情况，因此对于boundaries的定义是不一样的
            while len(last_regions)>1: # 反复尝试合并邻居
                m = len(last_regions)
                merged = False
                large_region_id = [] # 用来存储过程中发现的和所有邻居直接合并都太大的region
                for i in range(m):
                    for j in range(i+1,m):
                        n1,n2 = len(last_regions[i]),len(last_regions[j])
                        if len(set(list(last_regions[i])).update(list(last_regions[j])))>N:# 直接拼起来都不行直接跳过
                            continue
                        # 搜索g_left构成的连接(和公共点)
                        connection = []
                        for a in boundaries[i]:
                            for b in boundaries[j]:
                                if a==b:
                                    connection.append([])
                                elif nx.has_path(g_left,a,b):
                                    connection += list(nx.all_simple_paths(g_left, n1, n2, N-n1-n2+2))
                        
                        if len(connection)>1: # 有两条及以上的路径则试着合并
                            # 通过最短的两条路径合并
                            connection = sorted(connection,key = len)
                            g_test = nx.union(last_regions[i],last_regions[j])
                            add_path(g_test,connection[0])
                            add_path(g_test,connection[1])
                            if len(g_test)<=N: # 足够小
                                # 此处实际上是不严谨的，通过最短的两条边连接也超过N但是通过第三个Region却能合并是存在的，但是当前找到的最小反例也需要N=8，属于是无伤大雅的情况
                                # 现在连三区合并都不考虑了，就是当且仅当能合并的时候才保留
                                # No_small_neighbor = False
                                region = g_test
                                merged = True
                                remove_path(g_left,connection[0])
                                remove_path(g_left,connection[1])
                                # 尝试合并其他的路径(相当于是新region的自环),哪怕恰好是N也要试一下，避免第三小长度为2（直连）,因此==N没有单独处理
                                connection = connection[2:]
                                for path in connection:
                                    g_test = nx.Graph(region)
                                    add_path(g_test,path)
                                    if len(g_test)<=N:
                                        region = g_test
                                        remove_path(g_left,path)
                                    else: # 因为是从小到大排列的，所以短的不行长的更不行
                                        break
                                break
                    if merged:
                        break
                    else:# 跟所有都不能merge，反正下一个i不会考虑到此处的连接，所以放在最后删除
                        large_region_id.append(i)

                if merged: #合并相应变量
                    del last_regions[j]# 这里的j一定比i小也比所有的large_region_id小
                    del boundaries[j]
                    if len(region)==N:# 同上，不考虑极小概率的合并情况
                        large_region_id.append(i)
                    else:
                        last_regions[i] = region
                        boundaries[i] = [node for node in list(region) if len(region[node])<len(g[node])]

                # 每次开始下一次循环前把没朋友的region送进diction，为了确保merge的过程，需要在merge之后操作，因此把break单独放在后面
                large_region_id = sorted(large_region_id,reverse=True)# 得从后往前删，其实直接翻转好像就行，为了可读性还是这样吧
                for k in large_region_id: 
                    local_dict[R].append(last_regions[k])
                    del last_regions[k]
                    del boundaries[k]

                if not merged:# 遍历最新的last_regions后也没有成功
                    break
            
            # 不考虑三区合并的问题了，最小的情况也有N=7,完全是浪费时间
            local_dict[R] += last_regions # 直接加上就完事

            # # 邻居之间尝试合并完成，搭建剩余的不太大的region构成的graph用来搜索可多区块合并的结构
            # '''
            # 好麻烦啊...试试看常见的图到这一步复不复杂，实在不行就遍历
            # '''
            # last_regions = sorted(last_regions,key = len)# 从小到大排列
            # m = len(last_regions)
            # connections = [[False for _ in range(m)] for __ in range(m)]
            # region_graph = nx.Graph()
            # for i in range(m):
            #     for j in range(i+1,m):
            #         n1,n2 = len(last_regions[i]),len(last_regions[j])
            #         if len(set(list(last_regions[i])).update(list(last_regions[j])))>N:# 直接拼起来都不行直接跳过
            #             continue
            #         connection = []
            #         for a in boundaries[i]:
            #             for b in boundaries[j]:
            #                 if a==b:
            #                     connection.append([])
            #                 elif nx.has_path(G,a,b):
            #                     connection += list(nx.all_simple_paths(G, n1, n2, N-n1-n2+2))
            #         if len(connection) > 0:
            #             connections[i,j] = connection
            #             region_graph.add_edge(i,j)

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
    elif G_type == 'square_chain':
        num = parameter[0]
        g = qubic_graph(num)
        n = len(list(g))
    elif G_type == 'cross_square':
        num = parameter[0]
        g = cross_square(num)
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
    else:
        print('This graph is not defined in this function')

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

def cross_square(num):
    elist = [(0,1),(0,2),(1,3),(2,3)]
    for i in range(1,num):
        n = i*3
        elist += [(n,n+1),(n,n+2),(n+1,n+3),(n+2,n+3)]
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
    G,gname= graph('rrg',[200,3,1])
    # G,gname= graph('USA')
    # G,gname = graph('random_tree',[150,35])
    # n = 12
    # k = 3
    # seed_clique = 48
    # G = add_clique(G,n,k,seed_clique)
    # G = cut(G)
    # print(len(G))
    # G,gname= graph('cross_square',[5])

    R = 6
    N = 6
    t = time.time()
    region_dict = get_partition(G,R,N)
    print(time.time()-t)

    for r in range(3,R+1):
        print('R='+str(r))
        for region in region_dict[r]:
            print(list(region))
        print()
    

