import numpy as np
import glob
import shutil
import math
import time
import random
from collections import defaultdict

class convert_uppertriangle:
    """
    numpy 2d matrix를 numpy 2d uppertriangle matrix로 변환 시킨다.
    """
    def __init__(self):
        self._matrix=np.array([])
        self._uppermatrix=np.array([])
        
    @property
    def matrix(self):
        return self._matrix
    
    @matrix.setter
    def matrix(self,matrix):
        if len(matrix.shape) != 2:
            raise Exception("The dimension of array must be 2")
        else:
            self._matrix=matrix

    @property
    def uppermatrix(self):
        return self._uppermatrix
    
    def convert(self):
        l=len(self._matrix)
        for i in range(l):
            tmp=self._matrix[i]
            for j in range(l):
                if i>j:
                    tmp[j]=0
            self._matrix[i]=tmp
        self._uppermatrix=self._matrix

        
def vec_operation_mod_p(v,w,mod): #a를 b에다 더함 주의 (방향성 있음)
    n=mod
    a=v.copy()
    b=w.copy()
    for i in a:
        if (i in b):
            if (b[i]-a[i])%(n) ==0:
                del b[i]
            else:
                b[i]=(b[i]-a[i])%(n)
        else:
            b[i]=a[i]
    return b

def col_operation_mod_p(matrix,i,j,mod): #i를 j에다 더함 주의 (방향성 있음)
    a={}
    i_col=matrix[i]
    for k in i_col:
        a[k]=1
    b={}
    j_col=matrix[j]
    for k in j_col:
        b[k]=1
    c=vec_operation_mod_p(a,b,mod)
    d=[]
    for k in c:
        d.append(k)
    d.sort()
    return(tuple(d))

def link1(v,w):
    a=v.copy()
    b=w.copy()
    if a[0]==b[0]:
        del a[0]
        a.reverse()
        c=a+b
    elif a[0]==b[-1]:
        del a[0]
        c=b+a
    elif a[-1]==b[0]:
        del a[-1]
        c=a+b
    elif a[-1]==b[-1]:
        del a[-1]
        a.reverse()
        c=b+a
    else:
        return 'None'
    return c

def cycle_link(W):
    if len(W)==1:
        return W
    else:
        L=W.copy()
        a=L[0]
        del L[0]
        i=0
        c=link1(a,L[0])
        if c=='None':
            while link1(a,L[i])=='None':
                i+=1
                c=link1(a,L[i])
            del L[i]
            L=[c]+L
            return cycle_link(L)
        else:
            del L[0]
            L=[c]+L
            return cycle_link(L)

def Caculate_Vietoris_Ripscomplex_Standardalgorithm(distmat,max_homology_dim):
    """
    distance matrix와 계산하길 원하는 max homology dim을 입력받아 standard algorithm에 따라 Rips complex를 계산해낸다.
    이 코드로 barcode와 cycle을 계산해 낼수있다.
    """

    dist_mat=distmat
    max_dim=max_homology_dim+1

    if max_dim < 1 or not isinstance(max_dim,int):
        raise Exception("max_dim must be at least 1 or int type")


    for i in range(max_dim+1):
        globals()['S_{}'.format(i)]=[]

    l_1=len(dist_mat)
    for i in range(l_1):
        S_0.append((i,))


    S_1_info=[]
    for i in range(len(dist_mat)):
        for j in range(i,len(dist_mat)):
            if not dist_mat[i][j]==0:
                S_1_info.append([dist_mat[i][j],i,j])
    S_1_info.sort()

    l_2=len(S_0)
    matrix_index_dict={(i,):i for i in range (l_2)}
    matrix_col_nonzero={i:() for i in range (l_2)} 
    simplex_value={(i,):0 for i in range (l_2)}


    col_loca=l_2

    adjacent_matrix=defaultdict(dict)
    for i in range(len(S_0)):
        tmp=defaultdict(dict)
        adjacent_matrix[i]=tmp


    for info in S_1_info:
        birth=info[0] 
        row_loca=edge=(info[1],info[2]) # =row_loca
        new=(edge,birth) 

        birth_list=[(new,row_loca)] 

        while 1:

            if len(birth_list)!=0:

                birth_object=birth_list[0][0] #new=(edge,birth)
                A=birth_object[0] 
                B=birth_object[1]            

                l_4=len(A)
                globals()[f'S_{l_4-1}'].append(tuple(A)) #S_1 에 edge 넣고, S_2에 face 넣고, ...

                simplex_value[tuple(A)]=B
                L=0
                for d in range(max_dim+1):
                    L+=len(globals()['S_{}'.format(d)])
                matrix_index_dict[A]=L-1 #이부분을 col_loca로 안하는 이유는?
                matrix_col_nonzero[col_loca]=birth_list[0][1] 


                col_loca+=1
                birth_list.remove(birth_list[0])    

                adjacent_check_member=[]
                for i in range(l_4):
                    tmp=list(A)
                    tmp.remove(A[i])
                    adjacent_matrix[A[i]][tuple(tmp)]='a'
                    adjacent_check_member.append(tuple(tmp))

                if max_dim >= l_4: 

                    tmp_birth_list=[]

                    for k in S_0:
                        j=k[0]

                        if j not in A:
                            check_value=0

                            for z in adjacent_check_member:
                                check_value+=len(adjacent_matrix[j][z])

                            if check_value==l_4: 
                                row_loca=[matrix_index_dict[A]]
                                new_member=[]
                                #j 점과 연결된 a, b 점들의 엣지 만들기
                                for i in range(l_4):
                                    a=list(A)
                                    c=a.copy()
                                    c[i]=j
                                    c.sort()
                                    new_member.append(tuple(c))

                                for mem in new_member:

                                    row_loca.append(matrix_index_dict[mem])
                                row_loca.sort()
                                new=list(A)+[j]
                                new.sort()
                                tmp_birth_list.append(((tuple(new),B),tuple(row_loca)))

                    birth_list=tmp_birth_list+birth_list
            else:
                break



    D=matrix_col_nonzero
    R=D.copy()
    V={}
    low_R={}
    l_5=len(D)
    for i in range(l_5):
        V[i]=(i,)

    for j in range(l_5):
        if len(R[j])==0:
            continue
        else:
            while len(R[j])!=0:
                if max(R[j]) not in low_R:
                    low_R[max(R[j])]=j
                    break
                else:
                    i=low_R[max(R[j])]
                    R[j]=col_operation_mod_p(R,i,j,2)
                    V[j]=col_operation_mod_p(V,i,j,2)

    matrix_col_to_simplex={v:k for k,v in matrix_index_dict.items()}  
    col_to_value={}

    for i in range(l_5):
        col_to_value[i]=simplex_value[matrix_col_to_simplex[i]]


    for i in range(max_dim):
        globals()['barcode{}'.format(i)]=[]
        globals()['cycle{}'.format(i)]=[]

    for j in range(l_5):
        if len(R[j])==0:
            S=matrix_col_to_simplex[j]
            l=len(S)-1
            birth=col_to_value[j]

            if l==0: #0th homo
                if (j in low_R):
                    death=col_to_value[low_R[j]]
                    if birth != death:
                        a,b=matrix_col_to_simplex[low_R[j]]
                        globals()[f'barcode{l}'].append([(birth,death),(a,),(b,)])
                        globals()[f'cycle{l}'].append([a,b])
                else:
                    globals()[f'barcode{l}'].append([(birth,'infty'),S])
                    globals()[f'cycle{l}'].append(list(S))
            elif 1<= l <= max_homology_dim:
                if (j in low_R):
                    death=col_to_value[low_R[j]]           
                    if birth != death:
                        tmp_cycle=[]
                        for k in V[j]:
                            tmp_cycle.append(matrix_col_to_simplex[k])
                        e=''
                        for s in tmp_cycle:
                            e+=f' + {list(s)}'
                        globals()[f'barcode{l}'].append([  (birth,death)]+tmp_cycle )
                        globals()[f'cycle{l}'].append(list(map(list,tmp_cycle)))
                else:
                    tmp_cycle=[]
                    for k in V[j]:
                        tmp_cycle.append(matrix_col_to_simplex[k])
                    e=''
                    for s in tmp_cycle:
                        e+=f' + {list(s)}'            
                    globals()[f'barcode{l}'].append([  (birth,'infty')]+tmp_cycle )
                    globals()[f'cycle{l}'].append(list(map(list,tmp_cycle)))
            else:
                continue
        else:
            continue
    
    return matrix_index_dict,matrix_col_nonzero,simplex_value,R,D,V, low_R


def cycle1_index(show=None):
    """
    index version의 cycle을 보기좋게 보여주는 코드
    """
    Cycles=[]
    for i in range(1,len(cycle1)+1):
        a=cycle_link(cycle1[i-1])
        b=a[0][:-1]
        Cycles.append(b)
        if show =='True':
            print('C_{} : {}'.format(i,b))
    return Cycles

def cycle1_eng(show=None):
    """
    english version의 cycle을 보기좋게 보여주는 코드
    """
    Cycles=[]
    for i in range(1,len(cycle1)+1):
        a=cycle_link(cycle1[i-1])
        cycle_tmp=a[0][:-1]
        cycle=[nodelist[j] for j in cycle_tmp]
        Cycles.append(cycle)
        if show ==True:
            print('C_{} : {}'.format(i,cycle))
    return Cycles

def Overlapmatrix(Cycles, allnode,scale,normalized=False):
    """
    Cycles, allnode, nodelist, scale을 정해주면 overlapmatrix를 구해주는 코드
    """
    c=Cycles
    cyc_len=[len(ci) for ci in c]
    cyc_len_max=max(cyc_len)
    matrix = []; 
    for ci in c: 
        row = []
        for node in allnode:
            if node in ci:
                row.append(1)
            else:
                row.append(0)
        matrix.append(row)
    q=scale
    mat = [];
    for j,row in enumerate(matrix):
        if normalized==True:
            q=scale*cyc_len[j]/cyc_len_max         
        new = []; leng = 0
        for i,node in enumerate(row):
            if node == 0:
                if leng == 0:
                    new.append(0)
                else:
                    if leng >= q:
                        new += row[i-leng:i]+[0]
                    else:
                        new += [0]*leng+[0]
                    leng = 0
            else:
                leng += 1
        if leng != 0:
            if leng >= q:
                new += row[len(row)-leng:]
            else:
                new += [0]*leng
        mat.append(new)
    
    return mat



###정보의 계산
def CalculateInfo(folder,songname,max_homology_dim,unit_beat_length=11520,inverse_edge_weight=True,weight_finding_path=True,vertex_weight=False):
    info=PrepareInfo(folder,songname)
    distmat,allnode, nodelist = analysis_from_info(info,unit_beat_length,inverse_edge_weight,weight_finding_path,vertex_weight)
    matrix_index_dict,matrix_col_nonzero,simplex_value,R,D,V, low_R=Caculate_Vietoris_Ripscomplex_Standardalgorithm(distmat,max_homology_dim)
    Cycles=cycle1_eng()
    overlapmatrix=Overlapmatrix(Cycles, allnode,scale=4)
    result={'songname':songname,'distmat':distmat,'allnode':allnode, 'nodelist':nodelist,'Cycles':Cycles,'overlapmatrix':overlapmatrix,'matrix_index_dict':matrix_index_dict,'matrix_col_nonzero':matrix_col_nonzero,'simplex_value':simplex_value,'R':R,'D':D,'V':V, 'low_R':low_R}
    
    return result


###정보의 표시
def Show_barcode(homology_dimension=None,expression='index'):
    """
    계산해낸 barcode를 보기 좋게 보여주는 코드
    """
    if homology_dimension==None:
        dim_list=[0,1,2]
    else:
        dim_list=[homology_dimension]

    if 0 in dim_list:
        print(f'\n 0th homology: \n')
        for i in range(len(barcode0)):
            L=barcode0[i]
            birth = round(L[0][0],2)
            death = L[0][1]
            if type(death)!=str:
                death = round(L[0][1],2)
            print(f'B{i+1}= [{birth},{death}) : {L[1:]}')        
    if 1 in dim_list:
        print(f'\n 1th homology: \n')
        for i in range(len(barcode1)):
            L=barcode1[i]
            birth = L[0][0]
            death = L[0][1]
            if expression == 'index':
                Cycles=cycle1_index()
            elif expression == 'english':
                Cycles=cycle1_eng()
            cycle=Cycles[i]
            print(f'B{i+1}= [{round(birth,2)},{round(death,2)}) : {cycle}')
    if 2 in dim_list:
        print(f'\n 2th homology: \n')
        for i in range(len(barcode2)):
            L=barcode2[i]
            birth = L[0][0]
            death = L[0][1]
            print(f'B{i+1}= [{round(birth,2)},{round(death,2)}) : {L[1:]}')  
