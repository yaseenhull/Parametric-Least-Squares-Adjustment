# -*- coding: utf-8 -*-
"""
Created on Mon Jun 12 19:09:15 2017

@author: Yaseen Hull

Parametric Adjustment and Error ellipse calculation
"""
import sympy as sp
import numpy as np
import networkx as nx
from Functions import polar, polar1,edge_S, edge_D, edge_sD, edge_sS, get_X, get_Y, joinDir, joinDist
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
from networkx.drawing.nx_agraph import graphviz_layout
from scipy import stats




rho = 206264.8062

"""Read in points"""
""""""""""""""""""""
f2 = open("points.csv","r")
first2 = f2.readline()
#create points library
G = nx.DiGraph()

pnts = {}
unknowns = []
fixed = []
free = []

for line in f2:
    spl2 = line.split(";")
    name = spl2[0]
    Y = spl2[1]
    X = spl2[2]
    Type = spl2[3]
    
    if Type == '1':
        pnts[name]= [float(Y),float(X),True]
        fixed.append(name)
        
        
    elif Type == '0':
        pnts[name]= [float(Y),float(X),False]
        unknowns.append('Y'+name)
        unknowns.append('X'+name)
        free.append(name)

weights = []
obs = []

"""Read in observations """
""""""""""""""""""""""""""""""
f = open("data2.csv","r")
first = f.readline()

for line in f: 
    
    spl = line.split(";")
    From = spl[0]
    To = spl [1]
    Dist = (spl[2]) 
    Dir = (spl[3])
    Sig_Dist = (spl[4])
    Sig_Dir = (spl[5])
    Type = float(spl[6])
    
    if (Type == 1):
        obs.append([1,From, To,float(Dir), float(Sig_Dir)])
        G.add_edge(From,To, {'Direction':float(Dir),'Sigma Dir':float(Sig_Dir),'Type':1})
        weights.append(float(Sig_Dir))
        
    else:
        obs.append([0, From, To, float(Dist), float(Sig_Dist)])
        G.add_edge(From,To, {'Distance':float(Dist),'Direction':float(Dir), 'Sigma Dist': float(Sig_Dist),'Sigma Dir':float(Sig_Dir),'Type':0})
        weights.append(float(Sig_Dist))
      
P = np.diag(weights)

"""Read in orientation correction"""
""""""""""""""""""""""""""""""""""""
zlist =[]

f3 = open("Ori.csv","r")
firstline = f3.readline()

for line in f3:
    spl3 = line.split(";")
    zcorr = spl3[1]
    zlist.append(float(zcorr))
    unknowns.append(spl3[0])
        


"""Least Squares class - calculates partials for A matrix and misclosures for L matrix"""
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
        
class LeastSquares:
    def __init__(self):
        pass

    def distance(self, pi, pj):
        # pi and pj of the for [x, y, bool]
        dy = pj[0] - pi[0]
        dx = pj[1] - pi[1]
        return np.sqrt(dx * dx + dy * dy)
        
    def distanceEqn(self, pi, pj, obsv): #obs
        sijo = self.distance(pi, pj)       

        # calculate partial differentials
        dsdxi = dsdyi = dsdxj = dsdyj = None
        if not pi[2]:
            dsdyi, dsdxi = (pi[0] - pj[0]) / sijo, (pi[1] - pj[1]) / sijo

        if not pj[2]:
            dsdyj, dsdxj = (pj[0] - pi[0]) / sijo, (pj[1] - pi[1]) / sijo

        # calculate l - lo
        dl = obsv - sijo

        return [dsdxi, dsdyi, dsdxj, dsdyj, dl]
        
    def direction(self,pi,pj,obsv,zi):
        sijo = self.distance(pi,pj)
        
        dy = pj[0] - pi[0]
        dx = pj[1] - pi[1]
        
        sijo2 = np.arctan2(dy,dx)
        
        if sijo2<0:
            sijo2 = sijo2 + 2*np.pi
     
        drdxi = drdyi = drdxj = drdyj =  None
        
        if not pi[2]:
            drdxi, drdyi  = rho*((dy)/(sijo**2)), rho*(-(dx)/(sijo**2))
            #print(drdyi,drdxi)
        
        if not pj[2]:
            drdxj, drdyj = rho*(-(dy)/(sijo**2)), rho*((dx)/(sijo**2))
            #print(drdyj, drdxj)
        
        rad = rho*np.radians(obsv)
        
        for i in range(len(zlist)+1, len(unknowns)):
            if zi == unknowns[i]:
                dl = (rad) - rho*(sijo2) - rho*np.radians(zlist[i-6])
            else:
                continue
        
        return drdxi, drdyi, drdxj, drdyj, dl
            
ls = LeastSquares()

"""Iteration and indexing of values in Design matrix and Misclosure Matrix"""
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

for itera in range(1):
    A_list = []
    L_list = []
    num_unknown = len(unknowns)
    for o in obs:
        name_pi = o[1] #FROM
        name_pj = o[2] #TO
        Type = o[0]
        pi = pnts[name_pi] #TO -pi[0] = Y, pi[1] = X, pi[3] = bool
        pj = pnts[name_pj] #c3
        obsv = o[3]
        zi = 'z'+name_pi
        #zj = 'z'+name_pj
        
        rowA =[0.0] * num_unknown
        #print(name_pi,name_pj)
        
        if Type==1: #direction
            drdxi, drdyi, drdxj, drdyj, dl = ls.direction(pi, pj, obsv, zi)
            #print(ls.direction(pi,pj))           
            add_row = False
            
            if drdxi != None:
                index = unknowns.index('X'+name_pi)
                rowA[index] = drdxi
                add_row = True
                                
            if drdyi != None:
                index = unknowns.index('Y'+name_pi)
                rowA[index] = drdyi
                add_row = True
                
                
            if drdxj != None:
                index = unknowns.index('X'+name_pj)
                rowA[index] = drdxj 
                add_row = True
                                
            if drdyj != None:
                index = unknowns.index('Y'+name_pj)
                rowA[index] = drdyj
                add_row = True
                
            if zi in unknowns:
                #print('true') 
                index = unknowns.index('z'+name_pi)
                rowA[index] = -1
                add_row = True
                            
            if add_row:
                A_list += [rowA]
                L_list += [[dl]]
                
            else:
                A_list +=[rowA]
                L_list += [[dl]] 
                               
        
        elif o[0]==0:#
            
                       
            dsdxi, dsdyi, dsdxj, dsdyj, dl = ls.distanceEqn(pi, pj, obsv)
        #print(dsdxi, dsdyi, dsdxj, dsdyj, dl)   
            add_row = False
            
            if dsdxi != None:
                index = unknowns.index('X'+name_pi)
                rowA[index] = dsdxi
                add_row = True
        
            if dsdyi != None:
                index = unknowns.index('Y'+name_pi)
                rowA[index] = dsdyi
                add_row = True
        
            if dsdxj != None:
                index = unknowns.index('X'+name_pj)
                rowA[index] = dsdxj
                add_row = True
        
            if dsdyj != None:
                index = unknowns.index('Y'+name_pj)
                rowA[index] = dsdyj
                add_row = True
            
            if add_row:
                A_list += [rowA]
                L_list += [[dl]]
                
            else:
                A_list += [rowA]
                L_list += [[dl]]
                
    A = np.matrix(A_list)
    L = np.matrix(L_list)
    
    x = ((A.T) *P* A).I * (A.T) *P* L
    """Update unknowns"""
    
    for i in range(len(zlist)+1, len(unknowns)): 
        zn = zlist[i-6]+x[i,0]
        zlist[i-6] = zn

    for i in range(0,len(unknowns)-len(zlist)):
        if i%2 == 0:
            pnts[unknowns[i][1:]][0] = pnts[unknowns[i][1:]][0] + x[i,0]
        else:
            pnts[unknowns[i][1:]][1] = pnts[unknowns[i][1:]][1] + x[i,0]


"""Calculate residuals of observations"""
V = A*x-L

"""Calculate a posteriori poulation variance from covariance matrix"""

dof = len(A) - len(unknowns)
aPostriori = (V.T*P*V)/dof
var = np.float(aPostriori)
Ql = np.linalg.inv(P)
Qx = ((A.T*P*A).I)
Qv = Ql - A*Qx*A.T

El = var*Ql #covariance of obs a priori
Ex = var*Qx #covariance of unknowns
EL = var*A*Qx*A.T #covariance of adjusted observations
Ev = var*Qv #covariance of residuals

"""Index covariance matrix to calculate semi-major, semi-minor axis and orientation of ellipse"""
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
Eqn = []
elips = []
semiA = []
semiB = []

for i in range(Ex.shape[0]- len(zlist)):
    for j in range(Ex.shape[0]-len(zlist)):
        if i==j and i%2==0:
            varx = Ex[i,j]
            varxy= Ex[i,j+1]
            varyx = Ex[i+1,j]
            vary= Ex[i+1,j+1]
            """if varx<0:
                varx = varx**2
            if vary<0:
                vary = vary**2
            if varxy<0:
                varxy = varxy**2"""
            #print(varx,varxy,vary)
            a = 'a'+str(i)+str(j)
            b = 'b'+str(i)+str(j)
            alp = 'alpha'+str(i)+str(j)
            
            Eqn.append(varx)
            Eqn.append(vary)
            Eqn.append(varxy)
            Eqn.append(varyx)
            
            semiA.append((varx+vary+(((varx-vary)**2)+4*varxy*varxy)**0.5)/2)
            semiB.append((varx+vary-(((varx-vary)**2)+4*varxy*varxy)**0.5)/2)
            theta = np.arctan2((-2*varxy),(varx-vary))*0.5
            
            if theta < 0:
                theta = theta + 2*np.pi
            
            elips.append((theta))

nw = nx.MultiDiGraph()
#nx.draw_networkx(G)
plt.figure()
ax = plt.gca()
   
"""Plots error ellipse and draws network of nodes"""
""""""""""""""""""""""""""""""""""""""""""""""""""""""
for i in range(len(semiA)):
    
    ells = Ellipse(xy=[pnts[unknowns[i*2][1:]][1],pnts[unknowns[i*2][1:]][0]], width= semiB[i]*5000000, height= semiA[i]*5000000, angle= np.degrees(elips[i]), edgecolor='black', fc='white', lw=2 )
    ax.add_patch(ells)
    #print(pnts[unknowns[i*2][1:]][1],pnts[unknowns[i*2][1:]][0])
    pos = (pnts[unknowns[i*2][1:]][1],pnts[unknowns[i*2][1:]][0])
    nw.add_node(unknowns[i*2][1:],pos = (pnts[unknowns[i*2][1:]][1],pnts[unknowns[i*2][1:]][0]))
    
    #print(elips[i])
pos = nx.get_node_attributes(nw,'pos')
nw.add_edge('SUR10','RU4A')
nw.add_edge('RU4A', 'SUR11')
#nx.draw(G,pos, node_size = 30, node_color = 'blue', with_labels = True)
nx.draw_networkx_nodes(nw,pos,node_size =10, node_color = 'black')
nx.draw_networkx_edges(nw,pos,edgelist = nw.edges())
nx.draw_networkx_labels(nw,pos,font_size=10, font_family='arial', font_color='red')
plt.axis('scaled')
plt.show()


"""Plot network of nodes/points"""
""""""""""""""""""""""""""""""""""""
posit=nx.spring_layout(G)
#Nodes:
nx.draw_networkx_nodes(G,posit,nodelist = fixed,node_size=500,node_color='white',node_shape='^')
nx.draw_networkx_nodes(G,posit,nodelist = free,node_size=500,node_color='white')
#Edges:
nx.draw_networkx_edges(G,posit,edgelist=G.edges(),width=1, height=1)
#Labels:
nx.draw_networkx_labels(G,posit,font_size=10,font_family='calibri',font_color='red')
                  
plt.axis('off')
plt.savefig("traverse_network.png") # save as png
plt.show() # display""" 

"""Chi-square confidence interval"""

chi_left = stats.chi2.isf(q=0.025,df=11)
chi_right = stats.chi2.isf(q=0.975,df=11)
sigma_left = (dof*var)/chi_left
sigma_right = (dof*var)/chi_right



"""Output to file"""
""""""""""""""""""""        
wr = open('Output1.txt', 'w')
wt = open('Output2.txt', 'w')

print('Final Adjusted unknown points:')

wr.write('Design matrix' + '\n')
wr.write(str(np.round(A,3)) + '\n')
wr.write("\n ")
wr.write('Vector of misclosure'+'\n')
wr.write(str(np.round(L,3))+'\n')
wr.write("\n ")
wr.write('Weight matrix'+'\n')
wr.write(str(np.round(P,3))+'\n')
wr.write("\n ")
wr.write('Adjusted coordinates'+'\n')

for key in pnts:
    if pnts[key][2] != True:
        pnts[key][0] = round(pnts[key][0],2)
        pnts[key][1] = round(pnts[key][1],2)
        wr.write(key+'\t'+str(pnts[key][0]) +'\t'+str(pnts[key][1])+'\n')
        print(key+": "+'Y: '+str(pnts[key][0])+" "+ 'X: '+ str(pnts[key][1])+'\n')

print('Confidence interval of population variance'+'\n')
print(str(sigma_left)+"< sigma^2 <"+str(sigma_right))

wr.write("\n ")

wr.write('Residuals'+'\n')
wr.write(str(np.round(V,3)))
wr.write("\n ")

wr.write('Adjusted observaations' + '\n')
wr.write(str(L + V))
wr.write("\n ")

wr.write('Covariance of the observations a priori' + '\n')
wr.write(str(El)+'\n')
wr.write('\n')

wr.write('Covariance matrix of adjusted unknowns'+'\n')
wr.write(str(Ex))
wr.write("\n ")

wr.write('Covariance matrix of observations a posteriori'+'\n')
wr.write(str(EL)+'\n')
wr.write('\n')

wr.write('Covariance matrix of the residuals'+'\n')
wr.write(str(Ev)+'\n')
wr.write('\n')

wr.close()


wt.write('Solution vector x'+'\n')
wt.write(str(x)+'\n')
wt.write("\n ")
wt.write('a Posteriori variance factor'+ '\n')
wt.write(str(var)+ '\n')
wt.write("\n ")
wt.write('Confidence Interval for population variance' + '\n')
wt.write(str(sigma_left)+"< sigma^2 <"+str(sigma_right)) 
wt.close()


