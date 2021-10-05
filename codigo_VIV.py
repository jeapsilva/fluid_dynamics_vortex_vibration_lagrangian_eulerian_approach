# FORMULAÇÃO EF-LEA para HAGEN-POISSEUILLE.
# Autora: Jéssica Aparecida Silva
# Engenharia Mecânica - Universidade Federal do Rio de Janeiro

#------------------------------------------------------------------
#                                                         #
#           Bibliotecas utilizadas                        #
#                                                         #
#-------------------------------------------------------------------   

import numpy as np
import matplotlib.pyplot as plt
import meshio
import time

#------------------------------------------------------------------
#                                                         #
#           Seção com as funções utilizadas               #
#                                                         #
#-------------------------------------------------------------------   

def get_X():
    f = open("Nodes.txt","r")
    points = []
    for i in range(0,7359):
        line = f.readline()
        line = line.strip('\n')
        line = line.split(" ")
        line = line[1:]
        for i in range(0,3):
            line[i] = float(line[i])
        points.append(line)
    f.close()
    return points

def get_IEN():
    f = open("IEN.txt","r")
    points = []
    for i in range(0,14255):
        line = f.readline()
        line = line.strip('\n')
        line = line.split(" ")
        line = line[4:]
        for i in range(0,3):
            line[i] = int(line[i]) - 1
        points.append(line)
    f.close()
    return points
  
#------------------------------------------------------------------
#                                                         #
#           Seção para leitura de malha                   #
#                                                         #
#-------------------------------------------------------------------  
    
minY = 0
maxY = 2
minX = 0
maxX = 8
    
xc = (maxX-minX)/2.0
yc = (maxY-minY)/2.0

raio = 0.50
d = 0.01

IEN = np.array(get_IEN())
points = get_X()
npoints = len(points) # numero de pontos
ne = len(IEN)   

X = np.zeros(npoints,dtype='float') 
Y = np.zeros(npoints,dtype='float') 
for i in range(0,npoints):
    X[i] = points[i][0]
    Y[i] = points[i][1]
   
# Pontos de condição de contorno
cc = np.zeros((npoints,1),dtype='float')   # pontos que possuem uma condicao de contorno de psi (Top, Bottom, Inlet)
Fc = np.zeros((npoints,1),dtype='float')   # condicoes de contorno de psi inicial (Top, Bottom, Inlet)
circle = np.zeros((npoints,1),dtype='float')   # pontos que possuem uma condicao de contorno de psi (Top, Bottom, Inlet)
Fcircle = np.zeros((npoints,1),dtype='float')   # condicoes de contorno de psi inicial (Top, Bottom, Inlet)

ccoutlet = np.zeros((npoints,1),dtype='float')   # pontos que possuem uma condicao de contorno de psi no outlet
Fcoutlet = np.zeros((npoints,1),dtype='float') 

ccw = np.zeros((npoints,1),dtype='float')   # pontos de contorno para usar em omega

cc_velo = np.zeros((npoints,1),dtype='float')   # pontos que possuem uma condicao de contorno de velocidade (Top, Bottom, Inlet)
Fu = np.zeros((npoints,1),dtype='float')        # condições de contorno de u no topo e fundo 
Fv = np.zeros((npoints,1),dtype='float')         # codições de contorno de v no topo e fundo
Fyc = np.zeros((npoints,1),dtype='float')
Xc = np.zeros((npoints,1),dtype='float')  
Yc = np.zeros((npoints,1),dtype='float')  
fixo = np.zeros((npoints,1),dtype='float')
# Criando as matrizes de condição de contorno de Psi e de velocidades vx e vy
for i in range(0,npoints):                  
    # Verificando parte inferior do canal 
    if Y[i]== minY :
        fixo[i]=1.0
        cc[i] = 1.0
        Fc[i] = minY
        cc_velo[i] = 1.0
        Fu[i] = 0.0
        Fv[i] = 0.0
        ccw[i] = 1.0
    # Verificando parte superior do canal
    if Y[i] == maxY :
        fixo[i]=1.0
        cc[i] = 1.0
        Fc[i] = maxY
        cc_velo[i] = 1.0
        Fu[i] = 0.0
        Fv[i] = 0.0
        ccw[i] = 1.0
    # Verificando inlet 
    if X[i] ==  minX :
        fixo[i]=1.0
        cc[i] = 1.0
        Fc[i] = Y[i]
        cc_velo[i] = 1.0
        Fu[i] = 1.0   
        Fv[i] = 0.0
        ccw[i] = 1.0
    # Verificando outlet
    if X[i] == maxX:
        fixo[i]=1.0
        ccoutlet[i] = 1.0   #cc de psi
        Fcoutlet[i] = 0.0
        Fu[i]=0.0           #cc de velocidade
        Fv[i]=0.0
        cc_velo[i]=0.0
        ccw[i] = 1.0  #cc de omega
#        cc[i] = 1.0
#        Fc[i] = Y[i]
    # Verificando o cilindro
    if (raio -d) < ((xc-X[i])**2+ (Y[i]-yc)**2)**(0.5) < (raio + d) :
#        print(i)
        fixo[i]=1.0
        circle[i] = 1.0
        Fcircle[i] = yc
        cc_velo[i] = 1.0
        Fu[i] = 0.0   
        Fv[i] = 0.0
        Fyc[i] = yc
        Xc[i] = X[i]
        Yc[i] = Y[i]
        
#------------------------------------------------------------------
#                                                         #
#      Seção com criação das matrizes elementares         #
#                e matrizes globais                       #
#                                                         #
#------------------------------------------------------------------- 
        
# implementação do numero de iteracoes no tempo
nIter = 1000
t = 0
dt = 0.06
nu =0.001 

# Iniciando a velocidade        
vx = np.zeros( (npoints,1),dtype='float')
vy = np.zeros( (npoints,1),dtype='float')

    # ccs em vx e vy
for i in range(0,npoints):
    if cc_velo[i]==1.0:
        vx[i] = Fu[i]
        vy[i] = Fv[i]
        
#plt.figure(figsize=(12,6))
#plt.gca().set_aspect('equal')
#plt.triplot(X,Y,IEN)        
C1 = 1.5     
ISADORA = 0

GRAFICO=[]
DADOILA=[]
for j in range(0,nIter):  
     
    # inicio = time.time()        
    print('iteracao =',j) # visualizar a iteração 
    
    # inicializando as matrizes K e M, vetor condição de contorno F e Kest
    # e os gradientes na direção X e Y    
    K =  np.zeros( (npoints,npoints), dtype='float')
    M =  np.zeros( (npoints,npoints), dtype='float')
    Gx =  np.zeros( (npoints,npoints), dtype='float')
    Gy =  np.zeros( (npoints,npoints), dtype='float')
    Kest = np.zeros( (npoints,npoints), dtype='float')
    
    # loop dos elementos da malha
    for e in range(0,ne):
        v = IEN[e]
        det = X[v[2]]*( Y[v[0]]-Y[v[1]]) + X[v[0]]*( Y[v[1]]-Y[v[2]]) + X[v[1]]*(-Y[v[0]]+Y[v[2]])
        area = det/2.0
     
        # matriz de massa do elemento    
        m = (area/12.0) * np.array([ [2.0, 1.0, 1.0],
                                    [1.0, 2.0, 1.0],
                                    [1.0, 1.0, 2.0] ])        
        b1 = Y[v[1]]-Y[v[2]]
        b2 = Y[v[2]]-Y[v[0]]
        b3 = Y[v[0]]-Y[v[1]]
    
        c1 = X[v[2]]-X[v[1]]
        c2 = X[v[0]]-X[v[2]]
        c3 = X[v[1]]-X[v[0]]
        
        # matriz do gradiente
        B = (1.0/(2.0*area)) * np.array([ [b1, b2, b3],
                                          [c1, c2, c3] ])        
        # matriz do divergente
        BT = B.transpose()        
        # matriz de rigidez do elemento
        ke = area*np.dot(BT,B)
        
        gxe = (1.0/6.0)*np.array([ [b1, b2, b3],
                                   [b1, b2, b3],
                                   [b1, b2, b3] ])
        gye = (1.0/6.0)*np.array([ [c1, c2, c3],
                                   [c1, c2, c3],
                                   [c1, c2, c3] ])
        #   matriz de estabilização Kest
             #criacao da velocidade media
        v1 = IEN[e,0]
        v2 = IEN[e,1]
        v3 = IEN[e,2]
        
        vx_medio = ( vx[v1] + vx[v2] + vx[v3] )/ 3.0
        vy_medio = ( vy[v1] + vy[v2] + vy[v3] )/ 3.0
    
        kestx = ((dt/2.0)*(vx_medio/4*area))*np.array([ [vx_medio*b1*b1 + vy_medio*b1*c1 , vx_medio*b1*b2 + vy_medio*b1*c2 , vx_medio*b1*b3 + vy_medio*b1*c3],
                                                        [vx_medio*b2*b1 + vy_medio*b2*c1 , vx_medio*b2*b2 + vy_medio*b2*c2 , vx_medio*b2*b3 + vy_medio*b2*c3],
                                                        [vx_medio*b3*b1 + vy_medio*b3*c1 , vx_medio*b3*b2 + vy_medio*b3*c2 , vx_medio*b3*b3 + vy_medio*b3*c3] ])
        
        kesty = ((dt/2.0)*(vy_medio/4*area))*np.array([ [vx_medio*c1*b1 + vy_medio*c1*c1 , vx_medio*c1*b2 + vy_medio*c1*c2 , vx_medio*b1*b3 + vy_medio*c1*c3],
                                                        [vx_medio*c2*b1 + vy_medio*c2*c1 , vx_medio*c2*b2 + vy_medio*c2*c2 , vx_medio*b2*b3 + vy_medio*c2*c3],
                                                        [vx_medio*c3*b1 + vy_medio*c3*c1 , vx_medio*c3*b2 + vy_medio*c3*c2 , vx_medio*b3*b3 + vy_medio*c3*c3] ])         
        # matrizes globais
        for i in range(0,3):
            ii = IEN[e,i]
            for j in range(0,3):
                jj = IEN[e,j]
                K[ii,jj] = K[ii,jj] + ke[i,j]
                M[ii,jj] = M[ii,jj] + m[i,j]
                Gx[ii,jj] = Gx[ii,jj] + gxe[i,j]
                Gy[ii,jj] = Gy[ii,jj] + gye[i,j]   
                Kest[ii,jj] = Kest[ii,jj] + kestx[i,j] + kesty[i,j]     

    # CALCULO DA VELOCIDADE DA MALHA   
    R =[]    
    for NODE in range(0,npoints):
        V=[]
        for i in range(0,ne):                   
            if IEN[i,0] == NODE :
                V.append(IEN[i,1])
                V.append(IEN[i,2])
            if IEN[i,1] == NODE :
                V.append(IEN[i,0])
                V.append(IEN[i,2])
            if IEN[i,2] == NODE :
                V.append(IEN[i,0])
                V.append(IEN[i,1])
        V = list(set(V)) 
        R.append(V)
        
    p = np.zeros((npoints,1),dtype = 'float')          
    ui = np.zeros((npoints,1), dtype = 'float')
    vi = np.zeros((npoints,1), dtype = 'float')
    
    for i in range(0,len(R)):    
        VIZINHOS = R[i]
        sum_xj = 0.0
        sum_yj = 0.0
        sumdist = 0.0           
        for j in VIZINHOS:   
            dist = np.sqrt(((X[j]-X[i])**2)+((Y[j]-Y[i])**2))
            sumdist += dist
            sum_xj += X[j]*dist
            sum_yj += Y[j]*dist   
        ui[i] = c1*((sum_xj/sumdist)-X[i])/dt
        vi[i] = c1*((sum_yj/sumdist)-Y[i])/dt
        
#        # plotar a malha antes        
#        plt.figure(figsize=(14,8))
#        plt.gca().set_aspect('equal') 
#        plt.triplot(X.reshape(npoints),Y.reshape(npoints),IEN,color='blue')      
   
    for i in range(0,npoints): 
        if fixo[i] == 0.0:                              
            Y[i] = Y[i] + vi[i]*dt
            X[i] = X[i] + ui[i]*dt
        if circle[i]==1.0:
            Y[i] = Y[i] - (0.005)*np.sin(1.1*ISADORA)
            Fcircle[i] = Fcircle[i] + vi[i]*dt
            
    # R2D2 = np.sin(ISADORA*3) 
    # print(R2D2)   
    GRAFICO.append((0.05)*np.sin(1.1*ISADORA))/3
    DADOILA.append(ISADORA)
    ISADORA += dt     
    
#        # plotar a malha depois de mover
#        plt.triplot(X.reshape(npoints),Y.reshape(npoints),IEN,color='red') 
#        plt.savefig('Malha'+str(t)+'.png') 
#        plt.show()    
        
    # Calculando Gxvy - Gyvx
    b = np.dot(Gx,vy) - np.dot(Gy,vx)
    omega = np.linalg.solve(M,b)
    omegacc = omega.copy()   
            
    # Aplicando as condições de contorno de w 
    for i in range(0,len(cc)):
        if ccw[i] == 1.0:    
            omegacc[i] = omega[i] # so existe w na parede, no restante do dominio é nulo
            
    # nova velocidade corrigida
    unew = vx - ui      
    vnew = vy - vi    
      
    # Calculo de v.\nabla\omega
    VGO = np.diagflat(unew)*Gx + np.diagflat(vnew)*Gy
    
    # Iniciando solver transporte da vorticidade para omega n+1            
    LHSw = (1.0/dt)*M.copy()
    
    # Vetor do lado direito para eq. de transporte da vorticidade (omega)
    RHSw = (1.0/dt)*np.dot(M.copy(),omega) -np.dot(VGO,omega) + nu*np.dot(K.copy(),omega) #+ np.dot(Kest,omega)
    
    # imposicao das ccs para omega
    for i in range(0,len(cc)):
        if ccw[i] == 1.0:
            LHSw[i,:] = 0.0 
            LHSw[i,i] = 1.0 
            RHSw[i] = omegacc[i] 

    # SOLVER eq. transporte 
    omega = np.linalg.solve(LHSw,RHSw)
    
    # Funcao corrente K \psi = M* \omega
    RHSpsi = np.dot(M,omega)
    
    # Lado esquerdo da eq K\psi = M*\omega
    LHSpsi = K.copy()
    
    # Imposicao das ccs de \psi
    for i in range(0,len(cc)):
        if ccoutlet[i]==1.0:
            RHSpsi[i] = Fcoutlet[i]
        if cc[i] == 1.0:    
            LHSpsi[i,:] = 0.0 
            LHSpsi[i,i] = 1.0 
            RHSpsi[i] =  Fc[i] 
        if circle[i] == 1.0:
            LHSpsi[i,:] = 0.0 
            LHSpsi[i,i] = 1.0 
            RHSpsi[i] =  Fcircle[i] 
          
    # SOLVER da funcao corrente
    psi = np.linalg.solve(LHSpsi,RHSpsi)
    
    # encontrar as novas velocidades Mvx = Gy\psi, Mvy = -Gx\psi
    b_3 = np.dot(Gy,psi)   
    M3 = M.copy()
    vx = np.linalg.solve(M3,b_3)
    
    b_4 = np.dot(Gx,psi)
    M4 = M.copy()          
    vy = -np.linalg.solve(M4,b_4)

    # impor cc de velocidade nos vetores vx e vy
    for i in range(0,len(cc)):
        if cc[i] == 1.0:    
            vx[i] = Fu[i]
            vy[i] = Fv[i]
            
    # fim = time.time()
            
   
    msh = meshio.read('malhanova.msh')    
    IEN = msh.cells_dict['triangle']
    cells=[('triangle',IEN)]  
    
    xyz = msh.points
    
    for i in range(0,npoints):
        xyz[i]=[X[i],Y[i],0.0]
        
    # print('Tempo da interação:',(fim - inicio))
    # print(omega.max(),omega.min())  
    
    meshio.write_points_cells(               
        "./vtk/teste"+str(t)+".vtk",
        xyz,
        cells,
        # Optionally provide extra data on points, cells, etc.
         point_data={'vx':vx,'vy':vy,'psi':psi, 'omega':omega},
        # cell_data=cell_data,
        # field_data=field_data
        )
    
    # # GRAFICO 1 <j< nIter
    # #     Vorticidade Wz 
    # plt.figure(figsize=(12,6))
    # plt.gca().set_aspect('equal')
    # plt.tricontourf(X, Y,IEN, omega.reshape(npoints) , cmap='jet')
    # plt.colorbar()
    # plt.title('Vorticidade \n' )
    # plt.savefig('./plots/vortfig/vort'+str(t)+'.png' ) 
    # plt.clf()
            
    # #     Função corrente PSI 
    # plt.figure(figsize=(12,6))
    # plt.gca().set_aspect('equal')
    # plt.tricontourf(X, Y,IEN, psi.reshape(npoints) ,cmap='jet')
    # plt.colorbar()
    # plt.title('Função corrente \n')
    # plt.savefig('./plots/psifig/psi'+str(t)+'.png')
    # plt.clf()
    
    # #     Velocidade vx
    # plt.figure(figsize=(12,6))
    # plt.gca().set_aspect('equal')
    # plt.tricontourf(X, Y,IEN, vx.reshape(npoints) ,cmap='jet')
    # plt.colorbar()
    # plt.title('Velocidade vx \n')  
    # plt.savefig('./plots/vxfig/vx'+str(t)+'.png') 
    # plt.clf()
    
    # #     Velocidade vy
    # plt.figure(figsize=(12,6))
    # plt.gca().set_aspect('equal')
    # plt.tricontourf(X, Y,IEN, vy.reshape(npoints) ,cmap='jet')
    # plt.colorbar()
    # plt.title('Velocidade vy \n')
    # plt.savefig('./plots/vyfig/vy'+str(t)+'.png') 
    # plt.clf()
    
    t += 1

