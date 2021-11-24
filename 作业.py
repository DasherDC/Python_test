import numpy as np
# potential
sigma=1
epsilon=1
rcut=2.5*sigma
rcut2=rcut*rcut
mass=1
k=1


#lattice
alat=1.54*sigma
nx=3
ny=3
nz=3
xbox=nx*alat
ybox=ny*alat
zbox=nz*alat
box=np.array([xbox,ybox,zbox])

def Apply_PBC(xij):
    xij=xij/box
    xij=xij-np.ceil(xij-0.5)
    xij=xij*box
    return xij

#FCC crystal
natoms=nx*ny*nz*4
pos=np.zeros((natoms,3))
force=np.zeros((natoms,3))

iatom=0
for i in range(nx):
    for j in range(ny):
        for k in range(nz):
            pos[iatom]=np.array([i,j,k])
            pos[iatom+1]=np.array([i+0.5, j+0.5, k])
            pos[iatom+2]=np.array([i+0.5, j, k+0.5])
            pos[iatom+3]=np.array([i, j+0.5, k+0.5])
            iatom +=4
            pos *=alat

#Initialization
kb=1
temperature =0.1
vel=np.random.random((natoms,3)) -0.5
vel=vel-np.sum(vel, axis=0)/natoms
kin=0.5*mass*np.sum(vel*vel)
scale=np.sqrt(1.5*natoms*kb*temperature/kin)
vel=vel*scale

#Subroutine to calculate force
def Force():
    global Epot
    global force
    Epot =0
    force *=0
    for i in range (natoms):
        for j in range (i+1, natoms):
            xij=pos[j]-pos[i]
            xij=Apply_PBC(xij)
            rij2=np.sum(xij*xij)
            if (rij2>rcut2):
                continue
            r=np.sqrt(rij2)
            Epot+=0.5*k*xij*xij/(r*r)
            fij=k*xij/r
            force[i]+=fij
            force[j]-=fij

#Themo
def Thermo(istep):
   Ekin=0.5*mass*np.sum(vel*vel)
   Etot=Epot+Ekin
   T=Ekin/(1.5*natoms*kb)   
   fd.write("{} {} {} {} {} \n".format(istep*dt,T,Epot,Ekin,Etot))    
    print("{0:.4f} {1:.6f} {2:.6f} {3:.6f} {4:.6f}".format (istep*dt, T, Epot, Ekin, Etot))
 
#Dump xyz file
def Dump(istep):
    fd.write("{}\nStep: {}, Epot: {}, Ekin:{}\n".format (natoms,istep,Epot, Ekin))
    for i in range (natoms):
         fd.write("X {} {} {} \n".format(pos[i][0],pos[i][1], pos[i][2]))


#Setup                                       
#- -control parameter
timestep=0.001
dt=timestep
hdt=0.5*0.001
iout=5
idump=10
nsteps=1000

#-0th step
Epot=0
Ekin=0
Force()
acc=force/mass

#--Log file 
fp=open("log.dat","w")
fp.write("# time temperature Epot Ekin Etot\n")
print ("# time  temperature Epot Ekin Eto")
Thermo(0)

#--dump file 
fd=open("dump.xyz","w")
Dump(0)

#Main Loop:Velocity verlet

for istep in range (1, nsteps+1):
    vel=vel+acc*hdt
    pos=pos+vel*dt
    Force()
    acc=force/mass
if(istep%iout==0):
    Thermo(istep)
if(istep%idump==0):
    Dump(istep)

#Clean-up
fp.close()
fd.close()

