import os
import math
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


def surface(Bpoint1, Bpoint2, Bpoint3):
    set=[]
    point=[0]*3
    point=[Bpoint1,Bpoint2,Bpoint3]
    vector1=[0]*3
    vector2=[0]*3
    vector1[0]=Bpoint2[0]-Bpoint1[0]
    vector1[1]=Bpoint2[1]-Bpoint1[1]
    vector1[2]=Bpoint2[2]-Bpoint1[2]
    unitvector1=(vector1[0]**2+vector1[1]**2+vector1[2]**2)**0.5

    vector2[0]=Bpoint3[0]-Bpoint1[0]
    vector2[1]=Bpoint3[1]-Bpoint1[1]
    vector2[2]=Bpoint3[2]-Bpoint1[2]
    unitvector2=(vector2[0]**2+vector2[1]**2+vector2[2]**2)**0.5
    for i in range(int(unitvector1)):
        x= (i+1)/unitvector1*vector1[0]+Bpoint1[0]
        y= (i+1)/unitvector1*vector1[1]+Bpoint1[1]
        z= (i+1)/unitvector1*vector1[2]+Bpoint1[2]
        set.append([x,y,z])
    for i in range(int(unitvector2)):
        x= (i+1)/unitvector2*vector2[0]+Bpoint1[0]
        y= (i+1)/unitvector2*vector2[1]+Bpoint1[1]
        z= (i+1)/unitvector2*vector2[2]+Bpoint1[2]
        set.append([x,y,z])
    vector3=[0]*3
    vector3[0]=Bpoint3[0]-Bpoint2[0]
    vector3[1]=Bpoint3[1]-Bpoint2[1]
    vector3[2]=Bpoint3[2]-Bpoint2[2]
    unitvector3=(vector3[0]**2+vector3[1]**2+vector3[2]**2)**0.5
    for i in range(int(unitvector3)):
        x= (i+1)/unitvector3*vector3[0]+Bpoint2[0]
        y= (i+1)/unitvector3*vector3[1]+Bpoint2[1]
        z= (i+1)/unitvector3*vector3[2]+Bpoint2[2]
        set.append([x,y,z])
    return set

#------------------------------------------------------------
# define a function
def outersphere(point):
# input the list of point on the object
    #point=[[0,1,0],[1,1,0],[2,1,1],[3,0,1],[4,-2,-1],[2,-2,-1]]
    i = len(point)
# seperate the x and y position from each point for using later
    xarray=[0]*i
    yarray=[0]*i
    zarray=[0]*i
    for i in range(len(point)):
    # assign
        xarray[i]=point[i][0]
        yarray[i]=point[i][1]
        zarray[i]=point[i][2]

    a = 0
# set points to store the two points with the maximum distance
    point1=[0]*3
    point2=[0]*3
# calculate the distance
    for i in range(len(xarray)):
        for j in range(len(xarray)):
            Dist = ((xarray[i]-xarray[j])**2+(yarray[i]-yarray[j])**2+(zarray[i]-zarray[j])**2)**0.5
            if Dist >= a:
                a = Dist
            #update two point with max distance
                point1[0]=xarray[i]
                point1[1]=yarray[i]
                point1[2]=zarray[i]
                point2[0]=xarray[j]
                point2[1]=yarray[j]
                point2[2]=zarray[j]

            j+=1
        j=0
        i+=1
# get the center part of sphere and the radius
    center=[0]*3
    center[0]=(point1[0]+point2[0])/2
    center[1]=(point1[1]+point2[1])/2
    center[2]=(point1[2]+point2[2])/2
    radius=a/2

# check if all the point are in the sphere
    b=radius
    for i in range(len(xarray)):
        Dist2=(((xarray[i]-center[0])**2+(yarray[i]-center[1])**2+(zarray[i]-center[2])**2)**0.5)/2
    # if point is out of sphere, then update
        if Dist2> b:
            b=Dist2
        i+=1
    radius=b
    print('center:',center,'radius:',radius)
    return center,radius,point1,point2

def define_area(point1,point2,point3):
    point1=np.asarray(point1)
    point2=np.asarray(point2)
    point3=np.asarray(point3)
    AB=np.asmatrix(point2-point1)
    AC=np.asmatrix(point3-point1)
    N=np.cross(AB,AC)
    Ax=N[0,0]
    By=N[0,1]
    Cz=N[0,2]
    D=-(Ax*point1[0]+By*point1[1]+Cz*point1[2])
    return Ax,By,Cz,D

def point2area_distance(point1,point2,point3,point4):
    Ax,By,Cz,D=define_area(point1,point2,point3)
    mod_d=Ax*point4[0]+By*point4[1]+Cz*point4[2]+D
    mod_area=np.sqrt(np.sum(np.square([Ax,By,Cz])))
    d=abs(mod_d)/mod_area
    return d

def volume(point1,point2,point3,centerpoint):

    p1=point2-point1
    p2=point3-point1
    p3=point3-point2
    l1=math.hypot(math.hypot(p1[0],p1[1]),p1[2])
    l2=math.hypot(math.hypot(p2[0],p2[1]),p2[2])
    l3=math.hypot(math.hypot(p3[0],p3[1]),p3[2])
    d=(l1+l2+l3)/2
    area=math.sqrt(d*(d-l1)*(d-l2)*(d-l3))
    d1=point2area_distance(point1,point2,point3,centerpoint)
    #print("point to surface distance:", str(d1))
    volume=area*d1/3
    #print("area:",area)
    #print("volume:", volume)
    return volume

#def pairs(volume,Rradius,distance):
def pairs(v,R,d):
    rmin=((v-4/3*math.pi*R**3)/(4*math.pi/3))**(1/3)
    rmax=(v/(4/3*math.pi))**(1/3)
    i=0
    ra=rmin
    rb=rmax
    while i<60:
        rmid=(ra+rb)/2
        x = -v*R**3*rmid**3+math.pi*(1/3*R**6*rmid**3+1/3*R**3*rmid**6+7/3*R**5*rmid**2*((R**2+d**2-rmid**2)/2)-5/3*R**4*rmid*((R**2+d**2-rmid**2)/2)**2+1/3*R**3*((R**2+d**2-rmid**2)/2)**3+7/3*rmid**5*R**2*((rmid**2+d**2-R**2)/2)-5/3*rmid**4*R*((rmid**2+d**2-R**2)/2)**2+1/3*rmid**3*((rmid**2+d**2-R**2)/2)**3)
        #x = -v+math.pi*(1/3*R**3+1/3*rmid**3-((-7/3*R**2*(R**2+d**2-rmid**2)/2/rmid)+(5/3*R*((R**2+d**2-rmid**2)/2/rmid)**2)-(1/3*((R**2+d**2-rmid**2)/2/rmid)**3))-((-7/3*rmid**2*((rmid**2+d**2-R**2)/2/R))+(5/3*rmid*((rmid**2+d**2-R**2)/2/R)**2)-(1/3*((rmid**2+d**2-R**2)/2/R)**3)))
        if x>0:
            ra=ra
            rb=rmid
        if x<0:
            ra=rmid
            rb=rb
        i+=1
    return rmid
# for comparing the variance of different pairs, and choose the best pair
def split(objectvolume,point1,point2,set):
    Vari=999999999
    for i in range(49):
        unit = (point2-point1)/100
        Acenter=point1+(i+1)*unit
        R=(i+1)*(unit[0]**2+unit[1]**2+unit[2]**2)**0.5
        #dmax=((objectvolume-4/3*math.pi*R**3)/(4/3*math.pi))**(1/3)+R
        for j in range(100):
            #drange=dmax+R
            drange=R/(i+1)*100-R
            distance=drange/100*(j+1)
            Bcenter=[Acenter[0]+distance/((unit[0]**2+unit[1]**2+unit[2]**2)**0.5)*unit[0],
                     Acenter[1]+distance/((unit[0]**2+unit[1]**2+unit[2]**2)**0.5)*unit[1],
                     Acenter[2]+distance/((unit[0]**2+unit[1]**2+unit[2]**2)**0.5)*unit[2]]
            Bcenter=np.array(Bcenter)
            r=pairs(objectvolume,R,distance)
            
            if R>= distance+r:
                continue
            if r>=distance+R:
                continue
            #
            vari=0
            for k in range(len(set)):
                A=set[k]-Acenter
                dist1=(A[0]**2+A[1]**2+A[2]**2)**(0.5)
                variance1=dist1-R
                B=set[k]-Bcenter
                dist2=(B[0]**2+B[1]**2+B[2]**2)**(0.5)
                variance2=dist2-r
                if(abs(variance1) <= abs(variance2)):
                    variance=abs(variance1)
                else:
                    variance=abs(variance2)
                vari+=variance
            if abs(vari)<abs(Vari):
                spherepoint1=Acenter
                spherepoint2=Bcenter
                Radius=R
                radius=r
                Vari=vari
    return spherepoint1,spherepoint2,Radius,radius,Vari
#plot the sphere we create
def circle(r,origX,origY,origZ):
    "This function is used to draw the sphere. The input is only radius r and the origial position on x,y,z"
    #fig=plt.figure()
    #ax = fig.add_subplot(111,projection='3d')
    u=np.linspace(0,2*np.pi,100)
    v=np.linspace(0,np.pi,100)
    x=r*np.outer(np.cos(u),np.sin(v))+origX
    y=r*np.outer(np.sin(u),np.sin(v))+origY
    z=r*np.outer(np.ones(np.size(u)),np.cos(v))+origZ
    #ax.plot_surface(x,y,z,color='g')
    return x,y,z

#------------------------------------------------------------
#read the object file.
'''objFilePath = '/home/xin/Downloads/my_temp_part7(1).obj'
with open(objFilePath) as file:
    points = []
    face = []
    while 1:
        line = file.readline()
        if not line:
            break
        strs = line.split(" ")
        if strs[0] == "v":
            points.append((float(strs[1]),float(strs[2]),float(strs[3])))
        if strs[0] == "f":
            face.append((int(strs[1]),int(strs[2]),int(strs[3])))
#I get points set and surfaces set
points = np.array(points)
face = np.array(face)
#get the centerpoint of outersphere, and the two most distal points
centerpoint,radiusoutersphere,point1,point2=outersphere(points)
centerpoint=np.array(centerpoint)
point1=np.array(point1)
point2=np.array(point2)

print('pointa:',point1,'pointb:',point2)

Volume=0
print('radiusoutersphere:',radiusoutersphere)
#get the volume of object
for i in range(len(face)):
    delta=volume(points[face[i][0]-1],points[face[i][1]-1],points[face[i][2]-1],centerpoint)
    Volume+=delta
print('Volume:',Volume)
# interplate the point on each side for the preparation of computing variance
a=[]
for i in range(len(face)):
    set =surface(points[face[i][0]-1],points[face[i][1]-1],points[face[i][2]-1])
    a=a+set
# split the object with the same volume, and compare the variance to get the best one.
#spherepoint1,spherepoint2,Radius,radius,Vari=split(Volume,point1,point2,set)
spherepoint1,spherepoint2,Radius,radius,Vari=split(Volume,point1,point2,points)
print('point1:',spherepoint1,'point2:',spherepoint2)
print('Radius:',Radius,'radius:',radius,'Vari:',Vari)
fig=plt.figure()
ax = plt.axes(projection='3d')
x1,y1,z1=circle(Radius,spherepoint1[0],spherepoint1[1],spherepoint1[2])
ax.plot_surface(x1,y1,z1,color='g')
x2,y2,z2=circle(radius,spherepoint2[0],spherepoint2[1],spherepoint2[2])
ax.plot_surface(x2,y2,z2,color='y')
plt.show()'''