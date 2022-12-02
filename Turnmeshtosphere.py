import numpy as np

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

#sphere function

"""points=[[1,1,1],
        [2,1,2],
        [3,-1,-1],
        [1,4,4],
        [6,6,6]]
center, radius, point1, point2=outersphere(points)
print('center:',center,'radius:',radius)"""

"""points=np.array(points)
center=np.array(center)
unit=radius/10
vari=0
Vari=999999999999999999
for i in range(10):
    r=unit*(i+1)
    print('r:',r)
    for j in range(len(points)):
        dis=points[j]-center
        vari1=(dis[0]**2+dis[1]**2+dis[2]**2)**0.5-r
        vari=vari+vari1
    if abs(vari) < abs(Vari):
        Vari=vari
        Rad=r
    print('Vari:',vari)
    vari=0
    print('Rad:', Rad)"""

def Turnmeshtosphere(points, center, radius):
    points=np.array(points)
    center=np.array(center)
    unit=radius/10
    vari=0
    Vari=999999999999999999
    for i in range(10):
        r=unit*(i+1)
        print('r:',r)
        for j in range(len(points)):
            dis=points[j]-center
            vari1=(dis[0]**2+dis[1]**2+dis[2]**2)**0.5-r
            vari=vari+vari1
        if abs(vari) < abs(Vari):
            Vari=vari
            Rad=r
        print('Vari:',vari)
        vari=0
        print('Rad:', Rad)
    return Rad

points=[[1,1,1],
        [2,1,2],
        [3,-1,-1],
        [1,4,4],
        [6,6,6]]
center, radius, point1, point2=outersphere(points)
print('center:',center,'radius:',radius)
Rad=Turnmeshtosphere(points,center,radius)
print('Rad:',Rad)