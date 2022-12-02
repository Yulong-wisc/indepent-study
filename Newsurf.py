import numpy as np

"""R=10
r=8
d=15
a=[0,0,0]
b=[6*2**0.5,6*2**0.5,9]
a=np.array(a)
b=np.array(b)
cosa=(R**2+d**2-r**2)/(2*R*d)
vector1=b-a
print('cosa:',cosa)
dis=(vector1[0]**2+vector1[1]**2+vector1[2]**2)**0.5
vector=vector1*(R*cosa/(dis))
print('vector:',vector)
origin=a+vector
print('origin:',origin)"""

def Newsurf(R,r,d,center1,center2):
    center1=np.array(center1)
    center2=np.array(center2)
    cosa=(R**2+d**2-r**2)/(2*R*d)
    vector1=center2-center1
    dist=(vector1[0]**2+vector1[1]**2+vector1[2]**2)**0.5
    vector=vector1*(R*cosa/(dist))
    origin=center1+vector
    normal=vector1
    return origin, normal

R=1.3599742384857518
r=2.7759750362659563

a=[10.10190834,1.34563084,2.1240758 ]
b=[-3.55489543,1.21790984,-0.17204387]
d=((a[0]-b[0])**2+(a[1]-b[1])**2+(a[2]-b[2])**2)**0.5

if d>R+r:
    a=np.array(a)
    b=np.array(b)
    origin=(a-b)*(r/d)+b
    normal=a-b
if d<=R+r:
    origin, normal=Newsurf(R,r,d,a,b)
normal=np.array(normal)
#normal=normal/((normal[0]**2+normal[1]**2+normal[2]**2)**0.5)

print('origin:',origin,'normal:', normal)