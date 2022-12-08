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
    cosb=(r**2+d**2-R**2)/(2*r*d)
    if cosa>0 and cosb>0:
        vector1=center2-center1
        dist=(vector1[0]**2+vector1[1]**2+vector1[2]**2)**0.5
        BalllackingH1=R*(1-cosa)
        V1=np.pi*BalllackingH1**2*(R-BalllackingH1/3)
        BalllackingH2=r-(dist-BalllackingH1)
        V2=np.pi*BalllackingH2**2*(r-BalllackingH2/3)
        vector=vector1*(V1/(V1+V2))
        #vector=vector1*(R*cosa/(dist))
        origin=center1+vector
        normal=vector1
    if cosa<0:
        vector1=center2-center1
        dist=(vector1[0]**2+vector1[1]**2+vector1[2]**2)**0.5
        LackingpartH1=R+R*cosa
        V1=4/3*np.pi*R**3-(np.pi*LackingpartH1**2*(R-LackingpartH1/3))
        BalllackingH2=r-(dist+R*(-cosa))
        V2=np.pi*BalllackingH2**2*(r-BalllackingH2/3)
        vector=vector1*(V1/(V1+V2))
        origin=center1+vector
        normal=vector1
    if cosb<0:
        vector1=center1-center2
        dist=(vector1[0]**2+vector1[1]**2+vector1[2]**2)**0.5
        LackingpartH1=r+r*cosb
        V1=4/3*np.pi*r**3-(np.pi*LackingpartH1**2*(r-LackingpartH1/3))
        BalllackingH2=R-(dist+r*(-cosb))
        V2=np.pi*BalllackingH2**2*(R-BalllackingH2/3)
        vector=vector1*(V2/(V1+V2))
        origin=center2+vector
        normal=vector1

    return origin, normal

R=57.052399997386466
r=59.434789271448935

a=[-3.91178225,6.03496759,-67.44251252]
b=[-7.9997777,1.95607724,47.01280199]
d=((a[0]-b[0])**2+(a[1]-b[1])**2+(a[2]-b[2])**2)**0.5

if d>R+r:
    a=np.array(a)
    b=np.array(b)
    origin=(a-b)*(r/d)+b
    normal=a-b
if d<=R+r:
    origin, normal=Newsurf(R,r,d,a,b)
normal=np.array(normal)
normal=normal/((normal[0]**2+normal[1]**2+normal[2]**2)**0.5)

print('origin:',origin,'normal:', normal)