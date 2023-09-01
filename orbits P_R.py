import matplotlib.pyplot as plt

d = 10
G = 6.6743*10**(-11)
Msun = 2*10**30
Lsun = 3.828*10*26
Rsun = 7*10**8
pi = 3.141592
c = 3*10**8

def fx(x, y):
    return -G*Msun*m*x * (x*x+y*y)**(-1.5)
def fy(x, y):
    return -G*Msun*m*y * (x*x+y*y)**(-1.5)

def PR_x(x, y, R, vx, vy):
    r = (x*x+y*y)**0.5
    v = (vx*vx+vy*vy)**0.5
    return Lsun/4 * (R/r)**2 * v/(c*c)  *  x/r

def PR_y(x, y, R, vx, vy):
    r = (x*x+y*y)**0.5
    v = (vx*vx+vy*vy)**0.5
    return -Lsun/4 * (R/r)**2 * v/(c*c)  *  y/r


rx = 0
ry = 3*10**9
vx = (1*G*Msun*(ry**2+rx**2)**(-0.5))**0.5
vy = vx*0.5
T = 1600000

R = 3*10**7
A = 0
m = 10**-12

plotX=[]
plotY=[]
limsun = (Rsun*1.05)**2
for i in range(int(T/100/d)):
    if rx**2+ry**2 < limsun:
        break
    for ii in range(0, 99):
        vx += (fx(rx,ry) + PR_y(rx, ry, R, vx, vy))/m*d
        rx += vx * d
        plotX.append(rx)
    
        vy += (fy(rx,ry) + PR_x(rx, ry, R, vx, vy))/m*d
        ry += vy * d
        plotY.append(ry)
    
    

print(PR_x(rx, ry, R, vx, vy)/fx(rx,ry))
print((vx**2+vy**2)**0.5*(G*Msun*(ry**2+rx**2)**(-0.5))**(-0.5))

fig, ax = plt.subplots()
x0,x1 = ax.get_xlim()
y0,y1 = ax.get_ylim()

ax.set_aspect(abs(x1-x0)/abs(y1-y0))

circle = plt.Circle((0, 0), Rsun, color='orange')
ax.add_patch(circle)
plt.plot(plotX, plotY, color='black')
plt.scatter([plotX[len(plotX)-1]], [plotY[len(plotY)-1]], s=50, color='red')
plt.scatter([-plotX[len(plotX)-1]], [-plotY[len(plotY)-1]], s=0, color='white')
plt.scatter([plotY[len(plotY)-1]], [0], s=0, color='white')
plt.scatter([-plotY[len(plotY)-1]], [0], s=0, color='white')
plt.show()


'''
def E_x(x, y, R):
    return -0.25*x*R*R*Lsun * (x*x+y*y)**(-1.5)
def E_y(x, y, R):
    return -0.25*y*R*R*Lsun * (x*x+y*y)**(-1.5)
    '''
