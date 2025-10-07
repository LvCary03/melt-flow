
import numpy as np
import matplotlib.pyplot as plt
import math

#import sys
#np.set_printoptions(threshold=sys.maxsize)

#Set vars
n = 100 #   #number of grid points
l = 10 #m   #total length of tube
cfl = .001 #s #time step
omega = (cfl + 1/cfl)/2

tf = .01 #s   #time to run

x0 = 5 #m
gamma = 1.4
dx = l / (n-1)

#Set initial values
P_left = 1 #Pa
P_right = 0.1 #Pa
roe_left = 1 #kg/m^3
roe_right = 0.125 #kg/m^3
u_left = 0 #m/s
u_right = 0 #m/s


ghost_cells = 1
list_length = n + 2*ghost_cells


W = np.zeros((3, list_length))
U = np.zeros((4, list_length))


#m = roe * u    momentum is density times velocity

#Set up the initial matricies
x = np.linspace(0, l, list_length)
for i in range(len(x)):
    if (x[i] < x0):
        U[0] = P_left
        U[1] = roe_left
        U[2] = u_left
    else:
        U[0] = P_right
        U[1] = roe_right
        U[2] = u_right


def primsToCons():
    for i in range(list_length):
        W[0, i] = U[1, i]
        W[1, i] = U[1, i] * U[2, i]
        #Equation to use for internal energy e = P / ((gamma - 1) * roe)
        W[2, i] = U[0, i] / ((gamma - 1) * U[1, i])
    return W


def consToPrims(W, gamma, list_length, No_blowup):
    for l in range(1, list_length-1): ###Should I include the boarder cases ???
        
        if (math.isnan(W[0, l])):
            print("code broke here roe")
            No_blowup = False
            break
        U[1, l] = W[0, l]
        
        if (math.isnan(W[1, l] / U[1, l])):
            print("code broke here u")
            No_blowup = False
            break
        U[2, l] = W[1, l] / U[1, l]
        
        #E = roe * (e + 1/2 u^2)
        # where e = P / ((gamma - 1) * roe)
        #U[1, l] * ((U[0, l] / ((gamma -1) * U[1, l])) + .5 * (U[2, l]**2))
    
        if (math.isnan(U[1, l] * ((U[0, l] / ((gamma -1) * U[1, l])) + .5 * (U[2, l]**2)))):
            print("code broke here E")
            No_blowup = False
            break
        U[3, l] = U[1, l] * ((U[0, l] / ((gamma -1) * U[1, l])) + .5 * (U[2, l]**2))
        
        if (math.isnan( W[2, l] * (gamma - 1) * U[1, l] )):
            print("code broke here P ")
            No_blowup = False
            break
        U[0, l] = W[2, l] * (gamma - 1) * U[1, l]
        return U, No_blowup



# sigma = max(abs(u)+c)*dt/dx
# dt = (sigma * dx) / (max(abs(u) + c))

# F(U) = [  m 
#          (m^2/roe) + P
#          (m/roe)(e + P)  ]

# alpha = omega * (dt/dx)*(u + c)


def calculate_alpha(list_length, cfl, dx, gamma):
    uc_max = 0
    alpha = np.zeros(list_length)
    for i in range(list_length):
        c = np.sqrt(gamma * U[0, i] / U[1, i])
        if ((abs(U[2, i]) + c ) > uc_max):
            uc_max = (abs(U[2, i] + c))

        dt = (cfl * dx) / (uc_max)
        #print("Time step = ", dt)

        alpha[i] = omega * (dt/dx)*(U[2, i] + c)
        return alpha, dt

def calculate_F(U, gamma, list_length):
    F = np.zeros((3, list_length))
    
    for j in range(list_length):
        #m = roe[j] * u[j]
        #e = (P[j]/(gamma - 1)) + (u[j])**2 ###Most definetely the wrong equaton
        #F[0, j] = m
        F[0, j] = U[1, j] * U[2, j]
        F[1, j] = (((U[1, j] * U[2, j])**2) / U[1, j]) + U[0, j]
        F[2, j] = (U[1, j] * U[2, j] / U[1, j]) * ((U[0, j]/((gamma - 1) * U[1, j])) + U[0, j])
        return F

No_blowup = True
total_time = tf
current_time = 0
while((current_time < total_time) and (No_blowup)):
#for step in range(1):


    #W = primsToCons(U, list_length)
    W = primsToCons()

    alpha, dt = calculate_alpha(list_length, cfl, dx, gamma)

    F = calculate_F(U, gamma, list_length)

    Wn = W.copy()

    for k in range(1, (len(alpha)-1)):
        Wn[:, k] = (W[:, k] -
                   (dt/(2*dx))*(F[:, k+1] - F[:, k-1]) +
                   (1/4)*(((alpha[k+1] + alpha[k])*(W[:, k+1] - W[:, k])) - 
                          (alpha[k] - alpha[k-1])*(W[:, k] - W[:, k-1])))
        
    W = Wn.copy()

    U, No_blowup = consToPrims(W, gamma, list_length, No_blowup)
  

    current_time += dt




#Print outputs
print("current_time = ", current_time)
print("Pressure matrix =", U[0])
if (No_blowup):
    print("The system did not blowup")
else:
    print("The system blew up")
#print("Density matrix =", roe)
#print("\nConserved variables\nDensity=", W[0],'\nMomentum=', W[1], '\nEnergy=', W[2])
#print("F = ", F)


#Give the plots for all the variables
fig, axs = plt.subplots(2, 2)

# Plot on the first subplot (top-left)
axs[0, 0].plot(W[1])
axs[0, 0].set_title('Density')

# Plot on the second subplot (top-right)
axs[0, 1].plot(U[2], color='orange')
axs[0, 1].set_title('Velocity')

# Plot on the third subplot (bottom-left)
axs[1, 0].plot(U[3], color='green')
axs[1, 0].set_title('Energy')

# Plot on the fourth subplot (bottom-right)
axs[1, 1].plot(U[0], color='red')
axs[1, 1].set_title('Pressure')

plt.show()


