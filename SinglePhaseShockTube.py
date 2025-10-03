
import numpy as np
import matplotlib.pyplot as plt
import math

#import sys
#np.set_printoptions(threshold=sys.maxsize)

#Set vars
n = 100 #   #number of grid points
l = 10 #m   #total length of tube
cfl = .1 #s #time step
omega = (cfl + 1/cfl)/2
tf = .1 #s   #time to run
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



'''
P = np.zeros((n+2))
u = np.zeros((n+2))
roe = np.zeros((n+2))
W = np.zeros((3, n+2))
E = np.zeros((n+2))'''

W = np.zeros((3, n+2))
U = np.zeros((4, n+2))


#m = roe * u    momentum is density times velocity
#E = p/(gamma-1) + u^2  ##Not sure what this equation is doing???

#Set up the initial matricies
x = np.linspace(0, l, n+2)
for i in range(len(x)):
    if (x[i] < x0):
        #P[i] = P_left
        U[0] = P_left
        #roe[i] = roe_left
        U[1] = roe_left
        #u[i] = u_left
        U[2] = u_left
    else:
        #P[i] = P_right
        #roe[i] = roe_right
        #u[i] = u_right
        U[0] = P_right
        U[1] = roe_right
        U[2] = u_right
    #W[0, i] = roe[i]
    #W[1, i] = roe[i] * u[i]
    #W[2, i] = (P[i] / gamma -1) + u[i]**2    #Potential break here. Wrong formula ???  #Most definetely the wrong equation 
    W[0, i] = U[1, i]
    W[1, i] = U[1, i] * U[2, i]
    #Equation to use for internal energy e = P / ((gamma - 1) * roe)
    W[2, i] = U[0, i] / ((gamma - 1) * U[1, i])


#roe^n+1_i = roe^n_i + ( roe^n_i+1 - roe^n_i-1 ) / dx**2

# sigma = max(abs(u)+c)*dt/dx
# dt = (sigma * dx) / (max(abs(u) + c))

# F(U) = [  m 
#          (m^2/roe) + P
#          (m/roe)(e + P)  ]

# alpha = omega * (dt/dx)*(u + c)

No_blowup = True
total_time = tf
current_time = 0
while((current_time < total_time) and (No_blowup)):
    uc_max = 0
    F = np.zeros((3, n+2))
    alpha = np.zeros(n+2)
    for j in range(len(W[1])):
        #m = roe[j] * u[j]
        #e = (P[j]/(gamma - 1)) + (u[j])**2 ###Most definetely the wrong equaton
        #F[0, j] = m
        F[0, j] = U[1, j] * U[2, j]
        #F[1, j] = ((m**2)/roe[j]) + P[j]
        F[1, j] = (((U[1, j] * U[2, j])**2) / U[1, j]) + U[0, j]
        #F[2, j] = (m/roe[j])*(e + P[j])
        # e = P/((gamma - 1) * roe)
        F[2, j] = (U[1, j] * U[2, j] / U[1, j]) * ((U[0, j]/((gamma - 1) * U[1, j])) + U[0, j])

        #c = np.sqrt(gamma*P[j]/roe[j])
        c = np.sqrt(gamma * U[0, j] / U[1, j])
        #if ((abs(u[j]) + c ) > uc_max):
        if ((abs(U[2, j]) + c ) > uc_max):
            #uc_max = (abs(u[j] + c))
            uc_max = (abs(U[2, j] + c))

        dt = (cfl * dx) / (uc_max)
        #print("Time step = ", dt)

        #alpha[j] = omega * (dt/dx)*(u[j] + c)
        alpha[j] = omega * (dt/dx)*(U[2, j] + c)

    #Fn = F.copy()
    Wn = W.copy()
    #alphan = alpha.copy()


    for k in range(1, (len(alpha)-1)):      #Potential error, should all matricies in the equation be 'new'
        Wn[:, k] = (W[:, k] -
                   (dt/(2*dx))*(F[:, k+1] - F[:, k-1]) +
                   (1/4)*(((alpha[k+1] + alpha[k])*(W[:, k+1] - W[:, k])) - 
                          (alpha[k] - alpha[k-1])*(W[:, k] - W[:, k-1])))

    #print(W)

    # E = e + V^2 / 2
    for l in range(1, (len(alpha)-1)):
        
        if (math.isnan(W[0, l])):
            print("code broke here roe")
            No_blowup = False
            break
        #roe[l] = W[0, l]            #Potential break. We just assign our old density to be out new density
        U[1, l] = W[0, l]
        
        if (math.isnan(W[1, l] / U[1, l])):
            print("code broke here u")
            No_blowup = False
            break
        #u[l] = W[1, l] / roe[l]
        U[2, l] = W[1, l] / U[1, l] #Here I use the just calculated value of roe to calculat U[2, l]
        
        #E = roe * (e + 1/2 u^2)
        # where e = P / ((gamma - 1) * roe)

        #U[1, l] * ((U[0, l] / ((gamma -1) * U[1, l])) + .5 * (U[2, l]**2))
    
        if (math.isnan(U[1, l] * ((U[0, l] / ((gamma -1) * U[1, l])) + .5 * (U[2, l]**2)))):
            print("code broke here E")
            No_blowup = False
            break
        #E[l] = W[2, l] + (u[l]**2)/2
        U[3, l] = U[1, l] * ((U[0, l] / ((gamma -1) * U[1, l])) + .5 * (U[2, l]**2))
        
        if (math.isnan( W[2, l] * (gamma - 1) * U[1, l] )):
            print("code broke here P ")
            No_blowup = False
            break
        #U[0, l] = (W[2, l] - U[2, l]**2) * gamma
        U[0, l] = W[2, l] * (gamma - 1) * U[1, l]

    current_time += dt




#Print outputs
print("current_time = ", current_time)
print("Pressure matrix =", U[0])
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