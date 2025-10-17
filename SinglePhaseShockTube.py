
import numpy as np
import matplotlib.pyplot as plt
import math

#import sys
#np.set_printoptions(threshold=sys.maxsize)

#Set vars
n = 100 #   #number of grid points
l = 10 #m   #total length of tube
cfl = .05 #s #time step
omega = (cfl + 1/cfl)/2

tf = .3 #s   #time to run

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
        U[0, i] = P_left
        U[1, i] = roe_left
        U[2, i] = u_left
        #U[3, i] = U[1, i] * ((U[0, i] / ((gamma -1) * U[1, i])) + .5 * (U[2, i]**2))
        U[3, i] = U[0, i] / ((gamma - 1) * U[1, i])
    else:
        U[0, i] = P_right
        U[1, i] = roe_right
        U[2, i] = u_right
        #U[3, i] = U[1, i] * ((U[0, i] / ((gamma -1) * U[1, i])) + .5 * (U[2, i]**2))
        U[3, i] = U[0, i] / ((gamma - 1) * U[1, i])


def primsToCons(U, list_length):
    for i in range(list_length):
        W[0, i] = U[1, i]
        W[1, i] = U[1, i] * U[2, i]
        #Equation to use for internal energy e = P / ((gamma - 1) * roe)
        #W[2, i] = U[0, i] / ((gamma - 1) * U[1, i])
        W[2, i] = U[1, i] * ((U[0, i] / ((gamma -1) * U[1, i])) + .5 * (U[2, i]**2))
    return W


def consToPrims(W, gamma, list_length, No_blowup):
    # Convert conserved W -> primitive U.
    # Here W[0] = rho, W[1] = momentum, W[2] = total energy per volume (rho*(e + 0.5 u^2))
    for idx in range(1, list_length-1):
        rho = W[0, idx]
        mom = W[1, idx]
        E   = W[2, idx]

        if rho <= 0 or math.isnan(rho):
            print("code broke here: rho nonpositive or NaN at", idx, rho)
            No_blowup = False
            break
        U[1, idx] = rho

        u = mom / rho        # THIS is the correct velocity calculation
        if math.isnan(u) or abs(u) > 1e12:
            print("code broke here: u invalid at", idx, u)
            No_blowup = False
            break
        U[2, idx] = u

        # internal energy per mass: e = (E/rho) - 0.5*u^2
        e = (E / rho) - 0.5 * u**2
        if math.isnan(e):
            print("code broke here: e invalid at", idx, e)
            No_blowup = False
            break

        # pressure: P = (gamma - 1) * rho * e
        P = (gamma - 1) * rho * e
        U[0, idx] = P

        # store e per mass optionally in U[3], or store total energy per mass: e + 0.5 u^2
        U[3, idx] = e #+ 0.5 * u**2 # I just want this one as the same internal energy
    return U, No_blowup

def mainLoop(Wn, W, alpha, F):
    for k in range(1, (len(alpha)-1)):
        Wn[:, k] = (W[:, k] -
                   (dt/(2*dx))*(F[:, k+1] - F[:, k-1]) +
                   (1/4)*(((alpha[k+1] + alpha[k])*(W[:, k+1] - W[:, k])) - 
                          (alpha[k] - alpha[k-1])*(W[:, k] - W[:, k-1])))
        
    # I think here is where I set boundary conditions
    Wn[:, 0] = Wn[:, 1].copy()
    Wn[:, -1] = Wn[:, -2].copy()
        
    return Wn, W, alpha, F

# sigma = max(abs(u)+c)*dt/dx
# dt = (sigma * dx) / (max(abs(u) + c))

# F(U) = [  m 
#          (m^2/roe) + P
#          (m/roe)(e + P)  ]

# alpha = omega * (dt/dx)*(u + c)


def calculate_alpha(U, list_length, cfl, dx, gamma):
    uc_max = 0
    alpha = np.zeros(list_length)
    for i in range(list_length):
        c = np.sqrt(gamma * U[0, i] / U[1, i])
        if ((abs(U[2, i]) + c ) > uc_max):
            uc_max = abs(U[2, i]) + c

    dt = (cfl * dx) / (uc_max)
        #print("Time step = ", dt)
    for j in range(list_length):
        c = np.sqrt(gamma * U[0, j] / U[1, j])
        alpha[j] = omega * (dt/dx)*(U[2, j] + c)
    return alpha, dt

def calculate_F(U, gamma, list_length):
    F = np.zeros((3, list_length))
    
    for j in range(list_length):
        #m = roe[j] * u[j]
        #F[0, j] = m
        F[0, j] = U[1, j] * U[2, j]
        F[1, j] = (((U[1, j] * U[2, j])**2) / U[1, j]) + U[0, j]
        #F[2, j] = (U[1, j] * U[2, j] / U[1, j]) * ((U[0, j]/((gamma - 1) * U[1, j])) + U[0, j])
        F[2, j] = (U[1, j] * U[2, j] / U[1, j]) * ((U[1, j] * ((U[0, j] / ((gamma -1) * U[1, j])) + .5 * (U[2, j]**2))) + U[0, j])
    return F


No_blowup = True
total_time = tf
current_time = 0
while((current_time < total_time) and (No_blowup)):
#for step in range(2):

    print(U[2])

    W = primsToCons(U, list_length)


    alpha, dt = calculate_alpha(U, list_length, cfl, dx, gamma)

    F = calculate_F(U, gamma, list_length)

    Wn = W.copy()


    #for k in range(1, (len(alpha)-1)):
     #   Wn[:, k] = (W[:, k] -
      #             (dt/(2*dx))*(F[:, k+1] - F[:, k-1]) +
       #            (1/4)*(((alpha[k+1] + alpha[k])*(W[:, k+1] - W[:, k])) - 
        #                  (alpha[k] - alpha[k-1])*(W[:, k] - W[:, k-1])))
    
    Wn, W, alpha, F = mainLoop(Wn, W, alpha, F)
        
    W = Wn.copy()

    U, No_blowup = consToPrims(W, gamma, list_length, No_blowup)


    print(U[2])
    #print("F = ", F)
  

    current_time += dt



#Print outputs
print("current_time = ", current_time)
#print("Pressure matrix =", U[0])
if (No_blowup):
    print("The system did not blowup")
else:
    print("The system blew up")
#print("\nConserved variables\nDensity=", W[0],'\nMomentum=', W[1], '\nEnergy=', W[2])
#print("F = ", F)


def plotElementary():
    #Give the plots for all the variables
    fig, axs = plt.subplots(2, 2)

    # Plot on the first subplot (top-left)
    axs[0, 0].plot(U[1])
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


def plotConserved():
    fig2, axs2 = plt.subplots(2, 2)

    # Plot on the first subplot (top-left)
    axs2[0, 0].plot(W[1])
    axs2[0, 0].set_title('W 1')

    # Plot on the second subplot (top-right)
    axs2[0, 1].plot(W[2], color='orange')
    axs2[0, 1].set_title('W 2')

    # Plot on the third subplot (bottom-left)
    #axs2[1, 0].plot(U[3], color='green')
    #axs2[1, 0].set_title('Energy')

    # Plot on the fourth subplot (bottom-right)
    axs2[1, 1].plot(W[0], color='red')
    axs2[1, 1].set_title('W 0')


plotElementary()
#plotConserved()

plt.show()

