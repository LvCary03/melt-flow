
import numpy as np
import matplotlib.pyplot as plt
import math
from scipy.optimize import fsolve

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


#U_ex_L = np.zeros((3, list_length/2))
#U_ex_R = np.zeros((3, list_length/2))

U_ex_L = [roe_left, u_left, P_left]
U_ex_R = [roe_right, u_right, P_right]


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

        #U_ex_L[0, i] = P_left/(RL*TL)
        #U_ex_L[1, i] = u_left
        #U_ex_L[2, i] = P_left
        
        #U_ex_R[0, i] = P_right/(RR*TR)
        #U_ex_R[1, i] = u_right
        #U_ex_R[2, i] = P_righ
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
    c_mat = np.zeros(list_length)
    print(c_mat)
    for j in range(list_length):
        c = np.sqrt(gamma * U[0, j] / U[1, j])
        c_mat[j] = c 
        alpha[j] = omega * (dt/dx)*(U[2, j] + c)
    return alpha, dt, c_mat

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
#for step in range(40*4):

    #print(U[2])

    W = primsToCons(U, list_length)


    alpha, dt, c_mat = calculate_alpha(U, list_length, cfl, dx, gamma)

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


    #print(U[2])
    #print("F = ", F)
  

    current_time += dt



def eqn(r_p, p4op1, gam_m, a_4, u_4, u_1, a_1, gam_p):
    return p4op1 - r_p * (1 + gam_m/(2*a_4) * (
        u_4 - u_1 - a_1/gamma * (r_p - 1) /
        np.sqrt(gam_p / (2*gamma) * (r_p - 1) + 1)
    )) ** (-2*gamma/gam_m)

def riemann_exact(gamma,x0,x,tf,U_L,U_R):
    rho_L = U_L[0]                          #% Density at driver end 
    u_L = U_L[1]                            #% Velocity at driver end 
    p_L = U_L[2]                            #% Pressure at driver end 
    rho_R = U_R[0]                          #% Density at driven end 
    u_R = U_R[1]                            #% Velocity at driven end 
    p_R = U_R[2]                            #% Pressure at driven end 
    x_i = min(x)  
    x_f = max(x)              #% Grid dimensions
    npt = len(x)               #          % # of grid points

    p_4 = p_L
    rho_4 = rho_L
    u_4 = u_L
    p_1 = p_R
    rho_1 = rho_R
    u_1 = u_R
    gam_m = gamma-1                          #% Specific heat ratio notation
    gam_p = gamma+1
    gam_r = gam_p/gam_m
    p4op1 = p_4/p_1                        #% Driver/driven end pressure ratio
    print(f"value of roe {rho_4}")
    a_4 = math.sqrt(gamma*p_4/rho_4)
    a_1 = math.sqrt(gamma*p_1/rho_1)



    p2op1 = fsolve(eqn, np.pi, args=(p4op1, gam_m, a_4, u_4, u_1, a_1, gam_p))[0] #Call the exact solution

    ### Resolve properties
    p_2 = p2op1*p_1
    u_2 = u_4 + 2*a_4/gam_m*(1-(p4op1**(-1)*p2op1)**(gam_m/(2*gamma)))                             #% Eqn 5.4
                 
    a_2 = a_1*math.sqrt(p2op1*(gam_r+p2op1)/(1+gam_r*p2op1))     #% Eqn 3.54
    u_s = u_1 + a_1*math.sqrt(gam_p/(2*gamma)*(p2op1-1)+1)        #% Eqn 3.56
    u_3 = u_2                                              #% Eqn 3.57
    p_3 = p_2                                              #% Eqn 3.58
    a_3 = a_4 + gam_m/2*(u_4-u_3)                          #% Eqn 5.2


    #Replaced all '&&' with 'and'
    x_4 = x0 + (u_4-a_4)*tf
    x_3 = x0 + (u_3-a_3)*tf
    x_2 = x0 + u_2*tf
    x_1 = x0 + u_s*tf

    u = np.zeros(npt)
    a = np.zeros(npt)
    p = np.zeros(npt)
    rho = np.zeros(npt)
    '''u = []
    a = []
    p = []
    rho = []
    for i in range(1, npt):
        u.append(0)
        a.append(0)
        p.append(0)
        rho.append(0)'''

    for i in range(1,npt):
        if (x[i] < x_4):
            u[i] = u_4
            a[i] = a_4
            p[i] = p_4         
            rho[i] = rho_4
        elif (x[i] >= x_4 and x[i] <= x_3):
            u[i] = 2/gam_p*((x[i]-x0)/tf + gam_m/2*u_4 + a_4)    #% Eqn 3.47
            a[i] = u[i] - (x[i]-x0)/tf                          #% Eqn 3.48
            p[i] = p_4*(a[i]/a_4) **(2*gamma/gam_m)                 #% Eqn 3.49
            rho[i] = gamma*p[i]/a[i] **2
        elif (x[i] >= x_3 and x[i] <= x_2):
            u[i] = u_3
            a[i] = a_3
            p[i] = p_3  
            rho[i] = gamma*p[i]/a[i] **2
        elif (x[i] >= x_2 and x[i] <= x_1):
            u[i] = u_2
            a[i] = a_2
            p[i] = p_2  
            rho[i] = gamma*p[i]/a[i] **2
        elif (x[i] > x_1):
            u[i] = u_1
            a[i] = a_1
            p[i] = p_1
            rho[i] = rho_1

    print("important new variables")
    print(u, a, p)
    plt.plot(u)
    #plt.figure()
    #plt.plot(a)
    plt.figure()
    plt.plot(p)

print("Working on exact solution")
#Call the function
riemann_exact(gamma,x0,x,tf,U_ex_L,U_ex_R)


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
    fig, axs = plt.subplots(2, 3)

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

    # Plot on the fourth subplot (bottom-right)
    #axs[0, 2].plot(c_mat, color='blue')
    #axs[0, 2].set_title('sound')


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

