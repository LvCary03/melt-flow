
import numpy as np
import matplotlib.pyplot as plt


#Set vars
n = 1000 #   #number of grid points
l = 10 #m   #total length of tube
cfl = .5 #s #time step
tf = 4 #s   #time to run
x0 = 5 #m
gamma = 1.4
dx = l / (n-1)
#print(dx)


#Set initial values
P_left = 1 #Pa
P_right = 0.1 #Pa
roe_left = 1 #kg/m^3
roe_right = 0.125 #kg/m^3
u_left = 0 #m/s
u_right = 0 #m/s

P = np.zeros((n+2))
u = np.zeros((n+2))
roe = np.zeros((n+2))
W = np.zeros((3, n+2))
E = np.zeros((n+2))

#P[0:n/2 +1] = p_left
#Set initials

#roe * u
#E = p/(gamma-1) + u^2

x = np.linspace(0, l, n+2)
for i in range(len(x)):
    if (x[i] < x0):
        P[i] = P_left
        roe[i] = roe_left
        u[i] = u_left
    else:
        P[i] = P_right
        roe[i] = roe_right
        u[i] = u_right
    W[0, i] = roe[i]
    W[1, i] = roe[i] * u[i]
    W[2, i] = P[i] / gamma + u[i]**2


#roe^n+1_i = roe^n_i + ( roe^n_i+1 - roe^n_i-1 ) / dx**2

# sigma = max(abs(u)+c)*dt/dx
# dt = (sigma * dx) / (max(abs(u) + c))


# F(U) = [  m 
#          (m^2/roe) + P
#          (m/roe)(e + P)  ]

# alpha = omega * (dt/dx)*(u + c)

run = True
total_time = tf
current_time = 0
while(current_time < total_time):
    omega = .75
    uc_max = 0
    F = np.zeros((3, n+2))
    alpha = np.zeros(n+2)
    for j in range(len(P)):
        m = roe[j] * u[j]
        e = P[j]/(gamma - 1) + (u[j])**2
        #print("e =", e)
        #print((m/roe)*(e + P[j]))
        F[0, j] = m
        F[1, j] = ((m**2)/roe[j]) + P[j]
        F[2, j] = (m/roe[j])*(e + P[j])


        c = np.sqrt(gamma*P[j]/roe[j])
        if ((abs(u[j]) + c ) > uc_max):
            uc_max = (abs(u[j] + c))

        dt = (cfl * dx) / (uc_max)
        #print("Time step = ", dt)

        alpha[j] = omega * (dt/dx)*(u[j] + c)




    Fn = F.copy()
    Wn = W.copy()

    alphan = alpha.copy()


    for k in range(1, (len(P)-1)):
        W[:, k] = (Wn[:, k] -
                   (dt/(2*dx))*(F[:, k+1] - F[:, k-1]) +
                   (1/4)*(((alpha[k+1] + alpha[k])*(Wn[:, k+1] - Wn[:, k])) - 
                          (alpha[k] - alpha[k-1])*(Wn[:, k] - Wn[:, k-1])))

    #print(W)

    # E = e + V^2 / 2
    for l in range(1, (len(P)-1)):
        E[l] = W[2, l] + (u[l]**2)/2
        u[l] = W[1, l] / roe[l]



    #run = False
    current_time += dt


print("current_time = ", current_time)
#print("Pressure matrix =", P)
#print("Density matrix =", roe)
#print("\nConserved variables\nDensity=", W[0],'\nMomentum=', W[1], '\nEnergy=', W[2])
#print("F = ", F)


fig, axs = plt.subplots(2, 2)

# Plot on the first subplot (top-left)
axs[0, 0].plot(W[0])
axs[0, 0].set_title('Density')

# Plot on the second subplot (top-right)
axs[0, 1].plot(u, color='orange')
axs[0, 1].set_title('Velocity')

# Plot on the third subplot (bottom-left)
axs[1, 0].plot(E, color='green')
axs[1, 0].set_title('Energy')

# Plot on the fourth subplot (bottom-right)
axs[1, 1].plot(P, color='red')
axs[1, 1].set_title('Pressure')






#plt.plot(W[0], label = 'density')
#plt.plot(W[1], label = 'momentum')
#plt.plot(W[2], label = 'energy')
plt.show()