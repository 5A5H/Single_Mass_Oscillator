# Single mass oscillator
# Newmark time integration
import matplotlib.pyplot as plt;

# Parameters
dt 		= 0.01;    # time incremenet
nstep = 2000;    # number of time steps

# Parameters:
k = 225; # spring stiffness
m = 1e0; # mass
d = 1.0; # damping parameter

# start values
# current position, velocity, acceleration 
x = 0.0; v = 0.0; a = 0.0;

# setup loading function
def f(t):
    f_ext = 0.0;
    if (t<5):
        f_ext = t;
    if (t>=5):
       f_ext = 0;
    return f_ext

# postprocessing field
time = [0.0];
position = [x];
velocity = [v];
acceleration = [a];
ext_force = [f(0)];  

# time integration loop
t = 0.0;
for step in range(nstep):
    # update time
    x_n, v_n, a_n = x, v, a;
    t = t + dt;
    if (step%(nstep/10)==0): print('Current time: %5.1f'%t);
    
    # predecessors
    f_ext = f(t);
    #print('f_ext',f_ext);
    f_int = m * a_n + d * v_n + k * x_n;
    f_NM  = m * ((4.0/dt)*v_n+2.0*a_n)+2*d*x_n;
    f_delta = (4.0/dt**2.0)*m+d*(2.0/dt)+k;
    
    u = (f_ext-f_int+f_NM)/f_delta;
    
    x = x + u;
    v = - v_n + (2.0/dt)*u;
    a = (4.0/dt**2.0) * (u-v_n*dt) - a_n;
        
    # postprocessing
    time.append(t);
    position.append(x);
    velocity.append(v);
    acceleration.append(a);
    ext_force.append(f(t));
    
print('Finish time integration!');

# Print results
fig, axs = plt.subplots(2, 2);
fig.set_size_inches(10,10);

axs[0, 0].plot(time, position ,color='black', linewidth=1.5);
axs[0, 0].set_title('Position', fontsize=10)	
axs[0, 0].set_xlabel('time', fontsize=10, color='black')
axs[0, 0].set_ylabel('position', fontsize=10, color='black')

axs[0, 1].plot(time, velocity ,color='green', linewidth=1.5);
axs[0, 1].set_title('Velocity', fontsize=10)	
axs[0, 1].set_xlabel('time', fontsize=10, color='black')
axs[0, 1].set_ylabel('velocity', fontsize=10, color='black')

axs[1, 0].plot(time, acceleration ,color='blue', linewidth=1.5);
axs[1, 0].set_title('Acceleration', fontsize=10)	
axs[1, 0].set_xlabel('time', fontsize=10, color='black')
axs[1, 0].set_ylabel('acceleration', fontsize=10, color='black')

axs[1, 1].plot(time, ext_force ,color='red', linewidth=1.5);
axs[1, 1].set_title('External Force', fontsize=10)	
axs[1, 1].set_xlabel('time', fontsize=10, color='black')
axs[1, 1].set_ylabel('load', fontsize=10, color='black')

fig.tight_layout();
plt.show();
