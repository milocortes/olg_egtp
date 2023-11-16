import random
import numpy as np
from multiprocessing import Pool,cpu_count, Value, Array

# Para hacer el muestreo por Latin Hypecube
from scipy.stats.qmc import LatinHypercube,scale
import math
import ctypes

'''
------------------------------------------
                    PSO
        Particle Swarm Optimization

-------------------------------------------
## Implemented as a minimization algorithm

# Inputs:
    * f_cost        - function to be minimized
    * pop_size      - number of individuals in the population
    * max_iters     - maximum number of optimization iterations
    * lb            - lower bounds
    * ub            - upper bounds
    * α             - cognitive scaling parameter
    * β             - social scaling parameter
    * w             - velocity inertia 
    * w_min         - minimum value for the velocity inertia
    * w_max         - maximum value for the velocity inertia
# Output
    * best_theta    - best solution found
    * best_score    - history of best score
'''

    
def init_pso(gbest_val_arg, gbest_pos_arg, position_arg, velocity_arg, pbest_val_arg, 
         pbest_pos_arg,f_optim,α_arg,β_arg,w_arg,vMax_arg,vMin_arg,
         u_bounds_arg,l_bounds_arg):
    global gbest_val
    global gbest_pos
    global position
    global velocity
    global pbest_val
    global pbest_pos
    global f
    global α
    global β
    global w
    global vMax
    global vMin
    global u_bounds
    global l_bounds

    gbest_val = gbest_val_arg
    gbest_pos = gbest_pos_arg
    position = position_arg
    velocity = velocity_arg
    pbest_val = pbest_val_arg
    pbest_pos = pbest_pos_arg
    f = f_optim
    
    # Cognitive scaling parameter
    α = α_arg
    # Social scaling parameter
    β = β_arg
    
    # velocity inertia
    w = w_arg
    
    vMax = vMax_arg
    vMin = vMin_arg
    u_bounds = u_bounds_arg
    l_bounds = l_bounds_arg
    
def evalua_f_pso(i):    
    # Actualiza velocidad de la partícula
    ϵ1,ϵ2 = np.random.RandomState().uniform(), np.random.RandomState().uniform()
    with gbest_pos.get_lock():
        velocity[i] = w.value*velocity[i] + α*ϵ1*(pbest_pos[i] -  position[i]) + β*ϵ2*(np.array(gbest_pos[:]) - position[i])

            
    # Ajusta velocidad de la partícula
    index_vMax = np.where(velocity[i] > vMax)
    index_vMin = np.where(velocity[i] < vMin)

    if np.array(index_vMax).size > 0:
        velocity[i][index_vMax] = vMax[index_vMax]
    if np.array(index_vMin).size > 0:
        velocity[i][index_vMin] = vMin[index_vMin]

    # Actualiza posición de la partícula
    position[i] = position[i] + velocity[i] 

    # Ajusta posición de la particula
    index_pMax = np.where(position[i] > u_bounds)
    index_pMin = np.where(position[i] < l_bounds)

    if np.array(index_pMax).size > 0:
        position[i][index_pMax] = u_bounds[index_pMax]
    if np.array(index_pMin).size > 0:
        position[i][index_pMin] = l_bounds[index_pMin]

    # Evaluamos la función
    y = f(position[i])
    with gbest_val.get_lock():
        if y < gbest_val.value:
            with gbest_pos.get_lock(): 
                gbest_pos[:] = np.copy(position[i]) 
                pbest_pos[i] = np.copy(position[i])
                gbest_val.value = y
        if y < pbest_val[i]:
            pbest_pos[i] = np.copy(position[i])

def PSO(f_cost,pop_size,max_iters,lb,ub,α,β,w,w_max,w_min):
    # Tamaño de la población
    n = pop_size
    maxiter = max_iters
    # Número de variables
    n_var = len(lb)

    # Cognitive scaling parameter
    α = α
    # Social scaling parameter
    β = β

    # velocity inertia
    w = Value(ctypes.c_double,w)
    # minimum value for the velocity inertia
    w_min = w_min
    # maximum value for the velocity inertia
    w_max = w_max

    # Usamos Latin Hypercube Sampling para muestrear puntos en el espacio de búsqueda
    engine = LatinHypercube(d=n_var)
    sample = engine.random(n=n)

    # Definimos los límites superiores e inferiores para las variables de decisión
    l_bounds = np.array(lb)
    u_bounds = np.array(ub)

    # Creamos un arreglo compartido para el vector de limites superiores
    mp_l_bounds = Array(ctypes.c_double,l_bounds)
    # Creamos un nuevo arreglo de numpy usando el arreglo compartido
    np_l_bounds = np.frombuffer(mp_l_bounds.get_obj(), ctypes.c_double) 

    # Creamos un arreglo compartido para el vector de limites superiores
    mp_u_bounds = Array(ctypes.c_double,u_bounds)
    # Creamos un nuevo arreglo de numpy usando el arreglo compartido
    np_u_bounds = np.frombuffer(mp_u_bounds.get_obj(), ctypes.c_double) 

    # Velocidad máxima
    vMax = np.multiply(u_bounds-l_bounds,0.2)
    # Creamos un arreglo compartido para el vector de velocidad máxima
    mp_vMax = Array(ctypes.c_double,vMax) 
    # Creamos un nuevo arreglo de numpy usando el arreglo compartido
    np_vMax = np.frombuffer(mp_vMax.get_obj(), ctypes.c_double) 

    # Velocidad mínima
    vMin = -vMax
    # Creamos un arreglo compartido para el vector de velocidad máxima
    mp_vMin = Array(ctypes.c_double,vMin) 
    # Creamos un nuevo arreglo de numpy usando el arreglo compartido
    np_vMin = np.frombuffer(mp_vMin.get_obj(), ctypes.c_double) 


    # Escalamos los valores muestreados de LHS
    sample_scaled = scale(sample,l_bounds, u_bounds)

    # Creamos un arreglo compartido para el vector de velocidad
    mp_vel = Array(ctypes.c_double,n*n_var)
    # Creamos un nuevo arreglo de numpy usando el arreglo compartido
    vel = np.frombuffer(mp_vel.get_obj(), ctypes.c_double)
    # Convertimos a un arreglo 2-dimensional
    vel_resh = vel.reshape((n,n_var))

    # Creamos un arreglo compartido para el vector de posición
    mp_pos = Array(ctypes.c_double,n*n_var)
    # Creamos un nuevo arreglo de numpy usando el arreglo compartido
    pos = np.frombuffer(mp_pos.get_obj(), ctypes.c_double)
    # Convertimos a un arreglo 2-dimensional
    pos_resh = pos.reshape((n,n_var))
    # Inicializamos el vector de posición con el vector muestreado por LHS
    for i,v in enumerate(sample_scaled):
        pos_resh[i] = v

    # Mejor valor global (compartido) de la función objetivo
    gbest_val = Value(ctypes.c_double,math.inf)
    # Mejor vector de posición global (compartido)
    gbest_pos = Array(ctypes.c_double, sample_scaled[0])

    # Mejor valor para cada partícula
    pbest_val_arg = Array(ctypes.c_double, [math.inf]*n )

    # Mejor vector de posición individual para cada partícula
    pbest_pos_mp = Array(ctypes.c_double,n*n_var)
    # Creamos un nuevo arreglo de numpy usando el arreglo compartido
    pbest_pos = np.frombuffer(pbest_pos_mp.get_obj(), ctypes.c_double)
    # Convertimos a un arreglo 2-dimensional
    pbest_pos_arg = pbest_pos.reshape((n,n_var))
    # Inicializamos el vector de posición con el vector muestreado por LHS
    for i,v in enumerate(sample_scaled):
        pbest_pos_arg[i] = v

    p = Pool(processes = int(cpu_count()),initializer=init_pso,
            initargs=(gbest_val,gbest_pos,pos_resh, vel_resh, 
                      pbest_val_arg, pbest_pos_arg, f_cost,α,β,w,
                      np_vMax,np_vMin,np_u_bounds,np_l_bounds,))

    fitness_values = []
    for k in range(maxiter):
        p.map(evalua_f_pso, range(n))
        print("It {}  gbest_val {}".format(k, gbest_val.value))

        # Actualizamos w
        w.value = w_max - k * ((w_max-w_min)/maxiter)

        fitness_values.append(gbest_val.value)

    return fitness_values, gbest_pos[:], gbest_val.value