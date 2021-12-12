import matplotlib.pyplot as plt
import numpy as np

t_v=np.loadtxt('tiempo.txt') #Cargo el vector de tiempo

N=len(t_v)
# Parámetros de demanda del SECR
p_umbral_inf_SECR = 10.7e6 # 10MPa
p_umbral_sup_SECR = 13.7e6 #13.7 MPa
T_umbral_sup_SECR = 380 + 273.15 # aprox 1.1*T_sat_nom_prim
p_inf_SECR_v = np.array([p_umbral_inf_SECR for n in range(N)])
p_sup_SECR_v = np.array([p_umbral_sup_SECR for n in range(N)])
T_sup_SECR_v = np.array([T_umbral_sup_SECR for n in range(N)])
	# Parámetros de demanda de SCRAM
p_umbral_inf_SCRAM = 11e6 # 11 MPa 
p_umbral_sup_SCRAM = 13e6 # 13 MPa 
p_inf_SCRAM_v = np.array([p_umbral_inf_SCRAM for n in range(N)])
p_sup_SCRAM_v = np.array([p_umbral_sup_SCRAM for n in range(N)])

	# Parámetros de demanda del SIE
p_N2_nom = 2e6 #Pa
dP_disco_SIE = 0.5e6 #Pa

p_umbral_SIE = p_N2_nom - dP_disco_SIE
p_SIE_v = np.array([p_umbral_SIE for n in range(N)])

t_v=np.loadtxt('tiempo.txt') #Cargo el vector de tiempo
p_prim_v_mt=np.loadtxt('masa_total.txt')
p_prim_v_m2=np.loadtxt('mitad_masa.txt')
p_prim_v_m3=np.loadtxt('masa40.txt')
p_prim_v_m4=np.loadtxt('masa60.txt')

plt.figure(0)
plt.grid()
plt.plot(t_v, p_prim_v_mt/1e6,            label = 'mt = 30 Ton')
plt.plot(t_v, p_prim_v_m2/1e6,  '--', label = 'mt = 20 Ton')
plt.plot(t_v, p_prim_v_m3/1e6,  '--', label = 'mt = 40 Ton')
plt.plot(t_v, p_prim_v_m4/1e6,  '--', label = 'mt = 60 Ton')

plt.plot(t_v, p_sup_SECR_v/1e6,  '--', label = 'demanda de SECR por alta p')
plt.plot(t_v, p_sup_SCRAM_v/1e6, '--', label = 'demanda de SCRAM por alta p')
#plt.plot(t_v, p_inf_SCRAM_v/1e6, '--', label = 'demanda de SCRAM por baja p')
#plt.plot(t_v, p_inf_SECR_v/1e6,  '--', label = 'demanda de SECR por baja p')
#plt.plot(t_v, p_SIE_v/1e6,       '--', label = 'demanda de SIE por baja p')

plt.legend(loc="best", shadow=True, fontsize = 'large')
##plt.xlabel('Tiempo (s)', fontsize = 'x-large')
plt.xlabel('Tiempo (h)', fontsize = 'x-large')
plt.ylabel('Presión (MPa)', fontsize = 'x-large')
plt.show()
