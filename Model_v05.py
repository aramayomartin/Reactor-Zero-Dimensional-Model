
#||||||||||||||||||||||||||||||||||||||||||||||||||
#----------------- MODELO DE CAREM ----------------
#||||||||||||||||||||||||||||||||||||||||||||||||||

#-------------------- Librerias ------------------
import numpy as np
import math
import matplotlib.pyplot as plt
from iapws import IAPWS97 as ip

#--------------------------------------------------
#---------------------- DATOS ---------------------
#--------------------------------------------------

#-------------------- Primario --------------------

##Presión nominal: 12.25 MPa
p_prim_nom = 12.25e6 #Pa

lsat_prim_nom = ip(P=p_prim_nom/1e6,x=0)
vsat_prim_nom = ip(P=p_prim_nom/1e6,x=1)

##Temperatura nominal de saturación (K)
T_sat_nom_prim = lsat_prim_nom.T
##Potencia generada: 100 MWth
Q_gen_nom= 1e8 #W
##Volumen neto total (líquido+vapor): 55.9 m3
Vol_tot_prim = 55.9 #m3
##Volumen de agua hasta tope de núcleo: 13.31 m3
Vol_l_tope_Nuc = 13.31 #m3
##Volumen de vapor: 10.4 m3
Vol_g_prim = 10.4 #m3

vf_nom = lsat_prim_nom.v
vg_nom = vsat_prim_nom.v

mg_nom = Vol_g_prim/vg_nom
mf_nom = (Vol_tot_prim - Vol_g_prim)/vf_nom

m_tot_prim_nom = mf_nom + mg_nom

nivel_tope_nuc = Vol_l_tope_Nuc/Vol_tot_prim

#---------------------- SIE ---------------------

##Volumen agua SIE: 50 m3
Vol_l_SIE_nom = 50 #m3
##Volumen N2 SIE: 19 m3
Vol_N2_nom = 19 #m3
##Temperatura inicial del SIE: 45 ºC
T_SIE_nom = 273.15+45 #K
##Factor de fricción descarga SIE: 2.47*104 (en un área de 4.53*10-3m2)
##K_SIE = 2.47e4 (original de la guía)
K_SIE = 2.47e4 #ajustar para que el caudal inyectado matenga el nucleo cubierto
##Area de tubería: 4.53e-3 m2
A_SIE = 4.53e-3 #m2
##Presión de SIE: 2 MPa
p_N2_nom = 2e6 #Pa
##Diferencia de presión para ruptura disco de ruptura: 0.5 MPa
dP_disco_SIE = 0.5e6 #Pa

#-----------------  Combustibles ------------------

##Temperatura promedio de combustible: 600 ºC
T_comb_nom = 273.15+600 #K
##Radio pastilla (r1): 3.8 mm
r1 = 3.8e-3 #m
##Radio interno vaina (r2): 3.875 mm
r2 = 3.875e-3 #m
##Radio externo vaina (r3): 4.5 mm
r3 = 4.5e-3 #m
##Número total de barras combustibles: 6583
N_b = 6583
##Longitud zona activa: 1.4 m
L_act = 1.4 #m
##Capacidad calorífica UO2: 3*106 J/(m3*K)
rhoCp_UO2 = 3e6 #J/(m3*K)
##Capacidad calorífica vaina: 2*106 J/(m3*K)
rhoCp_vaina = 2e6 #J/(m3*K)

mCp_comb = (rhoCp_UO2*r1**2 + rhoCp_vaina*(r3**2-r2**2))*np.pi*L_act*N_b #J/K
hA_comb = Q_gen_nom/(T_comb_nom-T_sat_nom_prim) #W/K


#------------------------ GV -------------------------
##Potencia nominal extraida: 100MW
Q_GV_nom = 1e8 #W
##Presión de secundario : 47 bar
p_sec_nom = 4.7e6 #Pa
##Masa de líquido: 600 kg
m_tot_GV = 600 #kg

T_sat_nom_sec = ip(P=p_sec_nom/1e6,x=0).T

hf_GV = ip(P=p_sec_nom/1e6,x=0).h*1e3 #J/kg
hg_GV = ip(P=p_sec_nom/1e6,x=1).h*1e3 #J/kg
hfg_GV = hg_GV - hf_GV

hA_GV = Q_GV_nom/(T_sat_nom_prim - T_sat_nom_sec) #W/K

#--------------------- Estructuras -------------------
##Masa RPV: 143 ton
m_RPV = 143e3 #kg
##Masa internos RPV: 53 ton
m_internos_RPV = 53e3 #kg

m_estr = m_RPV + m_internos_RPV #kg
Cp_estr = 480 #J/(kgK) Acero Inoxidable AISI302 (Incropera 4ta ed)
mCp_estr = m_estr*Cp_estr
h_conv_nat_gas_ref = 5 #W/(m2K) (valor de referencia Kern)
h_cnv_nat_liq_ref = 500 #W/(m2K) (valor de referencia Kern)
frac_vap_nom = Vol_g_prim/Vol_tot_prim
frac_liq_nom = 1 - frac_vap_nom
h_estr = frac_vap_nom*h_conv_nat_gas_ref + frac_liq_nom*h_cnv_nat_liq_ref 
#(Calculada como la sup de un cilindro de 3x11 + 20% de internos)
A_estr = (np.pi*(3)*11+2*np.pi*(1.5)**2)*1.2 #m2 
hA_estr = h_estr*A_estr

#--------------------- SECR -----------------------------------

#Potencia nominal extraida por el SECR: 2 MW
Q_SECR_nom = 2e6 #W 
#Presión del SECR es la P de la contensión: 1 atm
p_cont = 101325 #Pa
p_SECR = p_cont

T_sat_SECR = ip(P=p_SECR/1e6, x=0).T

hA_SECR = Q_SECR_nom/(T_sat_nom_prim - T_sat_SECR)


#------------------------ Funciones ---------------------------

def SCRAM(p_prim):
    global trip_SCRAM
    global t_SCRAM

    p_inf = p_umbral_inf_SCRAM
    p_sup = p_umbral_sup_SCRAM
    
    if not trip_SCRAM:
        if p_prim >= p_sup or  p_prim <= p_inf or t >= t_SCRAM_manual:
            trip_SCRAM = True
            t_SCRAM = t
            print("se demanda SCRAM a t = {:0.2f} s".format(t_SCRAM))
    return None

def LOCA(t_LOCA):
    global t
    global trip_LOCA
    if t_LOCA == None:
        return None
    if not trip_LOCA and t >= t_LOCA:
        trip_LOCA = True
        print("se inicia el evento LOCA a t = {:0.2f} s".format(t_LOCA))
    return None

def LOHS(t_LOHS):
    global t
    global trip_LOHS
    if t_LOHS == None:
        return None
    if not trip_LOHS and t >= t_LOHS:
        trip_LOHS = True
        print("se inicia el evento LOHS a t = {:0.2f} s".format(t_LOHS))
    return None

def SIE(p_prim):
    global trip_SIE
    dp = p_N2_nom - p_prim 
    if not trip_SIE and dp>dP_disco_SIE:
        trip_SIE = True
        print("se demanda el SIE a t = {:0.2f} s".format(t))
    return None

def SECR(p_prim, T_prim):
    global trip_SECR

    p_inf = p_umbral_inf_SECR
    p_sup = p_umbral_sup_SECR
    T_sup = T_umbral_sup_SECR
    
    if not trip_SECR:
        if p_prim >= p_sup or p_prim <= p_inf or T_prim >= T_sup:
            trip_SECR = True
            print("se demanda el SECR a t = {:0.2f} s".format(t))
    return None

def pot_nuc(t, p_prim):

    if not trip_SCRAM:
        SCRAM(p_prim)

    if trip_SCRAM:
        t = t - t_SCRAM
        if t < 9:
            Q_gen = 0.15*Q_gen_nom*math.e**(-0.1*t)
        elif 9 <= t <= 180:
            Q_gen = 0.09192*Q_gen_nom*t**(-0.181)
        elif 180 < t < 4e6:
            Q_gen = 0.156*Q_gen_nom*t**(-0.283)
    else:
        Q_gen = Q_gen_nom

    return Q_gen


def calc_A_B(p_prim):

    liq_sat = ip(P=p_prim/1e6,x=0)
    vap_sat = ip(P=p_prim/1e6,x=1)
    
    u_f = liq_sat.u*1000 # J/Kg
    u_g = vap_sat.u*1000 # J/Kg
    v_f = liq_sat.v      # m3/Kg
    v_g = vap_sat.v      # m3/Kg

    A = (u_f - u_g)/(v_f - v_g) # J/m3
    B =  (u_g*v_f - u_f*v_g)/(v_f - v_g) # J/Kg
    
    return (A,B)

def calc_dXdp(p, dp):
    X_0 = calc_A_B(p)
    X_1 = calc_A_B(p+dp)
    
    dAdp = (X_1[0] - X_0[0])/dp

    dBdp = (X_1[1] - X_0[1])/dp
    
    return (dAdp, dBdp)

def calc_mflow_LOCA(p_prim):
    d_LOCA = 1.5*2.54/100 # m (1.5 inch)
    A_LOCA = np.pi*(d_LOCA/2)**2 # m2 
    rhog = ip(P=p_prim/1e6,x=1).rho
    K_LOCA = 1.46 # Contracción y expansion abrupta (verificar!)
    K_bloq = 1.3
    exp_bloq = (K_bloq+1)/(K_bloq-1)

    mflow_Bloq   = A_LOCA*(K_bloq*p_prim*rhog*(2/(K_bloq+1))**exp_bloq)**0.5
    mflow_noBloq = A_LOCA*(2*(p_prim - p_cont)*rhog/K_LOCA)**0.5

    mflow_LOCA = min(mflow_noBloq, mflow_Bloq)

##    if p_prim >= 2*p_cont:
##        mflow_LOCA = mflow_Bloq
##    elif p_prim < 2*p_cont:
##        mflow_LOCA = mflow_noBloq
        
    return mflow_LOCA

def calc_nivel_prim(p_prim, m_tot_prim_):
    
    v_f = ip(P=p_prim/1e6,x=0).v   # m3/Kg
    v_g = ip(P=p_prim/1e6,x=1).v   # m3/Kg

    Volf_prim = v_f*(Vol_tot_prim-m_tot_prim_*v_g)/(v_f - v_g)

    nivel_prim = Volf_prim/Vol_tot_prim

    return nivel_prim

#--------------------------------------------------------------
#------------------------ Esquema Numerico --------------------
#--------------------------------------------------------------

# Valores de paso de tiempo
dt = 1 #s
t_sim = 36*3600+100 #s
N = int(t_sim/dt)

print("tiempo de simulación t = {} s".format(t_sim))
print("paso de tiempo dt = {} s".format(dt))

# Parámetros de demanda del SECR
p_umbral_inf_SECR = 0.5*p_prim_nom
p_umbral_sup_SECR = 13.7e6 #13.7 MPa
T_umbral_sup_SECR = 1.1*T_sat_nom_prim

p_inf_SECR_v = np.array([p_umbral_inf_SECR for n in range(N)])
p_sup_SECR_v = np.array([p_umbral_sup_SECR for n in range(N)])
T_sup_SECR_v = np.array([T_umbral_sup_SECR for n in range(N)])


# Parámetros de demanda de SCRAM
p_umbral_inf_SCRAM = 0.98*p_prim_nom 
p_umbral_sup_SCRAM = 13e6 # 13 MPa 

p_inf_SCRAM_v = np.array([p_umbral_inf_SCRAM for n in range(N)])
p_sup_SCRAM_v = np.array([p_umbral_sup_SCRAM for n in range(N)])

# Parámetros de demanda del SIE
p_umbral_SIE = p_N2_nom - dP_disco_SIE
p_SIE_v = np.array([p_umbral_SIE for n in range(N)])

# Nivel limite de descubrimiento de Núcleo
nivel_tope_nuc_v = np.array([nivel_tope_nuc for n in range(N)])

# Valores iniciales para la corrida
trip_SECR  = False
trip_SIE   = False
trip_LOHS  = False
t_LOHS     = 100 # Momento en el que se produce el LOHS (colocar "None" en caso de No LOHS)
trip_LOCA  = False
t_LOCA     = None # Momento al que se produce el LOCA (colocar "None" en caso de No LOCA)
trip_SCRAM = False
t_SCRAM_manual = 99999e99 # Demanda manual de SCRAM por tiempo


t = 0
Err = 1
m_tot_prim_0 = m_tot_prim_nom
T_comb_0     = T_comb_nom
p_prim_0     = p_prim_nom
T_prim       = lsat_prim_nom.T
T_estr_0     = T_prim
m_GV_0       = 0
Vol_N2_0     = Vol_N2_nom
p_N2         = p_N2_nom
rhof_SIE     = ip(T=T_SIE_nom,P=p_N2_nom/1e6).rho
nivel_prim   = (Vol_tot_prim - Vol_g_prim)/Vol_tot_prim


t_v            = []
p_prim_v       = []

m_GV_v         = []
m_tot_prim_v   = []

T_prim_v       = []
T_comb_v       = []
T_estr_v       = []
T_GV_v         = []

Q_gen_v        = []
Q_neto_v       = []
Q_transf_v     = []
Q_GV_v         = []
Q_comb_v       = []
Q_estr_v       = []
Q_SECR_v       = []

mflow_h_LOCA_v = []
mflow_h_SIE_v  = []

mflow_LOCA_v   = []
mflow_SIE_v    = []

nivel_prim_v   = []
nivel_GV_v     = []



for i in range(N):
    SECR(p_prim_0, T_prim)
    SIE(p_prim_0)
    LOCA(t_LOCA)
    LOHS(t_LOHS)
    Q_gen = pot_nuc(t, p_prim_0)
    
    while Err > 0.01:

        liq_sat_prim = ip(P=p_prim_0/1e6,x=0)
        vap_sat_prim = ip(P=p_prim_0/1e6,x=1)

        T_prim = liq_sat_prim.T #Carga la forma funcional de la temperatura en el primario con la biblioteca
        
        Q_transf = hA_comb*(T_comb_0 - T_prim)
        
        if trip_SCRAM or trip_LOHS:
            r = (m_tot_GV - m_GV_0)/m_tot_GV
            r = max(0, r)
            r = min(1, r)
        else:
            r = 1

        dT_prim_GV = T_sat_nom_sec - T_prim
        if dT_prim_GV < 0:
            Q_GV = hA_GV*r*(T_sat_nom_sec - T_prim)
        else:
            Q_GV = 0
            

        Q_estr = hA_estr*(T_estr_0 - T_prim)

        Q_comb = Q_gen - Q_transf

        if trip_SECR:
            Q_SECR = hA_SECR*(T_sat_SECR - T_prim)
        else:
            Q_SECR = 0

        if trip_SIE:
            gamma     = 1.3 # expansión adiabatica
##            gamma     = 1   # expansión isotérmica
            p_N2      = p_N2_nom*(Vol_N2_nom/Vol_N2_0)**gamma
            rhof_SIE  = ip(T=T_SIE_nom,P=p_N2/1e6).rho
            mflow_SIE = A_SIE*(2*(p_N2 - p_prim_0)*rhof_SIE/K_SIE)**0.5
            Vol_N2_1  = Vol_N2_0 + mflow_SIE*dt/rhof_SIE
            
        else:
            mflow_SIE = 0

        if trip_LOCA:
            mflow_LOCA = calc_mflow_LOCA(p_prim_0)
        else:
            mflow_LOCA = 0

        mflow_h_LOCA = -mflow_LOCA*vap_sat_prim.h*1e3             #Flujo entalpico de LOCA
        mflow_h_SIE  = mflow_SIE*ip(T=T_SIE_nom,P=p_N2/1e6).h*1e3 #Flujo entalpico de SIE
        
        Q_neto = Q_transf + Q_GV + Q_estr + Q_SECR + mflow_h_LOCA + mflow_h_SIE

        dAdp,dBdp = calc_dXdp(p_prim_0,p_prim_0/1000)

        mflow_tot_prim = mflow_SIE - mflow_LOCA

        B = calc_A_B(p_prim_0)[1]
       
        dpdt = (Q_neto - mflow_tot_prim*B)/(Vol_tot_prim*dAdp + m_tot_prim_0*dBdp)
        
        p_prim_1 = p_prim_0 + dpdt*dt

        T_estr_1 = T_estr_0 - Q_estr*dt/mCp_estr

        T_comb_1 = T_comb_0 + Q_comb*dt/mCp_comb

        m_tot_prim_1 = m_tot_prim_0 + mflow_tot_prim*dt

        m_tot_prim_1 = max(0, m_tot_prim_1)
        m_tot_prim_1 = min(m_tot_prim_nom, m_tot_prim_1)

        if trip_SCRAM or trip_LOHS:
            m_GV_1 = m_GV_0 - Q_GV*dt/hfg_GV
        else:
            m_GV_1 = 0

        Vol_N2_1  = Vol_N2_0 + mflow_SIE*dt/rhof_SIE

        Err = abs(p_prim_1 - p_prim_0)/p_prim_0

        p_prim_0     = p_prim_1
        T_estr_0     = T_estr_1
        T_comb_0     = T_comb_1
        m_tot_prim_0 = m_tot_prim_1
        m_GV_0       = m_GV_1
        Vol_N2_0     = Vol_N2_1

    nivel_prim = calc_nivel_prim(p_prim_1, m_tot_prim_1)
    
    t_v.append(t)
    p_prim_v.append(p_prim_1)

    m_GV_v.append(m_GV_1)
    m_tot_prim_v.append(m_tot_prim_1)
    
    T_prim_v.append(T_prim)
    T_comb_v.append(T_comb_1) 
    T_estr_v.append(T_estr_1)
    T_GV_v.append(T_sat_nom_sec)
    
    Q_gen_v.append(Q_gen)
    Q_neto_v.append(Q_neto)
    Q_transf_v.append(Q_transf)
    Q_GV_v.append(-Q_GV)
    Q_comb_v.append(-Q_comb)
    Q_estr_v.append(Q_estr)
    Q_SECR_v.append(-Q_SECR)

    mflow_h_LOCA_v.append(-mflow_h_LOCA)
    mflow_h_SIE_v.append(mflow_h_SIE)
    
    mflow_LOCA_v.append(mflow_LOCA)
    mflow_SIE_v.append(mflow_SIE)
    
    nivel_GV_v.append(r)
    nivel_prim_v.append(nivel_prim)
    
    t += dt
    Err = 1
    
t_v            = np.array(t_v)	#Aqui cambié una division por 3600
p_prim_v       = np.array(p_prim_v)

m_GV_v         = np.array(m_GV_v)
m_tot_prim_v   = np.array(m_tot_prim_v)

T_prim_v       = np.array(T_prim_v)
T_comb_v       = np.array(T_comb_v)
T_estr_v       = np.array(T_estr_v)
T_GV_v         = np.array(T_GV_v)

Q_gen_v        = np.array(Q_gen_v)
Q_neto_v       = np.array(Q_neto_v)
Q_transf_v     = np.array(Q_transf_v)
Q_GV_v         = np.array(Q_GV_v)
Q_comb_v       = np.array(Q_comb_v)
Q_estr_v       = np.array(Q_estr_v)
Q_SECR_v       = np.array(Q_SECR_v)

mflow_h_LOCA_v = np.array(mflow_h_LOCA_v)
mflow_h_SIE_v  = np.array(mflow_h_SIE_v)

#------------------------------ Gráficos ------------------------------------
plt.figure(0)
plt.grid()
plt.plot(t_v, p_prim_v/1e6,            label = 'Presión del primario')
plt.plot(t_v, p_sup_SECR_v/1e6,  '--', label = 'demanda de SECR por alta p')
plt.plot(t_v, p_sup_SCRAM_v/1e6, '--', label = 'demanda de SCRAM por alta p')
#plt.plot(t_v, p_inf_SCRAM_v/1e6, '--', label = 'demanda de SCRAM por baja p')
#plt.plot(t_v, p_inf_SECR_v/1e6,  '--', label = 'demanda de SECR por baja p')
#plt.plot(t_v, p_SIE_v/1e6,       '--', label = 'demanda de SIE por baja p')
plt.legend(loc="best", shadow=True, fontsize = 'large')
##plt.xlabel('Tiempo (s)', fontsize = 'x-large')
plt.xlabel('Tiempo (s)', fontsize = 'x-large')
plt.ylabel('Presión (MPa)', fontsize = 'x-large')

plt.figure(1)
plt.grid() 
plt.plot(t_v, Q_neto_v/1e6,           label = 'Q Neto')
plt.plot(t_v, Q_transf_v/1e6,         label = 'Q Nucleo transferido')
plt.plot(t_v, Q_GV_v/1e6, '--',       label = 'Q GV')
plt.plot(t_v, Q_estr_v/1e6, '--',     label = 'Q Estructuras')
plt.plot(t_v, Q_SECR_v/1e6,           label = 'Q SECR')
if trip_LOCA:
    plt.plot(t_v, mflow_h_LOCA_v/1e6,'--',label = 'Flujo entálpico LOCA')
plt.plot(t_v, mflow_h_SIE_v/1e6,'--', label = 'Flujo entálpico SIE')
plt.legend(loc="best", shadow=True, fontsize = 'large')
##plt.xlabel('Tiempo (s)', fontsize = 'x-large')
plt.xlabel('Tiempo (s)', fontsize = 'x-large')
plt.ylabel('Q (MW)', fontsize = 'x-large')

plt.figure(2)
plt.grid() 
plt.plot(t_v, Q_gen_v/1e6,    label = 'Q_gen')
plt.plot(t_v, Q_transf_v/1e6,'--', label = 'Q_transf')
plt.plot(t_v, Q_comb_v/1e6,   label = 'Q_comb')
plt.legend(loc="best", shadow=True, fontsize = 'large')
##plt.xlabel('Tiempo (s)', fontsize = 'x-large')
plt.xlabel('Tiempo (s)', fontsize = 'x-large')
plt.ylabel('Q (MW)', fontsize = 'x-large')

plt.figure(3)
plt.grid()
plt.plot(t_v, T_prim_v-273.15, label = 'T_prim')
plt.plot(t_v, T_comb_v-273.15, label = 'T_comb')
plt.plot(t_v, T_estr_v-273.15,'--', label = 'T_estr')
plt.plot(t_v, T_GV_v-273.15,   label = 'T_GV')
plt.legend(loc="best", shadow=True, fontsize = 'large')
##plt.xlabel('Tiempo (s)', fontsize = 'x-large')
plt.xlabel('Tiempo (s)', fontsize = 'x-large')
plt.ylabel('Temp (ºC)', fontsize = 'x-large')

plt.figure(4)
plt.grid()
plt.plot(t_v, mflow_LOCA_v, label = 'mflow_LOCA')
plt.plot(t_v, mflow_SIE_v, label = 'mflow_SIE')
plt.legend(loc="best", shadow=True, fontsize = 'large')
##plt.xlabel('Tiempo (s)', fontsize = 'x-large')
plt.xlabel('Tiempo (s)', fontsize = 'x-large')
plt.ylabel('mflow (Kg/s)', fontsize = 'x-large')

plt.figure(5)
plt.grid()
plt.plot(t_v, nivel_prim_v,          label = 'nivel del primario colapsado')
plt.plot(t_v, nivel_tope_nuc_v,'--', label = 'nivel límite de descubrimiento de Núcleo')
plt.legend(loc="best", shadow=True, fontsize = 'large')
##plt.xlabel('Tiempo (s)', fontsize = 'x-large')
plt.xlabel('Tiempo (s)', fontsize = 'x-large')
plt.ylabel('fracción de nivel', fontsize = 'x-large')

plt.figure(6)
plt.grid()
plt.plot(t_v, nivel_GV_v, label = 'nivel_GV')
plt.legend(loc="best", shadow=True, fontsize = 'large')
##plt.xlabel('Tiempo (s)', fontsize = 'x-large')
plt.xlabel('Tiempo (s)', fontsize = 'x-large')
plt.ylabel('fracción de nivel', fontsize = 'x-large')

plt.figure(7)
plt.grid()
plt.plot(t_v, m_GV_v, label = 'm_GV')
plt.legend(loc="best", shadow=True, fontsize = 'large')
##plt.xlabel('Tiempo (s)', fontsize = 'x-large')
plt.xlabel('Tiempo (s)', fontsize = 'x-large')
plt.ylabel('masa (kg)', fontsize = 'x-large')

plt.show()


    

