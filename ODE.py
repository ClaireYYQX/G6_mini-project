import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import odeint


Vmax1_B = 1.69 *10**(-12)
KM1_B = 1.13
kp1_B = 0.0005
kp2_B = 0.007
kp3_B = 0.007
ke1_B = 0
ke2_B = 0
ke3_B = 0
ke4_B = 0
k12_B = 0.0104
k23_B = 0.191
k34_B = 0.355
k5_B = 0
K0_B = 5*10**6
KM5_B = 276*10**3
Vmax51_B = 2.57
Vmax52_B = 4.04
Vmax53_B = 3.78
Vmax54_B = 4.24
Vmax1_BM = Vmax1_B
KM1_BM = 1.13
kp1_BM = kp1_B
kp2_BM = kp2_B
kp3_BM = kp3_B
ke1_BM = 0
ke2_BM = 0
ke3_BM = 0
ke4_BM = 0
k12_BM = k12_B
k23_BM = k23_B
k34_BM = k34_B
k5_BM = 0.0023
K0_BM = 5*10**6
KM5_BM = 276*10**3
Vmax51_BM = 2.57
Vmax52_BM = 4.04
Vmax53_BM = 3.78
Vmax54_BM = 4.24
k1_BBM = 0.4
k2_BBM = 0.4
k3_BBM = 0.4
k4_BBM = 0.4
kCD19_BBM = 0
k1_BMB = 0.15
k2_BMB = 0.15
k3_BMB = 0.15
k4_BMB = 0.15
kCD19_BMB = 0.001
          
#set tspan
tStart=0
tEnd=2000


params = (Vmax1_B,
          KM1_B,
          kp1_B,
          kp2_B,
          kp3_B,
          ke1_B,
          ke2_B,
          ke3_B,
          ke4_B,
          k12_B,
          k23_B,
          k34_B,
          k5_B,
          K0_B,
          KM5_B,
          Vmax51_B,
          Vmax52_B,
          Vmax53_B,
          Vmax54_B,
          Vmax1_BM,
          KM1_BM,
          kp1_BM,
          kp2_BM,
          kp3_BM,
          ke1_BM,
          ke2_BM,
          ke3_BM,
          ke4_BM,
          k12_BM,
          k23_BM,
          k34_BM,
          k5_BM,
          K0_BM,
          KM5_BM,
          Vmax51_BM,
          Vmax52_BM,
          Vmax53_BM,
          Vmax54_BM,
          k1_BBM,
          k2_BBM,
          k3_BBM,
          k4_BBM,
          kCD19_BBM,
          k1_BMB,
          k2_BMB,
          k3_BMB,
          k4_BMB,
          kCD19_BMB)



def CART_PBDK(Y, time, Vmax1_B, KM1_B, kp1_B, kp2_B, kp3_B, ke1_B, ke2_B, ke3_B, ke4_B, k12_B, k23_B, k34_B, k5_B, K0_B,
              KM5_B, Vmax51_B, Vmax52_B, Vmax53_B, Vmax54_B, Vmax1_BM, KM1_BM, kp1_BM, kp2_BM, kp3_BM, ke1_BM, ke2_BM,
              ke3_BM, ke4_BM, k12_BM, k23_BM, k34_BM, k5_BM, K0_BM, KM5_BM, Vmax51_BM, Vmax52_BM, Vmax53_BM, Vmax54_BM,
              k1_BBM, k2_BBM, k3_BBM, k4_BBM, kCD19_BBM, k1_BMB, k2_BMB, k3_BMB, k4_BMB, kCD19_BMB):
    
    TN_B, TCM_B, TEM_B, TE_B, CD19_B, TN_BM, TCM_BM, TEM_BM, TE_BM, CD19_BM = Y[0], Y[1], Y[2], Y[3], Y[4], Y[5], Y[6], Y[7], Y[8], Y[9]
    dTN_B = (Vmax1_B*CD19_B*TN_B)/(KM1_B+TN_B) + kp1_B*TN_B - k12_B*TN_B - ke1_B*TN_B + k1_BMB*TN_BM - k1_BBM*TN_B
        
    dTCM_B = (Vmax1_B*CD19_B*TCM_B)/(KM1_B+TCM_B) + kp2_B*TCM_B + k12_B*TN_B - k23_B*TCM_B - ke2_B*TCM_B 
    + k2_BMB*TCM_BM - k2_BBM*TCM_B
    
    dTEM_B = (Vmax1_B*CD19_B*TEM_B)/(KM1_B+TEM_B) + kp3_B*TEM_B + k23_B*TCM_B - k34_B*TEM_B - ke3_B*TEM_B 
    + k3_BMB*TEM_BM - k3_BBM*TEM_B
    
    dTE_B =  k34_B*TEM_B - ke4_B*TE_B + k4_BMB*TE_BM - k4_BBM*TE_B
    
    dCD19_B = k5_B*(1-CD19_B/K0_B)*CD19_B - (Vmax51_B*TN_B*CD19_B)/(KM5_B+CD19_B) 
    - (Vmax52_B*TCM_B*CD19_B)/(KM5_B+CD19_B) - (Vmax53_B*TEM_B*CD19_B)/(KM5_B+CD19_B) 
    - (Vmax54_B*TE_B*CD19_B)/(KM5_B+CD19_B) + kCD19_BMB*CD19_BM - kCD19_BBM*CD19_B
    
    dTN_BM = (Vmax1_BM*CD19_BM*TN_BM)/(KM1_BM+TN_BM) + kp1_BM*TN_BM - k12_BM*TN_BM - ke1_BM*TN_BM - k1_BMB*TN_BM + k1_BBM*TN_B
        
    dTCM_BM = (Vmax1_BM*CD19_BM*TCM_BM)/(KM1_BM+TCM_BM) + kp2_BM*TCM_BM + k12_BM*TN_BM - k23_BM*TCM_BM - ke2_BM*TCM_BM 
    - k2_BMB*TCM_BM + k2_BBM*TCM_B
    
    dTEM_BM = (Vmax1_BM*CD19_BM*TEM_BM)/(KM1_BM+TEM_BM) + kp3_BM*TEM_BM + k23_BM*TCM_BM - k34_BM*TEM_BM - ke3_BM*TEM_BM 
    - k3_BMB*TEM_BM + k3_BBM*TEM_B
    
    dTE_BM =  k34_BM*TEM_BM - ke4_BM*TE_BM - k4_BMB*TE_BM + k4_BBM*TE_B
    
    dCD19_BM = k5_BM*(1-CD19_BM/K0_BM)*CD19_BM - (Vmax51_BM*TN_BM*CD19_BM)/(KM5_BM+CD19_BM) 
    - (Vmax52_BM*TCM_BM*CD19_BM)/(KM5_BM+CD19_BM) - (Vmax53_BM*TEM_BM*CD19_BM)/(KM5_BM+CD19_BM) 
    - (Vmax54_BM*TE_BM*CD19_BM)/(KM5_BM+CD19_BM) - kCD19_BMB*CD19_BM + kCD19_BBM*CD19_B
    
    dY = np.array([dTN_B, dTCM_B, dTEM_B, dTE_B, dCD19_B, dTN_BM, dTCM_BM, dTEM_BM, dTE_BM, dCD19_BM])
    
    return dY

def plot_sol(subject_list, i):
    plt.plot(t, sol[:, i], 'b')
    plt.legend(loc='best')
    plt.xlabel('t')
    plt.ylabel(label=subject_list[i])
    plt.grid()
    plt.show()
    plt.savefig(f'/Users/weiheli/mini-project/G6_mini-project/plots/{subject_list[i]}')

TN_B_0 = 0
TCM_B_0 = 0
TEM_B_0 = 0
TE_B_0 = 0
CD19_B_0 = 0
TN_BM_0 = 0
TCM_BM_0 = 0
TEM_BM_0 = 0
TE_BM_0 = 0
CD19_BM_0 = 0

y0 = np.array([TN_B_0, TCM_B_0, TEM_B_0, TE_B_0, CD19_B_0, TN_BM_0, TCM_BM_0, TEM_BM_0, TE_BM_0, CD19_BM_0])
t = np.arange(tStart, tEnd, 0.1)
sol = odeint(CART_PBDK, y0, t, args=params)

subject_list = ['TN_B', 'TCM_B', 'TEM_B', 'TE_B', 'CD19_B', 'TN_BM', 'TCM_BM', 'TEM_BM', 'TE_BM', 'CD19_BM']

for i in range(10):
    plot_sol(subject_list,i)


