import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import odeint
import pandas as pd
import os
import re
import matplotlib.colors as mcolors

Vmax1_B = 8.46 *10**(-6)
KM1_B = 1.13
kp1_B = 0.0005
kp2_B = 0.007
kp3_B = 0.007
ke1_B = 0.0104
ke2_B = 0.0104
ke3_B = 0.0104
ke4_B = 0.518
k12_B = 0.140
k23_B = 0.191
k34_B = 0.355
k5_B = 0
K0_B = 3.85*10**6
KM5_B = 204*10**3
Vmax51_B = 2570
Vmax52_B = 4040
Vmax53_B = 3780
Vmax54_B = 4240
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
K0_BM = 1.15*10**6
KM5_BM = 72*10**3
Vmax51_BM = Vmax51_B
Vmax52_BM = Vmax52_B
Vmax53_BM = Vmax53_B
Vmax54_BM = Vmax54_B
k1_BBM = 0.11
k2_BBM = k1_BBM
k3_BBM = k1_BBM
k4_BBM = k1_BBM
kCD19_BBM = 0
k1_BMB = 0.176
k2_BMB = k1_BMB
k3_BMB = k1_BMB
k4_BMB = k1_BMB
kCD19_BMB = 0.001
          
#set tspan
tStart=0
tEnd=200

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
        
    dTCM_B = (Vmax1_B*CD19_B*TCM_B)/(KM1_B+TCM_B) + kp2_B*TCM_B + k12_B*TN_B - k23_B*TCM_B - ke2_B*TCM_B + k2_BMB*TCM_BM - k2_BBM*TCM_B
    
    dTEM_B = (Vmax1_B*CD19_B*TEM_B)/(KM1_B+TEM_B) + kp3_B*TEM_B + k23_B*TCM_B - k34_B*TEM_B - ke3_B*TEM_B + k3_BMB*TEM_BM - k3_BBM*TEM_B
    
    dTE_B =  k34_B*TEM_B - ke4_B*TE_B + k4_BMB*TE_BM - k4_BBM*TE_B
    
    dCD19_B = k5_B*(1-CD19_B/K0_B)*CD19_B - (Vmax51_B*TN_B*CD19_B)/(KM5_B+CD19_B) - (Vmax52_B*TCM_B*CD19_B)/(KM5_B+CD19_B) - (Vmax53_B*TEM_B*CD19_B)/(KM5_B+CD19_B) - (Vmax54_B*TE_B*CD19_B)/(KM5_B+CD19_B) + kCD19_BMB*CD19_BM - kCD19_BBM*CD19_B
    
    dTN_BM = (Vmax1_BM*CD19_BM*TN_BM)/(KM1_BM+TN_BM) + kp1_BM*TN_BM - k12_BM*TN_BM - ke1_BM*TN_BM - k1_BMB*TN_BM + k1_BBM*TN_B
        
    dTCM_BM = (Vmax1_BM*CD19_BM*TCM_BM)/(KM1_BM+TCM_BM) + kp2_BM*TCM_BM + k12_BM*TN_BM - k23_BM*TCM_BM - ke2_BM*TCM_BM - k2_BMB*TCM_BM + k2_BBM*TCM_B
    
    dTEM_BM = (Vmax1_BM*CD19_BM*TEM_BM)/(KM1_BM+TEM_BM) + kp3_BM*TEM_BM + k23_BM*TCM_BM - k34_BM*TEM_BM - ke3_BM*TEM_BM - k3_BMB*TEM_BM + k3_BBM*TEM_B
    
    dTE_BM =  k34_BM*TEM_BM - ke4_BM*TE_BM - k4_BMB*TE_BM + k4_BBM*TE_B
    
    dCD19_BM = k5_BM*(1-CD19_BM/K0_BM)*CD19_BM - (Vmax51_BM*TN_BM*CD19_BM)/(KM5_BM+CD19_BM) - (Vmax52_BM*TCM_BM*CD19_BM)/(KM5_BM+CD19_BM) - (Vmax53_BM*TEM_BM*CD19_BM)/(KM5_BM+CD19_BM) - (Vmax54_BM*TE_BM*CD19_BM)/(KM5_BM+CD19_BM) - kCD19_BMB*CD19_BM + kCD19_BBM*CD19_B
    
    dY = np.array([dTN_B, dTCM_B, dTEM_B, dTE_B, dCD19_B, dTN_BM, dTCM_BM, dTEM_BM, dTE_BM, dCD19_BM])
    return dY

def separate_dose(t_start, t_mid, t_end):
    TN_B_0 = 14
    TCM_B_0 = 0
    TEM_B_0 = 0
    TE_B_0 = 0
    CD19_B_0 = 0
    TN_BM_0 = 0
    TCM_BM_0 = 0
    TEM_BM_0 = 0
    TE_BM_0 = 0
    CD19_BM_0 = 100*10**3
    h = 0.1
    t = np.arange(t_start, t_mid, h)
    y0_dose1 = np.array([TN_B_0, TCM_B_0, TEM_B_0, TE_B_0, CD19_B_0, TN_BM_0, TCM_BM_0, TEM_BM_0, TE_BM_0, CD19_BM_0])
    sol_dose1 = odeint(CART_PBDK, y0_dose1, t, args=params)
    steps_num = int(t_mid//h)
    TN_B_2 = 14
    TN_B_1 = sol_dose1[steps_num, 0] + TN_B_2
    TCM_B_1 = sol_dose1[steps_num, 1]
    TEM_B_1 = sol_dose1[steps_num, 2]
    TE_B_1 = sol_dose1[steps_num, 3]
    CD19_B_1 = sol_dose1[steps_num, 4]
    TN_BM_1 = sol_dose1[steps_num, 5]
    TCM_BM_1 = sol_dose1[steps_num, 6]
    TEM_BM_1 = sol_dose1[steps_num, 7]
    TE_BM_1 = sol_dose1[steps_num, 8]
    CD19_BM_1 = sol_dose1[steps_num, 9]
    t = np.arange(t_mid, t_end, h)
    y0_dose2 = np.array([TN_B_1, TCM_B_1, TEM_B_1, TE_B_1, CD19_B_1, TN_BM_1, TCM_BM_1, TEM_BM_1, TE_BM_1, CD19_BM_1])
    sol_dose2 = odeint(CART_PBDK, y0_dose2, t, args=params)
    t = np.arange(t_start, t_end, h)
    print(t)
    print(sol_dose1)
    print(sol_dose2)
    sol = np.concatenate((sol_dose1, sol_dose2), axis=0)
    df_sol = pd.DataFrame(sol, columns=['TN_B', 'TCM_B', 'TEM_B', 'TE_B', 'CD19_B', 'TN_BM', 'TCM_BM', 'TEM_BM', 'TE_BM', 'CD19_BM'])
    df_t = pd.DataFrame(t, columns=['t'])
    df = pd.concat([df_t, df_sol], axis=1)
    print(df)
    df.to_csv(f'doubledose_1times-1/tumour_{CD19_BM_0}-TN_{TN_B_0}_{TN_B_2}-TCM_{TCM_B_0}-TEM_{TEM_B_0}-TE_{TE_B_0}-t_{t_start}_{t_mid}_{t_end}', index=False)
    return t, sol


def plot_from_csv():
    colour = mcolors.XKCD_COLORS
    print(colour)
    # legend_colour_list = ['grey', 'purple', 'green', 'red', 'cyan', 'orange', 'pink', 'blue']
    # legend_colour_list = ['blue', 'pink', 'orange', 'green', 'grey', 'cyan', 'purple', 'red']
    # legend_colour_list = ['cyan', 'pink', 'red', 'blue', 'orange', 'green', 'purple', 'grey']
    
    # legend_colour_list = ['blue', 'pink', 'cyan', 'red' , 'orange', 'purple', 'grey']
    legend_colour_list = ['orange', 'purple', 'red' , 'cyan', 'blue', 'pink', 'grey']
    legend_colour_list = ['red', 'cyan', 'blue' , 'orange', 'purple', 'pink', 'grey']
    file_list = os.listdir('doubledose_1times-1/')    
    legend_list = []
    i=0
    for file_name in file_list:
        print(file_name)
        df = pd.read_csv(f'doubledose_1times-1/{file_name}')
            # Regex pattern to extract the values
        pattern = r"tumour_(\d+)-TN_(\d+)_(\d+)-TCM_(\d+)-TEM_(\d+)-TE_(\d+)-t_(\d+)_(\d+)_(\d+)"

        match = re.match(pattern, file_name)
        print(match)
        extracted_values = match.groups()
    
        # Original variable names based on your format
        variable_names = [
            'CD19_BM_0', 'TN_B_0', 'TN_B_2', 'TCM_B_0',
            'TEM_B_0', 'TE_B_0', 't_start', 't_mid', 't_end'
            ]
    
        # Create a dictionary mapping variable names to their extracted values
        variables_dict = {name: int(value) for name, value in zip(variable_names, extracted_values)}
        print(variables_dict)
        
        t = df['t']
        CD19_BM = df['CD19_BM']
        CD19_B = df['CD19_B']
        TE_B = df['TE_B']
        TE_BM = df['TE_BM']
        
        line_name = f'{variables_dict['TN_B_0']}_{variables_dict['TCM_B_0']}_{variables_dict['TEM_B_0']}_{variables_dict['TE_B_0']}'
        print(variables_dict['TN_B_0'])
        legend_list.append(line_name)
        
        # if variables_dict['TN_B_0'] == 28 or variables_dict['TCM_B_0'] == 28 or variables_dict['TEM_B_0'] == 28 or variables_dict['TE_B_0'] == 28:
        #     plt.plot(t[:300], CD19_BM[:300]/1000, colour['xkcd:'+legend_colour_list[i]])
        # else:
        #     plt.plot(t[:300], CD19_BM[:300]/1000, colour['xkcd:'+legend_colour_list[i]], linestyle='dashed')
        
        if variables_dict['TN_B_0'] == 28:
            line_name = f'single dose'
        else:
            line_name = f'second dose given on day {variables_dict['t_mid']}'
        plt.plot(t[:500], CD19_BM[:500]/1000, colour['xkcd:'+legend_colour_list[i]])
        legend_list.append(line_name)
        
        i+=1
        
    # plt.legend(loc='best',labels=legend_list)
    plt.xlabel('Time (days)')
    plt.ylabel('CD19+ Tumour Volume in Bone Marrow (mL)')
    # plt.ylabel('CD19+ Tumour Volume in Blood (mL)')
    plt.grid()
    # plt.show()
    plt.savefig(f'/Users/weiheli/mini-project/G6_mini-project/plots/results_graph_1times_doubledose_CD19BM_nolegend1-1')
    plt.close()

def plot_sol(subject_list, i, sol_1):
    plt.plot(t, sol_1[:, i], 'b')
    # plt.plot(t, sol_2[:, i], 'r')
    # plt.plot(t, sol_3[:, i], 'g')
    plt.legend(loc='best',labels=['sol1', 'sol2', 'sol3'])
    plt.xlabel('t')
    plt.grid()
    # plt.show()
    plt.savefig(f'/Users/weiheli/mini-project/G6_mini-project/plots/{subject_list[i]}')
    plt.close()

subject_list = ['TN_B', 'TCM_B', 'TEM_B', 'TE_B', 'CD19_B', 'TN_BM', 'TCM_BM', 'TEM_BM', 'TE_BM', 'CD19_BM']

# t, sol = separate_dose(0, 25, 200)
plot_from_csv()

# for i in range(10):
#     plot_sol(subject_list, i, sol)


TN_B_0 = 30
TCM_B_0 = 0
TEM_B_0 = 0
TE_B_0 = 0
CD19_B_0 = 0
TN_BM_0 = 0
TCM_BM_0 = 0
TEM_BM_0 = 0
TE_BM_0 = 0
CD19_BM_0 = 1*10**3

# y0_1 = np.array([TN_B_0, TCM_B_0, TEM_B_0, TE_B_0, CD19_B_0, TN_BM_0, TCM_BM_0, TEM_BM_0, TE_BM_0, CD19_BM_0])
# y0_2 = np.array([30, 20, 20, 30, CD19_B_0, TN_BM_0, TCM_BM_0, TEM_BM_0, TE_BM_0, CD19_BM_0])
# y0_3 = np.array([1, 30, 30, 30, CD19_B_0, TN_BM_0, TCM_BM_0, TEM_BM_0, TE_BM_0, CD19_BM_0])
# t = np.arange(tStart, tEnd, 1)
# sol_1 = odeint(CART_PBDK, y0_1, t, args=params)
# sol_2 = odeint(CART_PBDK, y0_2, t, args=params)
# sol_3 = odeint(CART_PBDK, y0_3, t, args=params)

# print(sum(sol_1[100:1000, 4]))
# print(sol_1[800])
# print(sol_1[600])
# print(sol_2[800])
# print(sol_2[600])
# print(sol_3[800])
# print(sol_3[600])


# subject_list = ['TN_B', 'TCM_B', 'TEM_B', 'TE_B', 'CD19_B', 'TN_BM', 'TCM_BM', 'TEM_BM', 'TE_BM', 'CD19_BM']

# for i in range(10):
#     plot_sol(subject_list,i, sol_1, sol_2, sol_3)


