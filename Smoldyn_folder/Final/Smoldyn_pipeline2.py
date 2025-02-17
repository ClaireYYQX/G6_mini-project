#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar  7 17:58:56 2024

@author: james
"""
#Import packages
import pandas as pd
import matplotlib.pyplot as plt
import subprocess
import re
import shutil
import os


#Run Smoldyn from command line
def run_smoldyn(input_file, output_file_smoldyn):
    
    smoldyn_command = ['smoldyn', input_file]
    #process = subprocess.Popen(smoldyn_command, stdout=subprocess.PIPE, stderr = subprocess.PIPE)
    subprocess.run(smoldyn_command)
    #stdout, stderr = process.communicate()
    
    
    # Decode the output to string
    #output_text = stdout.decode('utf-8')
    smoldyn_output = output_file_smoldyn
    folder = "Output3/"
    shutil.move(smoldyn_output, folder)
   # with open(output_file_path, 'r') as file:
        #smoldyn_output = file.read()
    
    # Do something with the output
    return(f"Output3/{smoldyn_output}")
        
   

    

#Plotting graphs with matplot.lib
##Plot total
def plot_everything(output_file):
    #Get output file from smoldyn
    file = pd.read_csv(output_file, sep =' ')
    df = pd.DataFrame(file)
    
    #Isolate time points
    x = df.iloc[:, 0]
    
    #Define variable names
    variables = ["TnB","TcmB","TemB","TeffB","CD19B","Comp_TnB","Comp_TcmB","Comp_TemB","Comp_TeffB","TnP","TcmP","TemP","TeffP","CD19P","Comp_TnP","Comp_TcmP","Comp_TemP","Comp_TeffP","CD19_micro","Comp_Tn_micro","Comp_Tcm_micro","Comp_Tem_micro","Comp_Teff_micro"]
    fig, ax = plt.subplots()
    
    for i in range(1, len(variables)):
        y = df.iloc[:, i]
        label = variables[i-1]
        plt.plot(x, y, label = label)
    
    plt.xlabel("Time Points (days)")
    plt.ylabel("Abundance")
    plt.title("Smoldyn Simulation of CAR-T Therapy on ALL")
    
    
    ax.legend()
    
    plt.show()
    

def get_best_dose(output_path):
    file = pd.read_csv(output_path, sep =' ')
    df = pd.DataFrame(file)
    tumour_blood = df.iloc[:, 5]
    tumour_bm = df.iloc[:, 14] + df.iloc[:, 19]
    tumour_total = tumour_blood + tumour_bm
    min_value = min(tumour_total)
    return(min_value)





##Plot by compartment
def plot_tumour(output_path):
    #Get output file from smoldyn
    file = pd.read_csv(output_path, sep =' ')
    df = pd.DataFrame(file)
    
    #Isolate time points
    x = df.iloc[:, 0]
    
    #Isolate tumour changes
    tumour_blood = df.iloc[:, 5]
    tumour_bm = df.iloc[:, 14] + df.iloc[:, 19]
    tumour_total = tumour_blood + tumour_bm
    
    fig, ax = plt.subplots()
    plt.plot(x, tumour_blood, label = "Tumour Blood")
    plt.title("Change in tumour abundance during CAR-T therapy")
    # plt.xlabel("Time Points (days)")
    # plt.ylabel("Abundance")
    # plt.show()
    
    plt.plot(x, tumour_bm, label = "Tumour Bone Marrow")
    # plt.title("Change in tumour abundance in the bone marrow during CAR-T therapy")
    # plt.xlabel("Time Points (days)")
    # plt.ylabel("Abundance")
    # plt.show()
    
    plt.plot(x, tumour_total, label = "Tumour total")
    # plt.title("Total change in tumour abundance during CAR-T therapy")
    plt.xlabel("Time Points (days)")
    plt.ylabel("Abundance")
    plt.legend(loc = 'center right')
    filename = os.path.basename(output_path)
    fig.savefig(f"Plots3/{filename}.png", dpi = 100)
    

def plot_tcell_total(output_file):
    #Get output file from smoldyn
    file = pd.read_csv(output_file, sep =' ')
    df = pd.DataFrame(file)
    
    #Isolate time points
    x = df.iloc[:, 0]
    
    tcell_blood = df.iloc[:,1] + df.iloc[:,2] + df.iloc[:,3] + df.iloc[:,4]
    tcell_bm = df.iloc[:,11] + df.iloc[:,12] + df.iloc[:,13] + df.iloc[:,11]
    tcell_total = tcell_blood + tcell_bm
    
    plt.plot(x, tcell_blood)
    plt.title("Change in T-Cell abundance in the blood during CAR-T therapy")
    plt.xlabel("Time Points (days)")
    plt.ylabel("Abundance")
    plt.show()
    
    plt.plot(x, tcell_bm)
    plt.title("Change in T-Cell abundance in the bone marrow during CAR-T therapy")
    plt.xlabel("Time Points (days)")
    plt.ylabel("Abundance")
    plt.show()
    
    plt.plot(x, tcell_total)
    plt.title("Total change in T-Cell abundance during CAR-T therapy")
    plt.xlabel("Time Points (days)")
    plt.ylabel("Abundance")
    plt.show()
    
    
def plot_tcell_type(output_file):
    #Get output file from smoldyn
    file = pd.read_csv(output_file, sep =' ')
    df = pd.DataFrame(file)
    
    #Isolate time points
    x = df.iloc[:, 0]
    
    tcell_naive_blood = df.iloc[:,1]
    tcell_cm_blood = df.iloc[:,2]
    tcell_em_blood = df.iloc[:,3]
    tcell_eff_blood = df.iloc[:,3]
    tcell_naive_bm = df.iloc[:,3]
    tcell_cm_bm = df.iloc[:,3]
    tcell_em_bm = df.iloc[:,3]
    tcell_eff_bm = df.iloc[:,3]
    
    
    plt.plot(x, tcell_naive_blood, color = 'red', label = 'TnB' )
    plt.plot(x, tcell_naive_bm, color = 'blue', label = 'TnBM')
    plt.legend(loc = 'upper right')
    plt.title("Change in Naive T-Cell abundance during CAR-T therapy")
    plt.xlabel("Time Points (days)")
    plt.ylabel("Abundance")
    plt.show()
    
   
    plt.plot(x, tcell_cm_blood, color = 'red', label = 'TcmB' )
    plt.plot(x, tcell_cm_bm, color = 'blue', label = 'TcmBM')
    plt.legend(loc = 'upper right')
    plt.title("Change in Central Memory T-Cell abundance  during CAR-T therapy")
    plt.xlabel("Time Points (days)")
    plt.ylabel("Abundance")
    plt.show()
    
    plt.plot(x, tcell_em_blood, color = 'red', label = 'TemB' )
    plt.plot(x, tcell_em_bm, color = 'blue', label = 'TemBM')
    plt.legend(loc = 'upper right')
    plt.title("Change in Effector Memory T-Cell abundance during CAR-T therapy")
    plt.xlabel("Time Points (days)")
    plt.ylabel("Abundance")
    plt.show()
    
    plt.plot(x, tcell_eff_blood, color = 'red', label = 'TeffB' )
    plt.plot(x, tcell_eff_bm, color = 'blue', label = 'TeffBM')
    plt.legend(loc = 'upper right')
    plt.title("Change in Effector T-Cell abundance during CAR-T therapy")
    plt.xlabel("Time Points (days)")
    plt.ylabel("Abundance")
    plt.show()

# start_value = 0.000001
# end_value = 1000
geometric_sequence = [250]
min_tumour_values = []
# current_value = start_value

# Keep multiplying by 10 until the value reaches or exceeds the end value
# while current_value <= end_value:
#     geometric_sequence.append(current_value)
#     current_value *= 10

#Loop to alter K_ON
for i in range(0,2):
    input_file = "3D_vr_tumour_ball.txt"
    with open(input_file, 'r') as file:
        lines = file.readlines()
    
    for z,line in enumerate(lines):
        if "define K_ON" in line:
            lines[z] = "define K_ON " + str(200) +"\n"
    #Write back changes to file
    with open(input_file, 'w') as file:
        file.writelines(lines)
        
    #For each parameter try a different dose
    for y in range(0,100,12):
        output_file_name = "out_Teff" + str(y*0.25) + "_KON200_rep"+str(i+1)+".txt"
        num_Tn = 100-(0.75*y)
        num_Tcm = 0.25*y
        num_Tem = 0.25*y
        num_Teff = 0.25*y
        with open(input_file, 'r') as file:
            lines = file.readlines()
        
        for x,line in enumerate(lines):
            if "compartment_mol" and "TnB blood" in line:
                lines[x] = "compartment_mol " + str(num_Tn) + " TnB blood\n"
            
            elif "compartment_mol" and "TcmB blood" in line:
                lines[x] = "compartment_mol " + str(num_Tcm) + " TcmB blood\n"
            
            elif "compartment_mol" and "TemB blood" in line:
                lines[x] = "compartment_mol " + str(num_Tem) + " TemB blood\n"
            
            elif "compartment_mol" and "TeffB blood" in line:
                lines[x] = "compartment_mol " + str(num_Teff) + " TeffB blood\n"
            
            #Finally need to change save of file
            elif "output_files" in line:
                lines[x] = "output_files " + output_file_name +"\n"
            
            elif "cmd i 0 1000 0.01 molcount" in line:
                lines[x] = "cmd i 0 1000 0.01 molcount " + output_file_name +"\n"
                
        #Write back changes to file
        with open(input_file, 'w') as file:
            file.writelines(lines)
        
        smoldyn_out = run_smoldyn(input_file, output_file_name)
        plot_tumour(smoldyn_out)
        min_tumour_values.append(get_best_dose(smoldyn_out))
       
        


        
    
    

