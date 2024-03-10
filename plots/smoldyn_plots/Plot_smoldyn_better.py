from Smoldyn_pipeline2 import get_best_dose
import matplotlib.pyplot as plt
import os
import math
import statistics

file_path = "Output"
output_files = os.listdir(file_path)
file_rank = []

def extract_number(file_name):
    return float(file_name.split("out_Teff")[1].split("_")[0])

sorted_output_files = sorted(output_files, key = extract_number)

for file in sorted_output_files:
    file_rank.append(get_best_dose(f"{file_path}/{file}"))


print(len(file_rank))
#Order important 
dose = [100,91,82,73,64,55,46,37,28, 19, 10, 1]


#Score the doses
score = []
average =[]
deviation = []
error = []
num_iterations = len(file_rank) //20 
for i in range(num_iterations):
    start_index = i*20
    end_index = (i+1)*20

    sublist = file_rank[start_index:end_index]
    sublist_average = sum(sublist)/len(sublist)
    sublist_deviation = statistics.stdev(sublist)
    

    average.append(sublist_average)
    deviation.append(sublist_deviation)

#Calculate and transform standard error
for i in deviation:
    sub_error = i/(1500*math.sqrt(20))
    error.append(sub_error)


#Calculate the score
for i in average:
    item_score = (1500- i)/(1500)
    score.append(item_score)

print(score)

    



#Plot Smoldyn simulation results
fig, ax = plt.subplots(figsize = (12, 8))
plt.plot(dose, score, color = 'black')
plt.errorbar(dose, score, yerr = error, fmt = 'o', capsize =3, color = "black", elinewidth = 0.5)
ax.set_xlabel("Dose Strategy (number of naive T cells)", fontsize = 15, labelpad = 10)
ax.set_ylabel("Efficacy Score", fontsize = 15, labelpad = 10)
plt.title("Naïve CAR-T Cell Number Influences Therapeutic Efficacy", fontsize = 20, fontstyle = 'oblique', pad = 10)
plt.grid()
ax.set_ylim(0,1)
plt.show()
fig.savefig("Plots/Simulationfig.png", dpi = 100)





    


