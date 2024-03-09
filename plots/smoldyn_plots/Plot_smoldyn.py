from Smoldyn_pipeline2 import plot_tumour, get_best_dose
import pandas as pd
import matplotlib.pyplot as plt
import os
import math
import statistics

file_path = "Final_Output/G6_mini-project/Smoldyn_folder/Final/Output"
output_files = os.listdir(file_path)
output_files.sort()
file_rank = []

for file in output_files:
    file_rank.append(get_best_dose(f"{file_path}/{file}"))
    




print(len(file_rank))
#Order important as data_points done alphabetically
dose = [100, 64, 55, 46, 37, 28, 91, 82, 73]
dose = [100,91,82,73,64,55,46,37,28]


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

#Reorder scores to match dose
score_reorder = []
score_reorder.append(score[0])
score_reorder.append(score[6:10])
score_reorder.append(score[1:6])

average_reorder = []
average_reorder.append(average[0])
average_reorder.append(average[6:10])
average_reorder.append(average[1:6])

error_reorder = []
error_reorder.append(error[0])
error_reorder.append(error[6:10])
error_reorder.append(error[1:6])


def flatten_list(nested_list):
    flat_list = []
    for element in nested_list:
        if isinstance(element, list):
            flat_list.extend(flatten_list(element))
        else:
            flat_list.append(element)
    return flat_list

score_reorder = flatten_list(score_reorder)
average_reorder = flatten_list(average_reorder)
error_reorder = flatten_list(error_reorder)
print(average_reorder)
print(error_reorder)


fig, ax = plt.subplots(figsize = (12, 8))
plt.plot(dose, score_reorder, color = 'black')
plt.errorbar(dose, score_reorder, yerr = error_reorder, fmt = 'o', capsize =3, color = "black", elinewidth = 0.5)
ax.set_xlabel("Dose Strategy (number of naive T cells)", fontsize = 15, labelpad = 10)
ax.set_ylabel("Efficacy Score", fontsize = 15, labelpad = 10)
plt.title("Naïve CAR-T Cell Number Determines Therapeutic Efficacy", fontsize = 20, fontstyle = 'oblique', pad = 10)
plt.grid()
ax.set_ylim(0,1)
#ax.set_xlim(0,120)
plt.show()
fig.savefig("Final_Output/G6_mini-project/Smoldyn_folder/Final/Plots/Simulationfig.png", dpi = 100)




#Missing rep 3 and 4 from comp 5 - currently 18 repeats
# search = "rep4_Comp4"
# found = [string for string in output_files if search in string]
# for string in found:
#     print(string)

    


