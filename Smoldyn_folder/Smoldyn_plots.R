setwd("~/Desktop")
library(ggplot2)
analysis_df = as.data.frame.matrix(read.table("out_tumour_ball.txt", header = FALSE))

time_points = analysis_df$V1
#Blood compartment variables
TnB = analysis_df$V2
TcmB = analysis_df$V3
TemB = analysis_df$V4
TeffB = analysis_df$V5
CD19B = analysis_df$V6
Comp_TnB = analysis_df$V7
Comp_TcmB = analysis_df$V8
Comp_TemB = analysis_df$V9
Comp_TeffB = analysis_df$V10

#Peripheral tissue variables

TcmP = analysis_df$V12
TemP = analysis_df$V13
TeffP = analysis_df$V14
CD19P = analysis_df$V15
Comp_TnP = analysis_df$V16
Comp_TcmP = analysis_df$V17
Comp_TemP = analysis_df$V18
Comp_TeffP = analysis_df$V19

#Tumour ball variables
CD19_micro = analysis_df$V20
Comp_Tn_micro = analysis_df$V21
Comp_Tcm_micro = analysis_df$V22
Comp_Tem_micro = analysis_df$V23
Comp_Teff_micro = analysis_df$V24

#Totals
Comp_Tn = Comp_Tn_micro + Comp_TnP + Comp_TnB
Comp_Tcm = Comp_Tcm_micro + Comp_TcmP + Comp_TcmB
Comp_Tem = Comp_Tem_micro + Comp_TemP + Comp_TemB
Comp_Teff = Comp_Teff_micro + Comp_TeffP + Comp_TeffB

CD19 = CD19B + CD19P + CD19_micro

Tn = TnP + TnB
Tcm = TcmP + TcmB
Tem = TemP + TemB
Teff = TeffP + TeffB

#Plot Total 
plot(time_points, TnB, main = "Blood compartment", type = 'l', col = "red", xlab = "Time")
plot(time_points, TcmB, col = "orange", type = 'l')
plot(time_points, TemB, col = "blue", type = 'l')
plot(time_points, TeffB, col = "purple", type = 'l')
lines(time_points, CD19B, col = "black", type = 'l')

plot(time_points, CD19P+CD19_micro, main = "Bone Marrow Compartment", type = 'l', col = "black", xlab = "Time", ylab = "Cell Number")
lines(time_points, TcmP, col = "orange", type = 'l')
lines(time_points, TemP, col = "blue", type = 'l')
lines(time_points, TeffP, col = "purple", type = 'l')
lines(time_points, TnP, col = "red", type = 'l')

