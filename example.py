from pvhub import*
#-------------------------------------------------------------------------------------#
f0 = open("./inputs/example.csv","r")
#--------------------------------------------------------------------------------------#
for line in f0:
    if line[0]!='#':
        SNID = str(line.split(",")[0])
        folder01 = str(line.split(",")[1])
        ra00 = float(line.split(",")[5])
        dec00 = float(line.split(",")[6])
        z00 = float(line.split(",")[9])
        pv = calculate_pv(SNID,folder01,ra00,dec00,z00,extrapolation="Yes")
        print (SNID,pv)
f0.close()
#---------------------------------------------------------------------------------------#
