import csv
import matplotlib.pyplot as plt

with open('output_agent_ADMM.csv') as csvfile:
    readCSV = csv.reader(csvfile, delimiter=',')

    g_draw_seq_no = []
    train_id = []
    line_id = []
    station_id = []
    station_name = []
    x_time_in_min = []
    y_station_location = []
    
    for row in readCSV:
        g_draw_seq_no.append(row[0])
        train_id.append(row[1])
        line_id.append(row[2])
        station_id.append(row[3])
        station_name.append(row[4])
        x_time_in_min.append(row[5])
        y_station_location.append(row[6])
    
#remove the titles
g_draw_seq_no.remove(g_draw_seq_no[0])
train_id.remove(train_id[0])
line_id.remove(line_id[0])
station_id.remove(station_id[0])
station_name.remove(station_name[0])
x_time_in_min.remove(x_time_in_min[0])
y_station_location.remove(y_station_location[0])

color_value = {
        '1': 'midnightblue', 
        '2': 'mediumblue', 
        '3':'c',
        '4':'orangered',
        '5':'m',
        '6':'fuchsia',
        '7':'olive'
    }

xlist = []
ylist = []

for i in range(len(g_draw_seq_no)):
    next_line_no = min(i + 1, len(g_draw_seq_no) - 1)
    if train_id[i] == train_id[next_line_no]: #for the current train
        if g_draw_seq_no[i] == g_draw_seq_no[next_line_no]:
            if next_line_no == len(g_draw_seq_no) - 1:
                xlist.append(int(x_time_in_min[i]))
                ylist.append(int(y_station_location[i]))
                plt.plot(xlist, ylist, color = color_value[str(line_id[i])], linewidth = 1.5)
                plt.text(xlist[0] + 0.8, ylist[0] + 4, '%d' % int(line_id[i]), ha='center', va= 'bottom', color = color_value[str(line_id[i])], weight = 'bold', family = 'Times new roman', fontsize= 9)
            else:
                xlist.append(int(x_time_in_min[i]))
                ylist.append(int(y_station_location[i]))
        else:
            xlist.append(int(x_time_in_min[i]))
            ylist.append(int(y_station_location[i]))
            plt.plot(xlist, ylist, color = color_value[str(line_id[i])], linewidth = 1.5)
            plt.text(xlist[0] + 0.8, ylist[0] + 4, '%d' % int(line_id[i]), ha='center', va= 'bottom', color = color_value[str(line_id[i])], weight = 'bold', family = 'Times new roman', fontsize= 9)
            xlist = []
            ylist = []
    else:
        xlist.append(int(x_time_in_min[i]))
        ylist.append(int(y_station_location[i]))
        plt.plot(xlist, ylist, color = color_value[str(line_id[i])], linewidth = 1.5)
        plt.text(xlist[0] + 0.8, ylist[0] + 4, '%d' % int(line_id[i]), ha='center', va= 'bottom', color = color_value[str(line_id[i])], weight = 'bold', family = 'Times new roman', fontsize= 9)
        xlist = []
        ylist = []

plt.grid(True) #show the grid
plt.ylim(0, 360)  # y range

plt.xlim(0, 90)  # x range
plt.xticks([0, 10, 20, 30, 40, 50, 60, 70, 80, 90])

plt.yticks([0, 85, 135, 195, 270, 360], ['S1', 'S2', 'S3', 'S4', 'S5', 'S6'], family = 'Times new roman')
plt.xlabel('Time (min)', family = 'Times new roman')
plt.ylabel('Space (km)', family = 'Times new roman')
plt.show()