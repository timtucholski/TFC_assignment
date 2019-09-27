import matplotlib.pyplot as plt

labels = 'Correct genus', 'Wrong genus'
sizes = [1540, 68]
explode = (0, 0.15)

fig1, ax1 = plt.subplots()
ax1.pie(sizes, explode=explode, labels=labels, autopct='%1.1f%%',
        shadow=False, startangle=90, colors=['limegreen', 'red'], pctdistance=1.05, labeldistance=1.3)
ax1.axis('equal')  # Equal aspect ratio ensures that pie is drawn as a circle.

plt.show()
