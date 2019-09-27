import matplotlib.pyplot as plt
import numpy as np

objects = ('Category 1', 'Category 2', 'Category 3', 'Category 4', 'Category 5')
y_pos = np.arange(len(objects))
performance = [0.09115159826070314, 0.25663133788587067, 0.11329672866539017, 0.1424284388300628, 0.044612605470126236]

plt.bar(y_pos, performance, align='center', alpha=0.5, color=['chartreuse', 'gold', 'darkorange', 'red', 'gray'])
plt.xticks(y_pos, objects)
plt.ylabel('Average entropy')
plt.title('Average entropy for each category')

plt.show()
