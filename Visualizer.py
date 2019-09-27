import seaborn as sb
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np


def vis_best_evalues(ordered_list):
    df = pd.DataFrame(ordered_list, columns=['Descriptor', 'Best Hit', 'Distance Factor', 'Actual Distance',
                                             'Distance Factor - Best Domain', 'Actual Distance - Best Domain',
                                             'Category'])
    df.to_csv('/sybig/home/ttu/BSc/dfout.txt', header=None, index=None, sep='\t', mode='w')
    sb.set(style="ticks")
    sb.relplot(x='Distance Factor', y='Actual Distance', hue='Category', data=df)
    sb.relplot(x='Distance Factor - Best Domain', y='Actual Distance - Best Domain', hue='Category', data=df)
    plt.show()



# palette=['#51f542', '#42aaf5', '#4290f5', '#4266f5', '#302bad', '#f54242']