import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
input_file_address = input("Input the input address or filename:\n")
data_xy = pd.read_csv(input_file_address, header=None, index_col=0, delimiter="\s+")
output_file_address = input("Input the output file address or filename:\n")
data_order = pd.read_csv(output_file_address)

data_xy.columns = ["x_coord","y_coord"]
#print(data_xy.iloc[0])

data_order.columns = ["order"]
data_xy_aftersorted = data_xy.loc[data_order["order"]]
data_xy_aftersorted = data_xy_aftersorted.dropna()
data_xy_aftersorted = pd.concat([data_xy_aftersorted,data_xy_aftersorted.iloc[[0]]],ignore_index=True)
#print(data_xy_aftersorted,sep="\n\n")

graph_name = input("Input the graph name:\n")
data_xy_aftersorted.plot(x=data_xy_aftersorted.columns[0],y=data_xy_aftersorted.columns[1],kind="line", title=graph_name, xlabel=data_xy_aftersorted.columns[0],ylabel=data_xy_aftersorted.columns[1], legend=False)
plt.show()

