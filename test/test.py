import numpy as np
import math as m

x=[0.82,0.34,0.69,0.44,0.97,0.10,0.23,0.78,0.92,0.81]
y=0

for i in range(10):
    y += m.sin(x[i])/m.log(x[i]+2)

print(y/10)
