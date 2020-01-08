import numpy as np
import time

v1 = np.random.rand(1000000)
v2 = np.random.rand(1000000)
v3 = np.empty(1000000)

t = time.clock()
for i in range(0,len(v1)-1) :
    v3[i] = v1[i] * v2[i]

print("Time with using for : " + str(time.clock()- t))

t = time.clock()
v3 = v1.dot(v2)
print("Time with using .dot : " + str(time.clock()- t))
