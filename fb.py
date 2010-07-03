# print flour beetle data in R readable format
import numpy as np
logdose = (-1.61, -1.14, -0.69, -0.22)
data=\
[3  ,7  ,5  ,4,
11 ,10 ,8  ,10,
10 ,11 ,11 ,8 ,
7  ,16 ,15 ,14,
4  ,3  ,4  ,8 ,
3  ,2  ,2  ,2 ,
2  ,1  ,1  ,1 ,
1  ,0  ,1  ,0 ,
0  ,0  ,0  ,0 ,
0  ,0  ,0  ,1 ,
0  ,0  ,0  ,0 ,
1  ,0  ,0  ,0 ,
1  ,0  ,0  ,0 ,
101,19,7   ,2]

in_data = np.array(data).reshape((-1,4))

print 'ld t.l t.u ic rc'
# 1,x,tl,tu,ic,rc,
for d in xrange(len(logdose)):
    x = logdose[d]
    for t in xrange(len(in_data)-1):
        # interv censored
        for r in xrange(in_data[t,d]):
            row = [x] + [t, t+1] + [1] + [0]
            print ' '.join(map(str,row))
        # right censored
    for r in xrange(in_data[-1,d]):
        row = [x] + [len(in_data)-1, -1] + [0] + [1]
        print ' '.join(map(str,row))

