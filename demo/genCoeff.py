import sys
from scipy.signal import firwin

firf = firwin(61, [0.66, 4.0], pass_zero=False, fs=30.0)
data = ""
for val in firf: 
        data += str(val) + "\n"
f = open("./bCoeff", "w")
f.write(data)
f.close()

f = open("./aCoeff", "w")
f.write("1")
f.close()