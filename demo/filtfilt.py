import sys
from numpy import genfromtxt
from matplotlib.pyplot import figure, plot, legend
from scipy.signal import filtfilt, firwin
from datetime import datetime

delim = None
if(len(sys.argv)>2):
    delim = sys.argv[2]

column = 0
if(len(sys.argv)>3):
    column = int(sys.argv[3])

ppg_data = genfromtxt(sys.argv[1], delimiter=delim)
roi = None
if delim is None:
    roi = ppg_data[0:-1]
else:
    roi = ppg_data[:,column]

b = genfromtxt("./bCoeff")
a = genfromtxt("./aCoeff")
dt_start = datetime.now()
filtered = filtfilt(b, a, roi)
dt_end = datetime.now()

print("Python filtfilt took ", (dt_end.microsecond - dt_start.microsecond)/1000, "ms")
demoFiltFilt = genfromtxt("./demoFiltFilt")[:]
fig = figure(figsize=(24,4))
plot(roi, label="Original")
plot(filtered, label="Python Filtered")
plot(demoFiltFilt, label="C++ Filtered")
legend()
fig.savefig('../filtfilt.pdf')
print("Plot was saved to ./filtfilt.pdf")