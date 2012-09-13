import pylab as pl
import numpy as np
import numpy.fft as fft

def plot(x, y, filename = False, clear = True):
    if clear:
        pl.clf()
    pl.plot (x, y, "-")
    if filename:
        pl.savefig (filename)
    else:
        pl.show()

N = 1024
x = np.arange(N)
f = fft.fftfreq (N)

# Una delta
y = np.zeros(N)
y[N/2] = 1
plot (x, y, "delta")

Y = abs(fft.fft(y))
plot (f, Y, "f_delta")

# 2 deltas a)
y[N/2-1] = 1
plot (x, y, "2deltas_a")

Y = abs(fft.fft(y))
plot (f, Y, "f_2deltas_a")

# 2 deltas b)
y = np.zeros(N)
y[N/2] = 1
y[N/2-2] = 1
plot (x, y, "2deltas_b")

Y = abs(fft.fft(y))
plot (f, Y, "f_2deltas_b")

# 2 deltas b)
y = np.zeros(N)
y[N/2] = 1
y[N/2-8] = 1
plot (x, y, "2deltas_c")

Y = abs(fft.fft(y))
plot (f, Y, "f_2deltas_c")
