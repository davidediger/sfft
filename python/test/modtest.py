import sfft.sfft as sfft
import numpy as np

NOISE_THRESHOLD = 0.1
ERROR_THRESHOLD = 0.1

n = 16384
k = 50
freq = np.random.random_integers(0, n-1, k)
a = np.zeros(n, dtype=complex)
for x in freq:
  a[x] = 1.0
b = np.fft.ifft(a)
b = b * n
myfft = sfft.sfft(n, k, 1)
result = myfft.execute(b)

for x in xrange(n):
  if np.absolute(a[x]) > NOISE_THRESHOLD:
    if np.absolute(result[x]) == 0.0:
      print "Error: frequency not recovered"
      print x
      print result[x]
      exit()
    if np.absolute(result[x] - a[x]) > ERROR_THRESHOLD:
      print "Error: freqency error is too big"
      print x
      print "Expected: " + str(a[x])
      print "Actual: " + str(result[x])
      exit()

print "OK"

