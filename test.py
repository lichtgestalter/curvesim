import math

minpixel = 0.005
maxpixel = 0.05
rlist = [240000, 6300, 25000, 700000]
print(f'{rlist=}')
log_list = [math.log10(i) for i in rlist]
print(f'{log_list=}')
print(log_list)
min_ok = [i * minpixel / min(log_list) for i in log_list]
print(f'{min_ok=}')
x = math.log10((maxpixel-max(min_ok))/max(min_ok))
y = math.log10(max(log_list)-min(log_list))
exponent = x / y
print(x, y, exponent)
fertig = [i * (1 + (j - min(log_list)) ** exponent) for i, j in zip(min_ok, log_list)]  # H7*(1+(D7-$D$4)^$I$1)
print(f'{fertig=}')