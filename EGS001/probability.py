import math

m = 10 ^ 7
N = 10 ^ 13
sum = 0
drawing = 0

for k in range(20):
    choose = math.comb(m, k)
    num = (N - 1) ^ (m-k)
    den = N ^ (m-1)
    result = choose * num / den
    sum = sum + result
    drawing = drawing + result * k
    print("TCRs drawn ", k, " times: ", result)
    print("Sum unique: ", sum)
    print("Sum drawings: ", drawing)

k = 2
m = 10000
choose = math.comb(m, k)
num = (N - 1) ^ (m-k)
den = N ^ (m-1)
result = choose * num / den
sum = sum + result
drawing = drawing + result * k
print("TCRs drawn ", k, " times: ", result)
print("Sum unique: ", sum)
print("Sum drawings: ", drawing)
