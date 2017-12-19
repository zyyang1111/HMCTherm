import csv
import numpy as np
import sys
import matplotlib.pyplot as plt

file_name = sys.argv[1]
S_id = int(sys.argv[2])
results = []

with open(file_name) as csvFile:
	reader = csv.reader(csvFile)
	next(reader, None)
	for row in reader: 
		results.append(row)

S_n = 0
L_n = 0;
X_D = 0;
Y_D = 0;
for row in results:
	if int(row[0]) > S_n:
		S_n = int(row[0])

	if int(row[1]) > L_n:
		L_n = int(row[1])

	if int(row[2]) > X_D:
		X_D = int(row[2])

	if int(row[3]) > Y_D:
		Y_D = int(row[3])

print S_n, L_n, X_D, Y_D


L_size = (X_D+1)*(Y_D+1)
S_size = (L_n+1)*L_size

if S_id >= S_size:
	S_id = S_size-1


I = np.zeros(((X_D+1, Y_D+1, L_n+1)))
for i in range(0, S_size):
	z = int(results[i+S_id*S_size][1])
	x = int(results[i+S_id*S_size][2])
	y = int(results[i+S_id*S_size][3])
	p = float(results[i+S_id*S_size][4])

	I[x][y][z] = p

for layer in range(0, L_n+1):
	plt.figure()
	plt.imshow(I[:,:,layer], aspect='auto')
	plt.show()

