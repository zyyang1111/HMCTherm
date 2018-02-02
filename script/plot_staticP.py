import csv
import numpy as np
import sys
import matplotlib.pyplot as plt

file_name = sys.argv[1]
results = []

with open(file_name) as csvFile:
	reader = csv.reader(csvFile)
	next(reader, None)
	for row in reader: 
		results.append(row)

L_n = 0;
X_D = 0;
Y_D = 0;

for row in results:
	if int(row[1]) > L_n:
		L_n = int(row[1])

	if int(row[2]) > X_D:
		X_D = int(row[2])

	if int(row[3]) > Y_D:
		Y_D = int(row[3])

print L_n, X_D, Y_D

T = np.zeros(((X_D+1, Y_D+1, L_n+1)))
for i in range(0, len(results)):
	z = int(results[i][1])
	x = int(results[i][2])
	y = int(results[i][3])
	p = float(results[i][4])

	T[x][y][z] = p

for layer in range(0, L_n+1):
	plt.figure()
	plt.imshow(T[:,:,layer], aspect='auto')
	plt.show()

