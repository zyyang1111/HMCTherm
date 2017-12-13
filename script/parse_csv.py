import csv

maxT = 0; 

with open('temperature_trace.csv') as csvDataFile:
	csvReader = csv.reader(csvDataFile)
	cnt = 0; 
	for row in csvReader:
		if cnt > 0:
			ID = row[0]
			layer = row[1]
			x = row[2]
			y = row[3]
			data = float(row[4])
			if data > maxT:
				maxT = data


		cnt = cnt + 1

print maxT