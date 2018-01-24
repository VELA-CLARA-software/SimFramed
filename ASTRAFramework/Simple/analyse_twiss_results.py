import csv

results = []

with open('twiss_best_2013/best_solutions.csv', 'r') as csvfile:
  reader = csv.reader(csvfile, delimiter=',', quoting=csv.QUOTE_NONNUMERIC)
  for row in reader:
    results.append(row)

astra = ASTRAInjector('', overwrite=False)
astra.loadSettings('short_240.settings')
astra.fileSettings['test.2']['quad_K'] = results[0][0:6]
astra.fileSettings['test.3']['quad_K'] = results[0][6:14]
astra.fileSettings['test.4']['quad_K'] = results[0][14:16]
astra.fileSettings['test.5']['quad_K'] = results[0][16:]

print astra.settings
