import csv
from ASTRAInjector import *

results = []

with open('Short_210MeV_6.7MeVps/best_solutions.csv', 'r') as csvfile:
  reader = csv.reader(csvfile, delimiter=',', quoting=csv.QUOTE_NONNUMERIC)
  for row in reader:
    results.append(row)

astra = ASTRAInjector('', overwrite=False)
astra.loadSettings('short_240_12b3.settings')
linac1field, linac1phase, linac2field, linac2phase, linac3field, linac3phase, fhcfield, fhcphase, linac4field, linac4phase = results[0]
astra.modifySetting('linac1_field', abs(linac1field))
astra.modifySetting('linac1_phase', linac1phase)
astra.modifySetting('linac2_field', abs(linac2field))
astra.modifySetting('linac2_phase', linac2phase)
astra.modifySetting('linac3_field', abs(linac3field))
astra.modifySetting('linac3_phase', linac3phase)
astra.modifySetting('4hc_field', abs(fhcfield))
astra.modifySetting('4hc_phase', fhcphase)
astra.modifySetting('linac4_field', abs(linac4field))
astra.modifySetting('linac4_phase', linac4phase)
astra.saveSettings(filename='short_240_12b3.settings')

astra.loadSettings('short_240_12b3.settings')
print astra.settings
