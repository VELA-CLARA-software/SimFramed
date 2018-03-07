import csv
from ASTRAInjector import *

def between(value, minvalue, maxvalue, absolute=True):
    if absolute:
        result = max([minvalue,min([maxvalue,abs(value)])])
    else:
        result = np.sign(value)*max([minvalue,min([maxvalue,abs(value)])])
    return float(result)


results = []

with open('longitudinal_best/longitudinal_best_solutions.csv', 'r') as csvfile:
  reader = csv.reader(csvfile, delimiter=',', quoting=csv.QUOTE_NONNUMERIC)
  for row in reader:
    results.append(row)

astra = ASTRAInjector('', overwrite=False)
astra.loadSettings('short_240_12b3.settings')
linac1field, linac1phase, linac2field, linac2phase, linac3field, linac3phase, fhcfield, fhcphase, linac4field, linac4phase, bcangle = results[0]
astra.modifySetting('linac1_field', between(linac1field,10,32))
astra.modifySetting('linac1_phase', between(linac1phase,0,45, False))
astra.modifySetting('linac2_field', between(linac2field,10,32))
astra.modifySetting('linac2_phase', between(linac2phase,0,45, False))
astra.modifySetting('linac3_field', between(linac3field,10,32))
astra.modifySetting('linac3_phase', between(linac3phase,0,45, False))
astra.modifySetting('4hc_field', between(fhcfield,10,35))
astra.modifySetting('4hc_phase', between(fhcphase,170,200))
astra.modifySetting('linac4_field', between(linac4field,10,32))
astra.modifySetting('linac4_phase', between(linac4phase,0,45, False))
astra.fileSettings['test.4']['variable_bunch_compressor']['angle'] = between(bcangle,0.05,0.12,True)
astra.saveSettings(filename='short_240_12b3.settings')
astra.loadSettings('short_240_12b3.settings')
print astra.settings
