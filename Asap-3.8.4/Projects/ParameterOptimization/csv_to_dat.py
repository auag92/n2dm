"""Convert Excel CSV file to semi-colon separated - and fix the line endings.

Excel on a Mac cannot export a proper semi-colon separated file.
"""

import csv

infile = "properties_metals.csv"
outfile = "properties_metals.dat"

of = open(outfile, 'w')
outwriter = csv.writer(of, delimiter=';')

for line in csv.reader(open(infile, 'rU')):
    outwriter.writerow(line)
    
    
