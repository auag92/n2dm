"""Get EMT2013 parameters from fitting script output file."""

from ase.atoms import string2symbols

def get_parameters(file, number=2):
    file = open(file)
    text = file.read()
    file.close()

    # Find elements
    s = -1
    for i in range(number):
        s = text.find('Optimization', s + 1)
        if s < 0:
            raise ValueError("No results in file (keyword Optimization missing in file)")
    e = text.find('\n', s)
    errfunc = float(text[s:e].split()[6])
    result = {}

    s = e + 1
    e = text.find('Fitting values', s) - 4
    for line in text[s:e].split('\n'):
        words = line.strip().split()
        if words[1] == 'parameters':
            elements = tuple(string2symbols(words[0]))
            result[elements] = []
        else:
            result[elements].append(float(words[1]))
    
    return result, errfunc
