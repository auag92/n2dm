import os

for root, dirs, files in os.walk("."):
    for file in files:
        suffix = os.path.splitext(file)[1]
        for cand in ['.py', '.swig', '.c', '.cpp', '.h']:
            if suffix == cand:
                txt = open(os.path.join(root, file)).read()
                if txt:
                    mx = ord(max(txt))
                else:
                    mx = 0
                if mx > 127:
                    print os.path.join(root, file), "contains character", mx

