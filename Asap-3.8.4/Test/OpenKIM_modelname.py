import os

openkimmodel = "EMT_Asap_Standard_Jacobsen_Stoltze_Norskov_AlAgAuCuNiPdPt__MO_118428466217_002"

if 'ASAP_KIM_DIR' in os.environ:
    _d = os.path.join(os.environ['ASAP_KIM_DIR'], 'src/models')
elif 'ASAP_KIM_LIB' in os.environ:
    _d  = os.path.join(os.environ['ASAP_KIM_LIB'], 'kim-api-v1/models')
elif 'KIM_HOME' in os.environ:
    _d  = os.path.join(os.environ['KIM_HOME'], 'lib/kim-api-v1/models')
else:
    _d = None

if _d is not None:
    openkimmodels = [x for x in os.listdir(_d) if os.path.isdir(os.path.join(_d,x))]
    openkimmodels.sort()
else:
    openkimmodels = []
    
if __name__ == "__main__":
    print openkimmodels
    print "\n\nThis it not a test, but a module imported from a few tests."
    
    
