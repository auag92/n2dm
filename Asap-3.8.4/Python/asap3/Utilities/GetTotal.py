import pickle
import Scientific.MPI
from Numeric import *

nProcs = Scientific.MPI.world.size
ismaster = Scientific.MPI.world.rank == 0

if nProcs > 1:
    parallel = 1
else:
    parallel = 0


def GetTotalScalar(value):
    
    if parallel == 0:
        return value
    else:
        sourceArray = array([value])
        targetArray = array([0.])
        Scientific.MPI.world.allreduce(sourceArray, targetArray, Scientific.MPI.sum)
        return targetArray[0]

def GetTotalArray(serialArray):
    if parallel == 0:
        return serialArray
    else:
        targetArray = zeros(serialArray.shape, serialArray.typecode())
        Scientific.MPI.world.allreduce(serialArray, targetArray, Scientific.MPI.sum)
        return targetArray



def GetAverageScalar(value):
    return GetTotalScalar(value) / nProcs


def GetMaxScalar(value):
    if parallel == 0:
        return value
    else:
        sourceArray = array([value])
        targetArray = array([0.])
        Scientific.MPI.world.allreduce(sourceArray, targetArray, Scientific.MPI.max)
        return targetArray[0]


def GetUnionOfList(list):
    if parallel == 0:
        return list
    else:

        if ismaster:
            unionDic = {}
            for item in list:
                unionDic[item] = 1
            for pIdx in xrange(1, nProcs):
                # receive from batch from processor pIdx
                listStr, src, tag = Scientific.MPI.world.receiveString(pIdx,0) # ????
                assert(src == pIdx)
                # unpickle it to get the list of keys of the other proc
                otherList = pickle.loads(listStr)
                for item in otherList:
                    unionDic[item] = 1
            # create a string containing the complete key set
            unionKeys = unionDic.keys()
            unionKeyStr = pickle.dumps(unionKeys)
        else:
            listStr = pickle.dumps(list)
            Scientific.MPI.world.send(listStr, 0, 0)

        # send union key string to all processors
        if ismaster:
            for pIdx in xrange(1, nProcs):
                Scientific.MPI.world.send(unionKeyStr, pIdx, 1)
        else:
            unionKeyStr, src, tag = Scientific.MPI.world.receiveString(0, 1)
            assert(src == 0)
            unionKeys = pickle.loads(unionKeyStr)
    
        return unionKeys
