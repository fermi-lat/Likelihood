"""
SourceModel interface to allow for manipulation of fit parameters.

@author J. Chiang <jchiang@slac.stanford.edu>
"""
#
# $Header$
#

import pyLike

class SourceModel(object):
    def __init__(self, logLike):
        self.logLike = logLike
        srcNames = pyLike.StringVector()
        logLike.getSrcNames(srcNames)
        self.srcNames = tuple(srcNames)
    def __getitem__(self, srcName):
        return Source(self.logLike.getSource(srcName))

class Source(object):
    def __init__(self, src):
        self.src = src
        funcs = src.getSrcFuncs()
        self.funcs = {}
        for item in funcs.keys():
            self.funcs[item] = Function(funcs[item])
        self.__dict__.update(self.funcs)
    def __getitem__(self, name):
        return self.funcs[name]

class Function(object):
    def __init__(self, func):
        self.func = func
        names = pyLike.StringVector()
        func.getParamNames(names)
        self.paramNames = list(names)
    def __getitem__(self, name):
        return self.func.getParamValue(name)
    def __setitem__(self, name, value):
        self.func.setParam(name, value)
    def __call__(self, name):
        return self.func.getParam(name)
