"""
SourceModel interface to allow for manipulation of fit parameters.

@author J. Chiang <jchiang@slac.stanford.edu>
"""
#
# $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/python/SrcModel.py,v 1.2 2005/02/02 00:01:16 jchiang Exp $
#
import sys
import pyLike

def ids(istart=0):
    i = istart - 1
    while True:
        i += 1
        yield(i)

class SourceModel(object):
    def __init__(self, logLike):
        self.logLike = logLike
        srcNames = pyLike.StringVector()
        logLike.getSrcNames(srcNames)
        self.srcNames = tuple(srcNames)
        self.srcs = {}
        for name in srcNames:
            self.srcs[name] = Source(self.logLike.getSource(name))
        self._walk()
    def _walk(self):
        indx = ids()
        self.params = []
        for srcName in self.srcNames:
            src = self[srcName]
            for funcName in src.funcs:
                func = src.funcs[funcName]
                for param in func.paramNames:
                    self.params.append(func.getParam(param))
                    src.funcs[funcName].appendParId(indx.next())
    def __setitem__(self, indx, value):
        self.params[indx].setValue(value)
    def __getitem__(self, srcName):
        try:
            return self.params[srcName]
        except:
            try:
                return self.srcs[srcName]
            except:
                pass
    def __repr__(self):
        lines = []
        for src in self.srcNames:
            lines.append(src)
            lines.append(self[src].__repr__('   '))
        return "\n".join(lines)

class Source(object):
    def __init__(self, src):
        self.src = src
        funcs = src.getSrcFuncs()
        self.funcs = {}
        for item in funcs.keys():
            self.funcs[item] = Function(funcs[item])
    def __getitem__(self, name):
        return self.funcs[name]
    def __repr__(self, prefix=''):
        lines = [] 
        for item in self.funcs:
            lines.append(prefix + item + ": " + self[item].genericName())
            lines.append(self[item].__repr__(prefix))
        return "\n".join(lines)

class Function(object):
    def __init__(self, func):
        self.func = func
        names = pyLike.StringVector()
        func.getParamNames(names)
        self.paramNames = list(names)
        self.params = {}
        for name in self.paramNames:
            self.params[name] = Parameter(self.func.getParam(name))
        self._parIds = []
    def __getitem__(self, name):
        return self.func.getParamValue(name)
    def __setitem__(self, name, value):
        self.func.setParam(name, value)
    def getParam(self, name):
        return self.params[name]
    def appendParId(self, indx):
        self._parIds.append(indx)
    def __repr__(self, prefix=''):
        lines = []
        for indx, parName in zip(self._parIds, self.paramNames):
            par = self.getParam(parName)
            lines.append("%-3i%s%s" % (indx, prefix, par.__repr__()))
        return "\n".join(lines)
    def __getattr__(self, attrname):
        return getattr(self.func, attrname)

class Parameter(object):
    def __init__(self, parameter):
        self.parameter = parameter
    def __repr__(self):
        par = self.parameter
        desc = ("%10s: %10.3e " % (par.getName(), par.getValue()) +
                "%10.3e %10.3e " % par.getBounds() +
                "(%10.3e)" % par.getScale())
        if not par.isFree():
            desc += " fixed"
        return desc
    def __getattr__(self, attrname):
        return getattr(self.parameter, attrname)
