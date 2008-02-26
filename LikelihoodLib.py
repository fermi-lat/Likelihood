#$Id: LikelihoodLib.py,v 1.2 2008/02/20 16:55:44 golpa Exp $
def generate(env, **kw):
    if not kw.get('depsOnly',0):
        env.Tool('addLibrary', library=['Likelihood'])
    env.Tool('astroLib')
    env.Tool('xmlBaseLib')
    env.Tool('tipLib')
    env.Tool('evtbinLib')
    env.Tool('map_toolsLib')
    env.Tool('optimizersLib')
    env.Tool('irfLoaderLib')
    env.Tool('st_facilitiesLib')
    env.Tool('dataSubselectorLib')
    env.Tool('hoopsLib')
    env.Tool('st_appLib')
    env.Tool('st_graphLib')
    env.Tool('healpixLib')
    env.Tool('addLibrary', library=env['cfitsioLibs'])
    env.Tool('addLibrary', library=env['fftwLibs'])
    
def exists(env):
    return 1
