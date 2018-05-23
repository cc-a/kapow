from simtk import openmm
# from simtk.openmm import app
from functools import partial
from six import with_metaclass
from collections import namedtuple, defaultdict
import inspect
import sys
from .pythonize import Pythonize
from .wrappers import DictWrapper, ValueWrapper, build_ArrayWrapper
from .class_map import class_map
from . import app

# TODO

# test reload logic

# possible for getters with a getNum and a getXXXName could build a dictwrapper

# Individual metaclass for integrator special cases


def print_args(func):
    """For debugging"""
    def wrapped(*args):
        print(args)
        return func(*args)
    return wrapped


class Context(with_metaclass(Pythonize, openmm.Context)):
    exclude = ['getParameter', 'setParameter', 'getParameters']

    @DictWrapper
    def parameters(self, key):
        return openmm.Context.getParameter(self, key)

    @parameters.setter
    def parameters(self, key, val):
        return openmm.Context.setParameter(self, key, val)

    @parameters.ranger
    def parameters(self):
        return openmm.Context.getParameters(self)


class System(with_metaclass(Pythonize, openmm.System)):
    exclude = ['getForces', 'getVirtualSite',
               'setVirtualSite', 'isVirtualSite']

    @partial(DictWrapper, ffilter=openmm.System.isVirtualSite,
             frange=lambda x: range(openmm.System.getNumParticles(x)))
    def virtualSites(self, key):
        return openmm.System.getVirtualSite(self, key)

    @virtualSites.setter
    def virtualSites(self, key, value):
        return openmm.System.setVirtualSite(self, key, value)


class Platform(with_metaclass(Pythonize, openmm.Platform)):
    exclude = ['getPropertyDefaultValue', 'setPropertyDefaultValue',
               'getNumPlatforms', 'getPlatform', 'getPlatformByName']
    preserve = ['getPropertyValue', 'setPropertyValue']

    @DictWrapper
    def propertyDefaults(self, key):
        return openmm.Platform.getPropertyDefaultValue(self, key)

    @propertyDefaults.setter
    def propertyDefaults(self, key, value):
        return openmm.Platform.setPropertyDefaultValue(self, key, value)

    propertyDefaults = propertyDefaults.ranger(
        openmm.Platform.getPropertyNames)

    @DictWrapper
    def platformsByName(key):
        return openmm.Platform.getPlatformByName(key)

    @platformsByName.ranger
    def platformsByName():
        return [openmm.Platform.getPlatform(i).getName()
                for i in range(openmm.Platform.getNumPlatforms())]


class AmoebaMultipoleForce(
        with_metaclass(Pythonize, openmm.AmoebaMultipoleForce)):
    exclude = ['getCovalentMap', 'setCovalentMap', 'getCovalentMaps']

    @DictWrapper
    def covalentMaps(self, key):
        return openmm.AmoebaMultipoleForce.getCovalentMaps(self, key)

    @covalentMaps.setter
    def covalentMaps(self, key, value):
        if len(value) != 8:
            raise TypeError('Covalent Map must be a tuple of length 8')
        for i, val in enumerate(value):
            openmm.AmoebaMultipoleForce.setCovalentMap(self, key, i, val)

    @covalentMaps.ranger
    def covalentMaps(self):
        return range(openmm.AmoebaMultipoleForce.getNumMultipoles(self))


class TwoParticleAverageSite(
        with_metaclass(Pythonize, openmm.TwoParticleAverageSite)):
    exclude = ['getNumParticles', 'getParticle', 'getWeight']
    particles = build_ArrayWrapper(
        'particles',
        openmm.TwoParticleAverageSite,
        'getNumParticles',
        ['getParticle', 'getWeight'],
        member_wrapper=namedtuple('Particle', ('index', 'weight')))


class CustomCentroidBondForce(
        with_metaclass(Pythonize, openmm.CustomCentroidBondForce)):
    exclude = ['getNumGroupsPerBond']
    numGroupsPerBond = ValueWrapper(
        openmm.CustomCentroidBondForce, 'getNumGroupsPerBond')


class CustomCompoundBondForce(
        with_metaclass(Pythonize, openmm.CustomCompoundBondForce)):
    exclude = ['getNumParticlesPerBond']
    numParticlesPerBond = ValueWrapper(
        openmm.CustomCompoundBondForce, 'getNumParticlesPerBond')


class CustomManyParticleForce(
        with_metaclass(Pythonize, openmm.CustomManyParticleForce)):
    exclude = ['getNumParticlesPerSet', 'getTypeFilter']
    numParticlesPerSet = ValueWrapper(
        openmm.CustomManyParticleForce, 'getNumParticlesPerSet')

    typeFilters = build_ArrayWrapper(
        'typeFilters',
        openmm.CustomManyParticleForce,
        'getNumParticlesPerSet',
        ['getTypeFilter'],
        ['setTypeFilter']
    )


class CustomIntegrator(
        with_metaclass(Pythonize, openmm.CustomIntegrator)):
    exclude = ['getPerDofVariable', 'getPerDofVariableName',
               'getPerDofVariableByName', 'setPerDofVariableByName',
               'addPerDofVariable', 'getNumPerDofVariables',
               'addGlobalVariable', 'getGlobalVariable',
               'getGlobalVariableName', 'setGlobalVariable',
               'getGlobalVariableByName', 'setGlobalVariableByName']

    perDofVariables = DictWrapper(
        getattr(openmm.CustomIntegrator, 'getPerDofVariableByName'),
        getattr(openmm.CustomIntegrator, 'setPerDofVariableByName'),
        fdel=None,
        doc=None,
        frange=lambda x: xrange(
            getattr(openmm.CustomIntegrator, 'getNumPerDofVariables')(x)),
        ffilter=None,
        fchange=getattr(openmm.CustomIntegrator, 'getPerDofVariableName'),
        fadd=getattr(openmm.CustomIntegrator, 'addPerDofVariable')
    )

    globalVariables = DictWrapper(
        getattr(openmm.CustomIntegrator, 'getGlobalVariableByName'),
        getattr(openmm.CustomIntegrator, 'setGlobalVariableByName'),
        fdel=None,
        doc=None,
        frange=lambda x: xrange(
            getattr(openmm.CustomIntegrator, 'getNumGlobalVariables')(x)),
        ffilter=None,
        fchange=getattr(openmm.CustomIntegrator, 'getGlobalVariableName'),
        fadd=getattr(openmm.CustomIntegrator, 'addGlobalVariable')
    )


class DualAMDIntegrator(
        with_metaclass(Pythonize, openmm.DualAMDIntegrator)):
    exclude = ['getPerDofVariable', 'getPerDofVariableName',
               'getPerDofVariableByName', 'setPerDofVariableByName',
               'addPerDofVariable', 'getNumPerDofVariables',
               'addGlobalVariable', 'getGlobalVariable',
               'getGlobalVariableName', 'setGlobalVariable',
               'getGlobalVariableByName', 'setGlobalVariableByName'
    ]

    perDofVariables = DictWrapper(
        getattr(openmm.DualAMDIntegrator, 'getPerDofVariableByName'),
        getattr(openmm.DualAMDIntegrator, 'setPerDofVariableByName'),
        fdel=None,
        doc=None,
        frange=lambda x: xrange(
            getattr(openmm.DualAMDIntegrator, 'getNumPerDofVariables')(x)),
        ffilter=None,
        fchange=getattr(openmm.DualAMDIntegrator, 'getPerDofVariableName'),
        fadd=getattr(openmm.DualAMDIntegrator, 'addPerDofVariable')
    )

    globalVariables = DictWrapper(
        getattr(openmm.DualAMDIntegrator, 'getGlobalVariableByName'),
        getattr(openmm.DualAMDIntegrator, 'setGlobalVariableByName'),
        fdel=None,
        doc=None,
        frange=lambda x: xrange(
            getattr(openmm.DualAMDIntegrator, 'getNumGlobalVariables')(x)),
        ffilter=None,
        fchange=getattr(openmm.DualAMDIntegrator, 'getGlobalVariableName'),
        fadd=getattr(openmm.DualAMDIntegrator, 'addGlobalVariable')
    )


class AMDIntegrator(
        with_metaclass(Pythonize, openmm.AMDIntegrator)):
    exclude = ['getPerDofVariable', 'getPerDofVariableName',
               'getPerDofVariableByName', 'setPerDofVariableByName',
               'addPerDofVariable', 'getNumPerDofVariables',
               'addGlobalVariable', 'getGlobalVariable',
               'getGlobalVariableName', 'setGlobalVariable',
               'getGlobalVariableByName', 'setGlobalVariableByName'
    ]

    perDofVariables = DictWrapper(
        getattr(openmm.AMDIntegrator, 'getPerDofVariableByName'),
        getattr(openmm.AMDIntegrator, 'setPerDofVariableByName'),
        fdel=None,
        doc=None,
        frange=lambda x: xrange(
            getattr(openmm.AMDIntegrator, 'getNumPerDofVariables')(x)),
        ffilter=None,
        fchange=getattr(openmm.AMDIntegrator, 'getPerDofVariableName'),
        fadd=getattr(openmm.AMDIntegrator, 'addPerDofVariable')
    )

    globalVariables = DictWrapper(
        getattr(openmm.AMDIntegrator, 'getGlobalVariableByName'),
        getattr(openmm.AMDIntegrator, 'setGlobalVariableByName'),
        fdel=None,
        doc=None,
        frange=lambda x: xrange(
            getattr(openmm.AMDIntegrator, 'getNumGlobalVariables')(x)),
        ffilter=None,
        fchange=getattr(openmm.AMDIntegrator, 'getGlobalVariableName'),
        fadd=getattr(openmm.AMDIntegrator, 'addGlobalVariable')
    )


class MTSIntegrator(
        with_metaclass(Pythonize, openmm.MTSIntegrator)):
    exclude = ['getPerDofVariable', 'getPerDofVariableName',
               'getPerDofVariableByName', 'setPerDofVariableByName',
               'addPerDofVariable', 'getNumPerDofVariables',
               'addGlobalVariable', 'getGlobalVariable',
               'getGlobalVariableName', 'setGlobalVariable',
               'getGlobalVariableByName', 'setGlobalVariableByName'
    ]

    perDofVariables = DictWrapper(
        getattr(openmm.MTSIntegrator, 'getPerDofVariableByName'),
        getattr(openmm.MTSIntegrator, 'setPerDofVariableByName'),
        fdel=None,
        doc=None,
        frange=lambda x: xrange(
            getattr(openmm.MTSIntegrator, 'getNumPerDofVariables')(x)),
        ffilter=None,
        fchange=getattr(openmm.MTSIntegrator, 'getPerDofVariableName'),
        fadd=getattr(openmm.MTSIntegrator, 'addPerDofVariable')
    )

    globalVariables = DictWrapper(
        getattr(openmm.MTSIntegrator, 'getGlobalVariableByName'),
        getattr(openmm.MTSIntegrator, 'setGlobalVariableByName'),
        fdel=None,
        doc=None,
        frange=lambda x: xrange(
            getattr(openmm.MTSIntegrator, 'getNumGlobalVariables')(x)),
        ffilter=None,
        fchange=getattr(openmm.MTSIntegrator, 'getGlobalVariableName'),
        fadd=getattr(openmm.MTSIntegrator, 'addGlobalVariable')
    )

    
class AMDForceGroupIntegrator(
        with_metaclass(Pythonize, openmm.AMDForceGroupIntegrator)):
    exclude = ['getPerDofVariable', 'getPerDofVariableName',
               'getPerDofVariableByName', 'setPerDofVariableByName',
               'addPerDofVariable', 'getNumPerDofVariables',
               'addGlobalVariable', 'getGlobalVariable',
               'getGlobalVariableName', 'setGlobalVariable',
               'getGlobalVariableByName', 'setGlobalVariableByName'
    ]

    perDofVariables = DictWrapper(
        getattr(openmm.AMDForceGroupIntegrator, 'getPerDofVariableByName'),
        getattr(openmm.AMDForceGroupIntegrator, 'setPerDofVariableByName'),
        fdel=None,
        doc=None,
        frange=lambda x: xrange(
            getattr(openmm.AMDForceGroupIntegrator,
                    'getNumPerDofVariables')(x)),
        ffilter=None,
        fchange=getattr(openmm.AMDForceGroupIntegrator,
                        'getPerDofVariableName'),
        fadd=getattr(openmm.AMDForceGroupIntegrator, 'addPerDofVariable')
    )

    globalVariables = DictWrapper(
        getattr(openmm.AMDForceGroupIntegrator, 'getGlobalVariableByName'),
        getattr(openmm.AMDForceGroupIntegrator, 'setGlobalVariableByName'),
        fdel=None,
        doc=None,
        frange=lambda x: xrange(
            getattr(openmm.AMDForceGroupIntegrator,
                    'getNumGlobalVariables')(x)),
        ffilter=None,
        fchange=getattr(openmm.AMDForceGroupIntegrator,
                        'getGlobalVariableName'),
        fadd=getattr(openmm.AMDForceGroupIntegrator, 'addGlobalVariable')
    )


class AmoebaTorsionTorsionForce(
        with_metaclass(Pythonize, openmm.AmoebaTorsionTorsionForce)):
    exclude = ['getNumTorsionTorsionGrids', 'getTorsionTorsionGrid']
    torsionTorsionGrids = build_ArrayWrapper(
        'torsionTorsionGrids',
        openmm.AmoebaTorsionTorsionForce,
        'getNumTorsionTorsionGrids',
        ['getTorsionTorsionGrid'],
        ['setTorsionTorsionGrid']
    )


class AmoebaVdwForce(
        with_metaclass(Pythonize, openmm.AmoebaVdwForce)):
    exclude = ['getParticleExclusions', 'setParticleExclusions']
    exclusions = build_ArrayWrapper(
        'exclusions',
        openmm.AmoebaVdwForce,
        'getNumParticles',
        ['getParticleExclusions'],
        ['setParticleExclusions']
    )


# specify methods to be excluded from adapter object
exclusions = defaultdict(
    lambda: [],
    {
        'System': ['getForces'],  # convenience function that is not needed
    }
)

# specify methods to be preserved in adapter objects without modification
preserves = defaultdict(
    lambda: [],
    {}
)

# specify classes to ignore
skip = []

module = sys.modules[__name__]
for name in [n for n in dir(openmm) if inspect.isclass(getattr(openmm, n))]:
    if name in skip:
        continue
    # does this logic cause problems for reloads?
    if getattr(module, name, None) is None:
        openmm_class = getattr(openmm, name)
        new_class = Pythonize(name,
                              [openmm_class],
                              {'exclude': exclusions[name],
                               'preserve': preserves[name]})
        setattr(module, name, new_class)
    class_map[getattr(openmm, name)] = getattr(module, name)
