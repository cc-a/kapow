#!/usr/bin/env python2
import unittest
import adapter
from simtk import openmm
from simtk.openmm import app

from simtk.unit import nanometer, elementary_charge, kilojoule_per_mole
from simtk.unit import dalton, radian


# TODO

# Test DictWrapper on objects that already contains items, e.g. a system object
# with virtual sites already created when wrapped.

# Unit tests for wrapper objects

class TestWrappedClasses(unittest.TestCase):
    def inst_test(self, inst, attrs):
        # for attr in inst.__class__.__dict__:
        #     self.assertNotIn('get', attr,
        #                      'Attribute %s starts with get' % attr)
        all_attrs = dir(inst)
        for attr in attrs:
            self.assertTrue(hasattr(inst, attr),
                            msg='Unable to find attribute "%s"' % attr)
            all_attrs.pop(all_attrs.index(attr))

        # for attr in all_attrs:
        #     if not attr.startswith('_'):
        #         raise Exception('Found unexpected non-magic method %s' % attr)
            

    def assignValues(self, inst, attr, val1, val2=None):
        setattr(inst, attr, val1)
        self.assertEqual(getattr(inst, attr), val1)
        if val2 is not None:
            setattr(inst, attr, val2)
            self.assertEqual(getattr(inst, attr), val2)

    def testNonbondedForce(self):
        base_inst = openmm.NonbondedForce()
        inst = adapter.NonbondedForce(base_inst)

        attrs_to_check = ['CutoffNonPeriodic',
                          'CutoffPeriodic',
                          'Ewald',
                          'NoCutoff',
                          'PME',
                          'exceptions',
                          'particles',
                          'createExceptionsFromBonds',
                          'cutoffDistance',
                          'forceGroup',
                          'ewaldErrorTolerance',
                          'nonbondedMethod',
                          'PMEParameters',
                          'getPMEParametersInContext',
                          'reactionFieldDielectric',
                          'reciprocalSpaceForceGroup',
                          'switchingDistance',
                          'useDispersionCorrection',
                          'useSwitchingFunction',
                          'updateParametersInContext',
                          'usesPeriodicBoundaryConditions']
        self.inst_test(inst, attrs_to_check)

        self.assertTrue(callable(inst.createExceptionsFromBonds))
        self.assignValues(inst, 'cutoffDistance', 2. * nanometer)
        self.assignValues(inst, 'ewaldErrorTolerance', 0.005)
        self.assignValues(inst, 'forceGroup', 4)
        self.assignValues(inst, 'nonbondedMethod', inst.Ewald, inst.NoCutoff)
        self.assignValues(inst, 'PMEParameters', [1. / nanometer, 1, 1, 1])
        self.assignValues(inst, 'reactionFieldDielectric', 1.)
        self.assignValues(inst, 'reciprocalSpaceForceGroup', 5)
        self.assignValues(inst, 'cutoffDistance', 5. * nanometer)
        self.assignValues(inst, 'useDispersionCorrection', True, False)
        self.assignValues(inst, 'useSwitchingFunction', True, False)

        inst.particles.append((0.0, 0.0, 0.0))
        self.assertEqual(inst.particles[0],
                         [0. * elementary_charge,
                          0. * nanometer,
                          0. * kilojoule_per_mole])
        self.assertEqual(inst.particles[0],
                         inst.wrapped_object.getParticleParameters(0))
        inst.exceptions.append((1, 2, 0.0, 0.0, 0.0))
        self.assertEqual(inst.exceptions[0],
                         [1, 2,
                          0. * elementary_charge**2,
                          0. * nanometer,
                          0. * kilojoule_per_mole])
        self.assertEqual(inst.exceptions[0],
                         inst.wrapped_object.getExceptionParameters(0))

    def testHarmonicBondForce(self):
        """A very comprehensive test for a very simpler class,
        exhaustively checks that all methods are defined as expected
        and can be accessed or added to as required
        """

        base_inst = openmm.HarmonicBondForce()
        inst = adapter.HarmonicBondForce(base_inst)
        attrs_to_check = ['bonds',
                          'forceGroup',
                          'updateParametersInContext',
                          'usesPeriodicBoundaryConditions']
        self.inst_test(inst, attrs_to_check)

        inst.bonds.append((0, 1, 0. * nanometer,
                          0. * kilojoule_per_mole / nanometer**2))
        self.assertEqual(inst.bonds[0],
                         [0, 1, 0. * nanometer,
                          0. * kilojoule_per_mole / nanometer**2])
        self.assertEqual(inst.bonds[0],
                         inst.wrapped_object.getBondParameters(0))

        inst.forceGroup = 2
        self.assertEqual(inst.forceGroup, 2)
        self.assignValues(inst, 'usesSwitchingFunction', True, False)
        self.assertRaisesRegexp(TypeError,
                                'updateParametersInContext\(\) takes exactly '
                                '2 arguments \(1 given\)',
                                callable_obj=inst.updateParametersInContext)

    def testHarmonicAngleForce(self):
        """A very comprehensive test for a very simpler class,
        exhaustively checks that all methods are defined as expected
        and can be accessed or added to as required
        """

        base_inst = openmm.HarmonicAngleForce()
        inst = adapter.HarmonicAngleForce(base_inst)
        attrs_to_check = ['angles',
                          'forceGroup',
                          'updateParametersInContext',
                          'usesPeriodicBoundaryConditions']
        self.inst_test(inst, attrs_to_check)

        inst.angles.append((0, 1, 2, 0. * radian,
                            0. * kilojoule_per_mole / radian**2))
        self.assertEqual(inst.angles[0],
                         [0, 1, 2, 0. * radian,
                          0. * kilojoule_per_mole / radian**2])
        self.assertEqual(inst.angles[0],
                         inst.wrapped_object.getAngleParameters(0))

        inst.forceGroup = 2
        self.assertEqual(inst.forceGroup, 2)
        self.assignValues(inst, 'usesSwitchingFunction', True, False)
        self.assertRaisesRegexp(TypeError,
                                'updateParametersInContext\(\) takes exactly '
                                '2 arguments \(1 given\)',
                                callable_obj=inst.updateParametersInContext)

    def testWrappingAlreadyCreated(self):
        base_inst = openmm.System()
        base_inst.addParticle(1.)
        base_inst.addParticle(2.)
        base_inst.addParticle(3.)
        vs = openmm.TwoParticleAverageSite(0, 1, 0.5, 0.5)
        base_inst.setVirtualSite(2, vs)
        inst = adapter.System(base_inst)

        self.assertEqual(len(inst.particles), 3)
        self.assertEqual(
            [p.value_in_unit(dalton) for p in inst.particles], [1., 2, 3.])
        self.assertTrue(inst.wrapped_object.isVirtualSite(2))
        self.assertEqual(inst.virtualSites.keys(), [2])

    def testSystem(self):
        base_inst = openmm.System()
        inst = adapter.System(base_inst)

        attrs_to_check = ['constraints',
                          'particles',
                          'forces',
                          'defaultPeriodicBoxVectors',
                          'virtualSites',
                          'usesPeriodicBoundaryConditions']
        self.inst_test(inst, attrs_to_check)

        self.assertEqual(len(inst.particles), 0)
        inst.particles.append((1.,))
        inst.particles.append((2.,))
        self.assertEqual(len(inst.particles), 2)
        self.assertEqual(inst.particles[0], 1. * dalton)
        self.assertEqual(inst.particles[1], 2. * dalton)
        inst.particles[0] = 3.
        self.assertEqual(inst.particles[0], 3. * dalton)

        self.assertEqual(len(inst.constraints), 0)
        inst.constraints.append((0, 1, 1. * nanometer))
        self.assertEqual(inst.constraints[0],
                         [0, 1, 1. * nanometer])
        self.assertEqual(len(inst.constraints), 1)
        inst.constraints[0] = (0, 1, 2. * nanometer)
        self.assertEqual(inst.constraints[0],
                         [0, 1, 2. * nanometer])

        base_hbf = openmm.HarmonicBondForce()
        hbf = adapter.HarmonicBondForce(base_hbf)
        self.assertEqual(len(inst.forces), 0)
        inst.forces.append(hbf)
        self.assertTrue(isinstance(inst.forces[0], adapter.HarmonicBondForce))
        self.assertEqual(len(inst.forces), 1)
        inst.forces.pop(0)
        self.assertEqual(len(inst.forces), 0)

        # Add a third particle to be a virtual site for the first two
        inst.particles.append((3.,))
        vs = openmm.TwoParticleAverageSite(0, 1, 0.5, 0.5)
        wrapped_vs = adapter.TwoParticleAverageSite(vs)
        self.assertEqual(len(inst.virtualSites), 0)
        inst.virtualSites[2] = wrapped_vs
        self.assertEqual(len(inst.virtualSites), 1)
        self.assertEqual(inst.virtualSites.keys(), [2])
        vs_out = inst.virtualSites[2]
        self.assertEqual(wrapped_vs.particles, vs_out.particles)

    def testCustomNonbondedForce(self):
        base_inst = openmm.CustomNonbondedForce('r')
        inst = adapter.CustomNonbondedForce(base_inst)
        attrs_to_check = [
            'exclusions',
            'functions',
            'globalParameters',
            'interactionGroups',
            'particles',
            'perParticleParameters',
            'tabulatedFunctions',
            'createExclusionsFromBonds',
            'cutoffDistance',
            'energyFunction',
            'forceGroup',
            'switchingDistance',
            'useLongRangeCorrection',
            'useSwitchingFunction',
            'updateParametersInContext',
            'usesPeriodicBoundaryConditions'
        ]
        self.inst_test(inst, attrs_to_check)

        self.assertEqual(len(inst.globalParameters), 0)
        inst.globalParameters.append(('jon', 0))
        self.assertEqual(len(inst.globalParameters), 1)
        self.assertEqual(inst.globalParameters[0], ('jon', 0))
        inst.globalParameters[0] = ('chris', 1)

    def testPlatform(self):
        # where is getPlatformByName?
        base_inst = openmm.Platform.getPlatformByName('CPU')
        inst = adapter.Platform(base_inst)
        attrs_to_check = [
            'findPlatform',
            'defaultPluginsDirectory',
            'name',
            'platforms',
            # 'getPlatformByName', # don't know why this is absent
            'getPropertyValue',
            'setPropertyValue',
            'propertyDefaults',
            'loadPluginLibrary',
            'loadPluginsFromDirectory',
            'registerPlatform',
            'supportsDoublePrecision',
            'supportsKernels'
        ]
        self.inst_test(inst, attrs_to_check)
        self.assertEqual(inst.propertyDefaults.keys(), ['CpuThreads'])
        self.assertEqual(inst.propertyDefaults['CpuThreads'],
                         base_inst.getPropertyDefaultValue('CpuThreads'))
        inst.propertyDefaults['CpuThreads'] = '1'
        self.assertEqual(inst.propertyDefaults['CpuThreads'], '1')

    def testAmoebaMultipoleForce(self):
        base_inst = openmm.AmoebaMultipoleForce()
        inst = adapter.AmoebaMultipoleForce(base_inst)
        attrs_to_check = [
            'multipoles',
            'AEwald',
            'covalentMaps',
            'cutoffDistance',
            'getElectrostaticPotential',
            'ewaldErrorTolerance',
            'extrapolationCoefficients',
            'forceGroup',
            'getInducedDipoles',
            'mutualInducedMaxIterations',
            'mutualInducedTargetEpsilon',
            'nonbondedMethod',
            'getPMEParametersInContext',
            'pmeBSplineOrder',
            'pmeGridDimensions',
            'polarizationType',
            'updateParametersInContext',
            'usesPeriodicBoundaryConditions',
        ]
        self.inst_test(inst, attrs_to_check)

        self.assertEqual(len(inst.multipoles), 0)
        inst.multipoles.append((0., (0.,)*3, (0.,)*9, 0, 0, 0, 0, 0, 0, 0,))
        self.assertEqual(len(inst.multipoles), 1)
        self.assertEqual(len(inst.covalentMaps), 1)
        self.assertEqual(inst.covalentMaps.keys(), [0])
        self.assertEqual(
            inst.covalentMaps[0],
            base_inst.getCovalentMaps(0))
        inst.covalentMaps[0] = ((1,),) * 8
        self.assertEqual(
            inst.covalentMaps[0], ((1,),) * 8)
        with self.assertRaises(TypeError):
            inst.covalentMaps[0] = ((1,),) * 7

    def testTwoParticleAverageSite(self):
        base_inst = openmm.TwoParticleAverageSite(2, 3, 0.2, 0.8)
        inst = adapter.TwoParticleAverageSite(base_inst)
        attrs_to_check = ['particles']
        self.inst_test(inst, attrs_to_check)
        self.assertEqual(len(inst.particles), 2)
        self.assertEqual(inst.particles[0], (2, 0.2))
        self.assertEqual(inst.particles[1], (3, 0.8))
        with self.assertRaises(AttributeError):
            inst.particles[0] = (0, 0.5)

if __name__ == '__main__':
    unittest.main()
