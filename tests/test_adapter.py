#!/usr/bin/env python
import unittest
import kapow
from simtk import openmm
from simtk.unit import nanometer, elementary_charge, kilojoule_per_mole
from simtk.unit import dalton, radian


# TODO

# Unit tests for wrapper objects

# AmoebaTorsionTorsionForce - the grids are weirdc

# AmoebaVdwForce - the exclusions are weird

# finish documentation!!!

class TestWrappedClasses(unittest.TestCase):
    def inst_test(self, inst, attrs):
        all_attrs = dir(inst)
        for attr in attrs:
            self.assertTrue(hasattr(inst, attr),
                            msg='Unable to find attribute "%s"' % attr)
            all_attrs.pop(all_attrs.index(attr))

    def assignValues(self, inst, attr, val1, val2=None):
        setattr(inst, attr, val1)
        self.assertEqual(getattr(inst, attr), val1)
        if val2 is not None:
            setattr(inst, attr, val2)
            self.assertEqual(getattr(inst, attr), val2)

    def assignArray(self, array, val, val2=None):
        self.assertEqual(len(array), 0)
        array.append(val)
        self.assertEqual(len(array), 1)
        self.assertEqual(
            array[0],
            array.member_wrapper(*val)
        )

        if val2 is not None:
            array[0] = val2
            self.assertEqual(
                array[0],
                array.member_wrapper(*val2)
            )

    def assignDict(self, array, key, val, val2=None):
        self.assertEqual(len(array), 0)
        array[key] = (val)
        self.assertEqual(len(array), 1)
        self.assertEqual(
            array[key],
            val
        )

        if val2 is not None:
            array[key] = val2
            self.assertEqual(
                array[key],
                val2
            )

    def testNonbondedForce(self):
        inst = kapow.NonbondedForce()

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

        self.assignArray(inst.particles,
                       (0. * elementary_charge,
                        0. * nanometer,
                        0. * kilojoule_per_mole),
                       (1. * elementary_charge,
                        1. * nanometer,
                        1. * kilojoule_per_mole))
        self.assignArray(inst.exceptions,
                       (0, 0, 0. * elementary_charge**2,
                        0. * nanometer, 0. * kilojoule_per_mole),
                       (1, 1, 1. * elementary_charge**2,
                        1. * nanometer, 1. * kilojoule_per_mole))

    def testHarmonicBondForce(self):
        """A very comprehensive test for a very simpler class,
        exhaustively checks that all methods are defined as expected
        and can be accessed or added to as required
        """

        inst = kapow.HarmonicBondForce()
        attrs_to_check = ['bonds',
                          'forceGroup',
                          'updateParametersInContext',
                          'usesPeriodicBoundaryConditions']
        self.inst_test(inst, attrs_to_check)
        self.assignArray(
            inst.bonds,
            (0, 1, 0. * nanometer, 0. * kilojoule_per_mole / nanometer**2))

        inst.forceGroup = 2
        self.assertEqual(inst.forceGroup, 2)
        self.assignValues(inst, 'usesSwitchingFunction', True, False)
        # self.assertRaisesRegexp(TypeError,
        #                         'updateParametersInContext\(\) takes exactly'
        #                         ' 2 arguments \(1 given\)',
        #                         callable_obj=inst.updateParametersInContext)

    def testHarmonicAngleForce(self):
        """A very comprehensive test for a very simpler class,
        exhaustively checks that all methods are defined as expected
        and can be accessed or added to as required
        """

        inst = kapow.HarmonicAngleForce()
        attrs_to_check = ['angles',
                          'forceGroup',
                          'updateParametersInContext',
                          'usesPeriodicBoundaryConditions']
        self.inst_test(inst, attrs_to_check)
        self.assignArray(
            inst.angles,
            (0, 0, 0, 0. * radian, 0. * kilojoule_per_mole / radian**2),
            (1, 1, 1, 1. * radian, 1. * kilojoule_per_mole / radian**2))

        inst.forceGroup = 2
        self.assertEqual(inst.forceGroup, 2)
        self.assignValues(inst, 'usesSwitchingFunction', True, False)

    def testWrappingAlreadyCreated(self):
        base_inst = openmm.System()
        base_inst.addParticle(1.)
        base_inst.addParticle(2.)
        base_inst.addParticle(3.)
        vs = openmm.TwoParticleAverageSite(0, 1, 0.5, 0.5)
        base_inst.setVirtualSite(2, vs)
        inst = kapow.System.Wrap(base_inst)

        self.assertEqual(len(inst.particles), 3)
        self.assertEqual(
            [p.mass.value_in_unit(dalton) for p in inst.particles],
            [1., 2, 3.])
        self.assertTrue(inst.wrapped_object.isVirtualSite(2))
        self.assertEqual(inst.virtualSites.keys(), [2])

    def testSystem(self):
        inst = kapow.System()
        attrs_to_check = ['constraints',
                          'particles',
                          'forces',
                          'defaultPeriodicBoxVectors',
                          'virtualSites',
                          'usesPeriodicBoundaryConditions']
        self.inst_test(inst, attrs_to_check)
        self.assignArray(inst.particles, (1 * dalton,), (2 * dalton,))
        self.assertIsInstance(inst.particles[0], inst.particles.member_wrapper)

        self.assignArray(inst.constraints,
                         (0, 0, 0. * nanometer),
                         (1, 1, 1. * nanometer))

        hbf = kapow.HarmonicBondForce()
        self.assertEqual(len(inst.forces), 0)
        inst.forces.append(hbf)
        # this is not wrapped by a namedtuple
        self.assertIsInstance(inst.forces[0],
                              kapow.HarmonicBondForce)
        self.assertEqual(len(inst.forces), 1)
        inst.forces.pop(0)
        self.assertEqual(len(inst.forces), 0)

        # Add two particles to be a virtual site for the first two
        vs = kapow.TwoParticleAverageSite(0, 1, 0.5, 0.5)
        self.assertEqual(len(inst.virtualSites), 0)
        inst.virtualSites[0] = inst.virtualSites.member_wrapper(vs)
        self.assertEqual(len(inst.virtualSites), 1)
        self.assertEqual(inst.virtualSites.keys(), [0])
        vs_out = inst.virtualSites[0]
        self.assertEqual(vs.particles, vs_out.particles)

    def testCustomNonbondedForce(self):
        inst = kapow.CustomNonbondedForce('r')
        inst.energyFunction
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

        func_args = (0, 1), 0., 1.
        tfunc = kapow.Continuous1DFunction(*func_args)
        inst.tabulatedFunctions.append(('first', tfunc))

        # tabulatedFunctions returns any added function cast as a
        # TabulatedFunction object
        name, tfunc = inst.tabulatedFunctions[0]
        self.assertEqual(name,
                         inst.wrapped_object.getTabulatedFunctionName(0))
        self.assertIsInstance(tfunc, kapow.Continuous1DFunction)

        # functions returns the parameters associated with
        # a Continuous1DFunction
        func = list(inst.functions[0][1:])
        self.assertEqual(func,
                         inst.wrapped_object.getFunctionParameters(0)[1:])

    def testPlatform(self):
        # self.assertEqual(kapow.Platform.platformsByName.keys(),
        #                  ['Reference', 'CPU'])
        keys = kapow.Platform.platformsByName.keys()
        self.assertIn('Reference', keys)
        self.assertIn('CPU', keys)
        inst = kapow.Platform.platformsByName['CPU']
        attrs_to_check = [
            'findPlatform',
            'defaultPluginsDirectory',
            'name',
            'platformsByName',
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
        self.assertEqual(inst.propertyDefaults.keys(), 
                         ['Threads', 'DeterministicForces'])
        self.assertEqual(
            inst.propertyDefaults['Threads'],
            inst.wrapped_object.getPropertyDefaultValue('Threads'))
        inst.propertyDefaults['Threads'] = '1'
        self.assertEqual(inst.propertyDefaults['Threads'], '1')

    def testAmoebaMultipoleForce(self):
        inst = kapow.AmoebaMultipoleForce()
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
            inst.wrapped_object.getCovalentMaps(0))
        inst.covalentMaps[0] = ((1,),) * 8
        self.assertEqual(
            inst.covalentMaps[0], ((1,),) * 8)
        with self.assertRaises(TypeError):
            inst.covalentMaps[0] = ((1,),) * 7

    def testTwoParticleAverageSite(self):
        inst = kapow.TwoParticleAverageSite(2, 3, 0.2, 0.8)
        attrs_to_check = ['particles']
        self.inst_test(inst, attrs_to_check)
        self.assertEqual(len(inst.particles), 2)
        self.assertEqual(inst.particles[0], (2, 0.2))
        self.assertEqual(inst.particles[1], (3, 0.8))
        with self.assertRaises(AttributeError):
            inst.particles[0] = (0, 0.5)

    def testContext(self):
        """This test checks for a basic level of suitability
        in combining app classes with openmm classes"""
        prmtop = kapow.app.AmberPrmtopFile('../tests/prot_lig1.prmtop')
        system = prmtop.createSystem()

        self.assertIsInstance(system, kapow.System)
        integrator = kapow.VerletIntegrator(0.001)
        context = kapow.Context(system, integrator)

    def testCustomManyParticleForce(self):
        cmpf = kapow.CustomManyParticleForce(3, 'r')
        self.assertEqual(3, cmpf.numParticlesPerSet)
        self.assertEqual(3, len(cmpf.typeFilters))
        self.assertEqual(cmpf.typeFilters[0], ())
        cmpf.typeFilters[0] = (0, 1)
        self.assertEqual(cmpf.typeFilters[0], (0, 1))

    def testCustomIntegrator(self):
        inst = kapow.CustomIntegrator(0.002)
        self.assertEqual(len(inst.perDofVariables), 0)
        inst.perDofVariables["oldx"] = 1
        self.assertEqual(len(inst.perDofVariables), 1)
        self.assertEqual(inst.perDofVariables["oldx"], (1., 1., 1.))
        inst.perDofVariables["oldx"] = ((2., 2., 2.),)

        self.assignDict(inst.globalVariables, "x", 1, 2)

    def testDualAMDIntegrator(self):
        # this class starts with some perDofVariables and globalVariables
        # already assigned
        inst = kapow.DualAMDIntegrator(0.002, 0, 0., 0., 0., 0.)
        self.assertEqual(len(inst.perDofVariables), 2)
        inst.perDofVariables["test"] = 1
        self.assertEqual(len(inst.perDofVariables), 3)
        self.assertEqual(inst.perDofVariables["test"], (1., 1., 1.))
        inst.perDofVariables["test"] = ((2., 2., 2.),)

        self.assertEqual(len(inst.globalVariables), 5)
        inst.globalVariables["test"] = 1
        self.assertEqual(len(inst.globalVariables), 6)
        self.assertEqual(inst.globalVariables["test"], 1.)
        inst.globalVariables["test"] = 2.

    def testAMDIntegrator(self):
        # this class starts with some perDofVariables and globalVariables
        # already assigned
        inst = kapow.AMDIntegrator(0.002, 0, 0.)
        self.assertEqual(len(inst.perDofVariables), 1)
        inst.perDofVariables["test"] = 1
        self.assertEqual(len(inst.perDofVariables), 2)
        self.assertEqual(inst.perDofVariables["test"], (1., 1., 1.))
        inst.perDofVariables["test"] = ((2., 2., 2.),)

        self.assertEqual(len(inst.globalVariables), 2)
        inst.globalVariables["test"] = 1
        self.assertEqual(len(inst.globalVariables), 3)
        self.assertEqual(inst.globalVariables["test"], 1.)
        inst.globalVariables["test"] = 2.

    def testMTSIntegrator(self):
        # this class starts with some perDofVariables and globalVariables
        # already assigned
        inst = kapow.MTSIntegrator(1, [(0, 1), (1, 2)])
        self.assertEqual(len(inst.perDofVariables), 1)
        inst.perDofVariables["test"] = 1
        self.assertEqual(len(inst.perDofVariables), 2)
        self.assertEqual(inst.perDofVariables["test"], (1., 1., 1.))
        inst.perDofVariables["test"] = ((2., 2., 2.),)

        self.assignDict(inst.globalVariables, "test", 1., 2.)

    def testAMDForceGroupIntegrator(self):
        # this class starts with some perDofVariables and globalVariables
        # already assigned
        inst = kapow.AMDForceGroupIntegrator(0.002, 0, 0., 0.)
        self.assertEqual(len(inst.perDofVariables), 2)
        inst.perDofVariables["test"] = 1
        self.assertEqual(len(inst.perDofVariables), 3)
        self.assertEqual(inst.perDofVariables["test"], (1., 1., 1.))
        inst.perDofVariables["test"] = ((2., 2., 2.),)

        self.assertEqual(len(inst.globalVariables), 3)
        inst.globalVariables["test"] = 1
        self.assertEqual(len(inst.globalVariables), 4)
        self.assertEqual(inst.globalVariables["test"], 1.)
        inst.globalVariables["test"] = 2.


if __name__ == '__main__':
    unittest.main()
