# This file is part of GridCal.
#
# GridCal is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# GridCal is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with GridCal.  If not, see <http://www.gnu.org/licenses/>.


from GridCal.Engine.Devices.meta_devices import EditableDevice, DeviceType, GCProp


class Load(EditableDevice):
    """
    The load object implements the so-called ZIP model, in which the load can be
    represented by a combination of power (P), current(I), and impedance (Z).

    The sign convention is: Positive to act as a load, negative to act as a generator.

    Arguments:

        **name** (str, "Load"): Name of the load

        **G** (float, 0.0): Conductance in equivalent MW

        **B** (float, 0.0): Susceptance in equivalent MVAr

        **Ir** (float, 0.0): Real current in equivalent MW

        **Ii** (float, 0.0): Imaginary current in equivalent MVAr

        **P** (float, 0.0): Active power in MW

        **Q** (float, 0.0): Reactive power in MVAr

        **a** (float, 1.0): phase A share of the declared power

        **b** (float, 1.0): phase B share of the declared power

        **c** (float, 1.0): phase C share of the declared power

        **G_prof** (array, None): Numpy array with the conductance profile in equivalent MW

        **B_prof** (array, None): Numpy array with the susceptance profile in equivalent MVAr

        **Ir_prof** (array, None): Numpy array with the real current profile in equivalent MW

        **Ii_prof** (array, None): Numpy array with the imaginary current profile in equivalent MVAr

        **P_prof** (array, None): Numpy array with the active power profile in equivalent MW

        **Q_prof** (array, None): Numpy array with the reactive power profile in equivalent MVAr

        **a_prof** (array, None): Numpy array with the phase A share of the declared power

        **b_prof** (array, None): Numpy array with the phase B share of the declared power

        **c_prof** (array, None): Numpy array with the phase C share of the declared power

        **active** (bool, True): Is the load active?

        **mttf** (float, 0.0): Mean time to failure in hours

        **mttr** (float, 0.0): Mean time to recovery in hours

    """

    def __init__(self, name='Load', G=0.0, B=0.0, Ir=0.0, Ii=0.0, P=0.0, Q=0.0, a=1.0, b=1.0, c=1.0, cost=0.0,
                 G_prof=None, B_prof=None, Ir_prof=None, Ii_prof=None, P_prof=None, Q_prof=None,
                 a_prof=None, b_prof=None, c_prof=None,  active=True, mttf=0.0, mttr=0.0):

        EditableDevice.__init__(self,
                                name=name,
                                active=active,
                                device_type=DeviceType.LoadDevice,
                                editable_headers={'name': GCProp('', str, 'Load name'),
                                                   'bus': GCProp('', None, 'Connection bus name'),
                                                   'active': GCProp('', bool, 'Is the load active?'),
                                                   'P': GCProp('MW', float, 'Active power'),
                                                   'Q': GCProp('MVAr', float, 'Reactive power'),
                                                   'Ir': GCProp('MW', float,
                                                                'Active power of the current component at V=1.0 p.u.'),
                                                   'Ii': GCProp('MVAr', float,
                                                                'Reactive power of the current component at V=1.0 p.u.'),
                                                   'G': GCProp('MW', float,
                                                               'Active power of the impedance component at V=1.0 p.u.'),
                                                   'B': GCProp('MVAr', float,
                                                               'Reactive power of the impedance component at V=1.0 p.u.'),
                                                   'mttf': GCProp('h', float, 'Mean time to failure'),
                                                   'mttr': GCProp('h', float, 'Mean time to recovery'),
                                                   'a': GCProp('p.u.', float,
                                                               'phase A share of the declared power'),
                                                   'b': GCProp('p.u.', float,
                                                               'phase B share of the declared power'),
                                                   'c': GCProp('p.u.', float,
                                                               'phase C share of the declared power'),
                                                   'Cost': GCProp('e/MWh', float,
                                                                  'Cost of not served energy. Used in OPF.')
                                                  },
                                non_editable_attributes=list(),
                                properties_with_profile={'P': 'P_prof',
                                                         'Q': 'Q_prof',
                                                         'Ir': 'Ir_prof',
                                                         'Ii': 'Ii_prof',
                                                         'G': 'G_prof',
                                                         'B': 'B_prof',
                                                         'a': 'a_prof',
                                                         'b': 'b_prof',
                                                         'c': 'c_prof',
                                                         'Cost': 'Cost_prof'})

        # connection bus
        self.bus = None

        # mean time to failure
        self.mttf = mttf

        # mean time to repair
        self.mttr = mttr

        self.Cost = cost

        self.Cost_prof = None

        # Impedance in equivalent MVA
        # active and reactive admittance in MW and MVAr at v=1.0 p.u
        self.G = G
        self.B = B

        # active and reactive current in MW and MVAr at v=1.0 p.u
        self.Ir = Ir
        self.Ii = Ii

        # active and reactive power in MW and MVAr
        self.P = P
        self.Q = Q

        # active and reactive admittance profiles in MW and MVAr at v=1.0 p.u
        self.G_prof = G_prof
        self.B_prof = B_prof

        # active and reactive current profiles in MW and MVAr at v=1.0 p.u
        self.Ir_prof = Ir_prof
        self.Ii_prof = Ii_prof

        # active and reactive power profiles in MW and MVAr
        self.P_prof = P_prof
        self.Q_prof = Q_prof

        # shape per phase
        self.a = a
        self.b = b
        self.c = c

        # share per phase profiles
        self.a_prof = a_prof
        self.b_prof = b_prof
        self.c_prof = c_prof

    def copy(self):

        load = Load(a=self.a,
                    b=self.b,
                    c=self.c,
                    a_prof=self.a_prof,
                    b_prof=self.b_prof,
                    c_prof=self.c_prof)

        load.name = self.name

        # Impedance (MVA)
        load.G = self.G
        load.B = self.B

        # Current (MVA)
        load.Ir = self.Ir
        load.Ii = self.Ii

        # Power (MVA)
        load.P = self.P
        load.Q = self.Q

        # Impedance (MVA)
        load.G_prof = self.G_prof
        load.B_prof = self.B_prof

        # Current (MVA)
        load.Ir_prof = self.Ir_prof
        load.Ii_prof = self.Ii_prof

        # Power (MVA)
        load.P_prof = self.P_prof
        load.Q_prof = self.Q_prof

        load.mttf = self.mttf

        load.mttr = self.mttr

        return load

    def get_json_dict(self, id, bus_dict):
        """
        Get json dictionary
        :param id: ID: Id for this object
        :param bus_dict: Dictionary of buses [object] -> ID
        :return:
        """
        return {'id': id,
                'type': 'load',
                'phases': 'ps',
                'name': self.name,
                'bus': bus_dict[self.bus],
                'active': self.active,
                'G': self.G,
                'B': self.B,
                'Ir': self.Ir,
                'Ii': self.Ii,
                'P': self.P,
                'Q': self.Q}

