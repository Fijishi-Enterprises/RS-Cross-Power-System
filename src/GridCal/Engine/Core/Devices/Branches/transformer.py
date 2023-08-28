# GridCal
# Copyright (C) 2015 - 2023 Santiago Peñate Vera
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 3 of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with this program; if not, write to the Free Software Foundation,
# Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

import pandas as pd
import numpy as np
from matplotlib import pyplot as plt

from GridCal.Engine.basic_structures import Logger
from GridCal.Engine.Core.Devices.Substation.bus import Bus
from GridCal.Engine.Core.Devices.enumerations import BranchType, TransformerControlType, WindingsConnection, BuildStatus
from GridCal.Engine.Core.Devices.Branches.templates.parent_branch import ParentBranch
from GridCal.Engine.Core.Devices.Branches.templates.transformer_type import TransformerType
from GridCal.Engine.Core.Devices.Branches.tap_changer import TapChanger
from GridCal.Engine.Core.Devices.editable_device import DeviceType


class Transformer2W(ParentBranch):
    """
    The **Branch** class represents the connections between nodes (i.e.
    :ref:`buses<bus>`) in **GridCal**. A branch is an element (cable, line, capacitor,
    transformer, etc.) with an electrical impedance. The basic **Branch** class
    includes basic electrical attributes for most passive elements, but other device
    types may be passed to the **Branch** constructor to configure it as a specific
    type.

    For example, a transformer may be created with the following code:

    .. code:: ipython3

        from GridCal.Engine.Core.multi_circuit import MultiCircuit
        from GridCal.Engine.Core.Devices import *
        from GridCal.Engine.Core.Devices.types import *

        # Create grid
        grid = MultiCircuit()

        # Create buses
        POI = Bus(name="POI",
                  vnom=100, #kV
                  is_slack=True)
        grid.add_bus(POI)

        B_C3 = Bus(name="B_C3",
                   vnom=10) #kV
        grid.add_bus(B_C3)

        # Create transformer types
        SS = TransformerType(name="SS",
                             hv_nominal_voltage=100, # kV
                             lv_nominal_voltage=10, # kV
                             nominal_power=100, # MVA
                             copper_losses=10000, # kW
                             iron_losses=125, # kW
                             no_load_current=0.5, # %
                             short_circuit_voltage=8) # %
        grid.add_transformer_type(SS)

        # Create transformer
        X_C3 = Branch(bus_from=POI,
                      bus_to=B_C3,
                      name="X_C3",
                      branch_type=BranchType.Transformer,
                      template=SS,
                      )

        # Add transformer to grid
        grid.add_branch(X_C3)

    Refer to the :class:`GridCal.Engine.Devices.branch.TapChanger` class for an example
    using a voltage regulator.

    Arguments:

        **bus_from** (:ref:`Bus`): "From" :ref:`bus<Bus>` object

        **bus_to** (:ref:`Bus`): "To" :ref:`bus<Bus>` object

        **name** (str, "Branch"): Name of the branch

        **r** (float, 1e-20): Branch resistance in per unit

        **x** (float, 1e-20): Branch reactance in per unit

        **g** (float, 1e-20): Branch shunt conductance in per unit

        **b** (float, 1e-20): Branch shunt susceptance in per unit

        **rate** (float, 1.0): Branch rate in MVA

        **tap** (float, 1.0): Branch tap module

        **shift_angle** (int, 0): Tap shift angle in radians

        **active** (bool, True): Is the branch active?

        **tolerance** (float, 0): Tolerance specified for the branch impedance in %

        **mttf** (float, 0.0): Mean time to failure in hours

        **mttr** (float, 0.0): Mean time to recovery in hours

        **r_fault** (float, 0.0): Mid-line fault resistance in per unit (SC only)

        **x_fault** (float, 0.0): Mid-line fault reactance in per unit (SC only)

        **fault_pos** (float, 0.0): Mid-line fault position in per unit (0.0 = `bus_from`, 0.5 = middle, 1.0 = `bus_to`)

        **branch_type** (BranchType, BranchType.Line): Device type enumeration (ex.: :class:`GridCal.Engine.Devices.transformer.TransformerType`)

        **length** (float, 0.0): Length of the branch in km

        **vset** (float, 1.0): Voltage set-point of the voltage controlled bus in per unit

        **temp_base** (float, 20.0): Base temperature at which `r` is measured in °C

        **temp_oper** (float, 20.0): Operating temperature in °C

        **alpha** (float, 0.0033): Thermal constant of the material in °C

        **bus_to_regulated** (bool, False): Is the `bus_to` voltage regulated by this branch?

        **template** (BranchTemplate, BranchTemplate()): Basic branch template
    """

    def __init__(self, bus_from: Bus = None, bus_to: Bus = None, HV=None, LV=None, name='Branch', idtag=None, code='',
                 r=1e-20, x=1e-20, g=1e-20, b=1e-20,
                 rate=1.0,
                 tap=1.0, tap_module_max=1.2, tap_module_min=0.5,
                 shift_angle=0.0, theta_max=6.28, theta_min=-6.28,
                 active=True, tolerance=0, cost=100.0,
                 mttf=0, mttr=0,
                 vset=1.0, Pset=0, bus_to_regulated=False,
                 temp_base=20, temp_oper=20, alpha=0.00330,
                 control_mode: TransformerControlType = TransformerControlType.fixed,
                 template: TransformerType = None,
                 rate_prof=None, Cost_prof=None, active_prof=None, temp_oper_prof=None,
                 tap_module_prof=None, angle_prof=None,
                 contingency_factor=1.0,
                 contingency_enabled=True, monitor_loading=True, contingency_factor_prof=None,
                 r0=1e-20, x0=1e-20, g0=1e-20, b0=1e-20,
                 r2=1e-20, x2=1e-20, g2=1e-20, b2=1e-20,
                 conn: WindingsConnection = WindingsConnection.GG,
                 capex=0, opex=0, build_status: BuildStatus = BuildStatus.Commissioned):

        ParentBranch.__init__(self,
                              name=name,
                              idtag=idtag,
                              code=code,
                              bus_from=bus_from,
                              bus_to=bus_to,
                              active=active,
                              active_prof=active_prof,
                              rate=rate,
                              rate_prof=rate_prof,
                              contingency_factor=contingency_factor,
                              contingency_factor_prof=contingency_factor_prof,
                              contingency_enabled=contingency_enabled,
                              monitor_loading=monitor_loading,
                              mttf=mttf,
                              mttr=mttr,
                              build_status=build_status,
                              capex=capex,
                              opex=opex,
                              Cost=cost,
                              Cost_prof=Cost_prof,
                              device_type=DeviceType.Transformer2WDevice,
                              branch_type=BranchType.Transformer)

        # set the high and low voltage values
        self.HV = 0
        self.LV = 0
        self.set_hv_and_lv(HV, LV)

        # List of measurements
        self.measurements = list()

        # branch impedance tolerance
        self.tolerance = tolerance

        # total impedance and admittance in p.u.
        self.R = r
        self.X = x
        self.G = g
        self.B = b

        self.R0 = r0
        self.X0 = x0
        self.G0 = g0
        self.B0 = b0

        self.R2 = r2
        self.X2 = x2
        self.G2 = g2
        self.B2 = b2

        self.conn = conn 

        # Conductor base and operating temperatures in ºC
        self.temp_base = temp_base
        self.temp_oper = temp_oper

        self.temp_oper_prof = temp_oper_prof

        # Conductor thermal constant (1/ºC)
        self.alpha = alpha

        # tap changer object
        self.tap_changer = TapChanger()

        # Tap module
        if tap != 0:
            self.tap_module = tap
            self.tap_changer.set_tap(self.tap_module)
        else:
            self.tap_module = self.tap_changer.get_tap()

        self.tap_module_prof = tap_module_prof

        # Tap angle
        self.angle = shift_angle
        self.angle_prof = angle_prof

        self.tap_module_max = tap_module_max
        self.tap_module_min = tap_module_min
        self.angle_max = theta_max
        self.angle_min = theta_min

        # type template
        self.template = template

        self.vset = vset
        self.Pset = Pset

        self.control_mode = control_mode

        self.bus_to_regulated = bus_to_regulated

        if bus_to_regulated and self.control_mode == TransformerControlType.fixed:
            print(self.name, self.idtag, 'Overriding to V controller')
            self.control_mode = TransformerControlType.Vt

        # converter for enumerations
        self.conv = {'branch': BranchType.Branch,
                     'line': BranchType.Line,
                     'transformer': BranchType.Transformer,
                     'switch': BranchType.Switch,
                     'reactance': BranchType.Reactance}

        self.inv_conv = {val: key for key, val in self.conv.items()}

        self.register(key='HV', units='kV', tpe=float, definition='High voltage rating')
        self.register(key='LV', units='kV', tpe=float, definition='Low voltage rating')

        self.register(key='R', units='p.u.', tpe=float, definition='Total positive sequence resistance.')
        self.register(key='X', units='p.u.', tpe=float, definition='Total positive sequence reactance.')
        self.register(key='G', units='p.u.', tpe=float, definition='Total positive sequence shunt conductance.')
        self.register(key='B', units='p.u.', tpe=float, definition='Total positive sequence shunt susceptance.')
        self.register(key='R0', units='p.u.', tpe=float, definition='Total zero sequence resistance.')
        self.register(key='X0', units='p.u.', tpe=float, definition='Total zero sequence reactance.')
        self.register(key='G0', units='p.u.', tpe=float, definition='Total zero sequence shunt conductance.')
        self.register(key='B0', units='p.u.', tpe=float, definition='Total zero sequence shunt susceptance.')
        self.register(key='R2', units='p.u.', tpe=float, definition='Total negative sequence resistance.')
        self.register(key='X2', units='p.u.', tpe=float, definition='Total negative sequence reactance.')
        self.register(key='G2', units='p.u.', tpe=float, definition='Total negative sequence shunt conductance.')
        self.register(key='B2', units='p.u.', tpe=float, definition='Total negative sequence shunt susceptance.')
        self.register(key='conn', units='', tpe=WindingsConnection,
                      definition='Windings connection (from, to):G: grounded starS: ungrounded starD: delta')
        self.register(key='tolerance', units='%', tpe=float,
                      definition='Tolerance expected for the impedance values7% is expected for transformers0% for lines.')
        self.register(key='tap_module', units='', tpe=float, definition='Tap changer module, it a value close to 1.0',
                      profile_name='tap_module_prof')
        self.register(key='tap_module_max', units='', tpe=float, definition='Tap changer module max value')
        self.register(key='tap_module_min', units='', tpe=float, definition='Tap changer module min value')
        self.register(key='angle', units='rad', tpe=float, definition='Angle shift of the tap changer.',
                      profile_name='angle_prof')
        self.register(key='angle_max', units='rad', tpe=float, definition='Max angle.')
        self.register(key='angle_min', units='rad', tpe=float, definition='Min angle.')
        self.register(key='control_mode', units='', tpe=TransformerControlType,
                      definition='Control type of the transformer')
        self.register(key='vset', units='p.u.', tpe=float,
                      definition='Objective voltage at the "to" side of the bus when regulating the tap.')
        self.register(key='Pset', units='p.u.', tpe=float,
                      definition='Objective power at the "from" side of when regulating the angle.')
        self.register(key='temp_base', units='ºC', tpe=float, definition='Base temperature at which R was measured.')
        self.register(key='temp_oper', units='ºC', tpe=float, definition='Operation temperature to modify R.',
                      profile_name='temp_oper_prof')
        self.register(key='alpha', units='1/ºC', tpe=float,
                      definition='Thermal coefficient to modify R,around a reference temperatureusing a linear '
                                 'approximation.For example:Copper @ 20ºC: 0.004041,Copper @ 75ºC: 0.00323,'
                                 'Annealed copper @ 20ºC: 0.00393,Aluminum @ 20ºC: 0.004308,Aluminum @ 75ºC: 0.00330')
        self.register(key='template', units='', tpe=DeviceType.TransformerTypeDevice, definition='', editable=False)

    def set_hv_and_lv(self, HV, LV):
        """
        set the high and low voltage values
        :param HV: higher voltage value (kV)
        :param LV: lower voltage value (kV)
        """
        if self.bus_from is not None:
            vh = max(self.bus_from.Vnom, self.bus_to.Vnom)
            vl = min(self.bus_from.Vnom, self.bus_to.Vnom)
        else:
            vh = 1.0
            vl = 1.0

        if HV is None:
            self.HV = vh
        else:
            self.HV = HV

        if LV is None:
            self.LV = vl
        else:
            self.LV = LV

    @property
    def R_corrected(self):
        """
        Returns a temperature corrected resistance based on a formula provided by:
        NFPA 70-2005, National Electrical Code, Table 8, footnote #2; and
        https://en.wikipedia.org/wiki/Electrical_resistivity_and_conductivity#Linear_approximation
        (version of 2019-01-03 at 15:20 EST).
        """
        return self.R * (1 + self.alpha * (self.temp_oper - self.temp_base))

    def change_base(self, Sbase_old, Sbase_new):

        b = Sbase_new / Sbase_old

        self.R *= b
        self.X *= b
        self.G *= b
        self.B *= b

    def get_weight(self):
        return np.sqrt(self.R * self.R + self.X * self.X)

    def branch_type_converter(self, val_string):
        """
        function to convert the branch type string into the BranchType
        :param val_string:
        :return: branch type conversion
        """
        return self.conv[val_string.lower()]

    def copy(self, bus_dict=None):
        """
        Returns a copy of the branch
        @return: A new  with the same content as this
        """

        if bus_dict is None:
            f = self.bus_from
            t = self.bus_to
        else:
            f = bus_dict[self.bus_from]
            t = bus_dict[self.bus_to]

        # z_series = complex(self.R, self.X)
        # y_shunt = complex(self.G, self.B)
        b = Transformer2W(bus_from=f,
                          bus_to=t,
                          name=self.name,
                          r=self.R,
                          x=self.X,
                          g=self.G,
                          b=self.B,
                          rate=self.rate,
                          tap=self.tap_module,
                          shift_angle=self.angle,
                          active=self.active,
                          mttf=self.mttf,
                          mttr=self.mttr,
                          bus_to_regulated=self.bus_to_regulated,
                          vset=self.vset,
                          temp_base=self.temp_base,
                          temp_oper=self.temp_oper,
                          alpha=self.alpha,
                          template=self.template,
                          opex=self.opex,
                          capex=self.capex)

        b.measurements = self.measurements

        b.active_prof = self.active_prof
        b.rate_prof = self.rate_prof
        b.Cost_prof = self.Cost_prof

        return b

    def flip(self):
        """
        Change the terminals' positions
        """
        F, T = self.bus_from, self.bus_to
        self.bus_to, self.bus_from = F, T

    def tap_up(self):
        """
        Move the tap changer one position up
        """
        self.tap_changer.tap_up()
        self.tap_module = self.tap_changer.get_tap()

    def tap_down(self):
        """
        Move the tap changer one position up
        """
        self.tap_changer.tap_down()
        self.tap_module = self.tap_changer.get_tap()

    def apply_tap_changer(self, tap_changer: TapChanger):
        """
        Apply a new tap changer

        Argument:

            **tap_changer** (:class:`GridCal.Engine.Devices.branch.TapChanger`): Tap changer object

        """
        self.tap_changer = tap_changer

        if self.tap_module != 0:
            self.tap_changer.set_tap(self.tap_module)
        else:
            self.tap_module = self.tap_changer.get_tap()

    def get_sorted_buses_voltages(self):
        """
        GEt the sorted bus voltages
        :return: high voltage, low voltage
        """
        bus_f_v = self.bus_from.Vnom
        bus_t_v = self.bus_to.Vnom
        if bus_f_v > bus_t_v:
            return bus_f_v, bus_t_v
        else:
            return bus_t_v, bus_f_v

    def get_max_bus_nominal_voltage(self):
        return max(self.bus_from.Vnom, self.bus_to.Vnom)

    def get_min_bus_nominal_voltage(self):
        return min(self.bus_from.Vnom, self.bus_to.Vnom)

    def get_from_to_nominal_voltages(self):

        bus_f_v = self.bus_from.Vnom
        bus_t_v = self.bus_to.Vnom

        dhf = abs(self.HV - bus_f_v)
        dht = abs(self.HV - bus_t_v)

        if dhf < dht:
            # the HV side is on the from side
            tpe_f_v = self.HV
            tpe_t_v = self.LV
        else:
            # the HV side is on the to side
            tpe_t_v = self.HV
            tpe_f_v = self.LV

        return tpe_f_v, tpe_t_v

    def get_virtual_taps(self):
        """
        Get the branch virtual taps

        The virtual taps generate when a transformer nominal winding voltage differs
        from the bus nominal voltage.

        Returns:

            **tap_f** (float, 1.0): Virtual tap at the *from* side

            **tap_t** (float, 1.0): Virtual tap at the *to* side

        """
        # resolve how the transformer is actually connected and set the virtual taps
        bus_f_v = self.bus_from.Vnom
        bus_t_v = self.bus_to.Vnom

        # obtain the nominal voltages at the from and to sides
        tpe_f_v, tpe_t_v = self.get_from_to_nominal_voltages()

        tap_f = tpe_f_v / bus_f_v if bus_f_v > 0 else 1.0
        tap_t = tpe_t_v / bus_t_v if bus_t_v > 0 else 1.0

        if tap_f == 0.0:
            tap_f = 1.0

        if tap_t == 0.0:
            tap_t = 1.0

        return tap_f, tap_t

    def apply_template(self, obj: TransformerType, Sbase, logger=Logger()):
        """
        Apply a branch template to this object

        Arguments:
            **obj**: TransformerType or Tower object
            **Sbase** (float): circuit base power in MVA
            **logger** (list, []): Log list
        """
        if isinstance(obj, TransformerType):

            VH, VL = self.get_sorted_buses_voltages()

            # get the transformer impedance in the base of the transformer
            z_series, y_shunt = obj.get_impedances(VH=VH, VL=VL, Sbase=Sbase)

            self.R = np.round(z_series.real, 6)
            self.X = np.round(z_series.imag, 6)
            self.G = np.round(y_shunt.real, 6)
            self.B = np.round(y_shunt.imag, 6)

            self.rate = obj.rating

            self.HV = obj.HV
            self.LV = obj.LV

            if self.template is not None:
                if obj != self.template:
                    self.template = obj
                else:
                    logger.add_error('Template not recognised', self.name)
            else:
                self.template = obj

    def get_save_data(self):
        """
        Return the data that matches the edit_headers
        :return:
        """
        data = list()
        for name, properties in self.editable_headers.items():
            obj = getattr(self, name)

            if properties.tpe == BranchType:
                obj = self.branch_type.value

            elif properties.tpe == DeviceType.BusDevice:
                obj = obj.idtag

            elif properties.tpe == DeviceType.TransformerTypeDevice:
                if obj is None:
                    obj = ''
                else:
                    obj = obj.idtag

            elif properties.tpe not in [str, float, int, bool]:
                obj = str(obj)

            data.append(obj)
        return data

    def get_properties_dict(self, version=3):
        """
        Get json dictionary
        :return:
        """
        # get the virtual taps
        tap_f, tap_t = self.get_virtual_taps()

        # get the nominal voltages
        v_from, v_to = self.get_from_to_nominal_voltages()

        '''
        TransformerControlType(Enum):
        fixed = '0:Fixed'
        Pt = '1:Pt'
        Qt = '2:Qt'
        PtQt = '3:Pt+Qt'
        Vt = '4:Vt'
        PtVt = '5:Pt+Vt'
        
        '''
        control_modes = {TransformerControlType.fixed: 0,
                         TransformerControlType.Vt: 1,
                         TransformerControlType.Pt: 2,
                         TransformerControlType.PtVt: 3,
                         TransformerControlType.Qt: 4,
                         TransformerControlType.PtQt: 5}
        if version == 2:
            d = {'id': self.idtag,
                 'type': 'transformer',
                 'phases': 'ps',
                 'name': self.name,
                 'name_code': self.code,
                 'bus_from': self.bus_from.idtag,
                 'bus_to': self.bus_to.idtag,
                 'active': self.active,
                 'rate': self.rate,
                 'Vnomf': v_from,
                 'Vnomt': v_to,
                 'hv': self.HV,
                 'lv': self.LV,
                 'r': self.R,
                 'x': self.X,
                 'g': self.G,
                 'b': self.B,
                 'tap_module': self.tap_module,
                 'min_tap_module': self.tap_module_min,
                 'max_tap_module': self.tap_module_max,
                 'id_tap_module_table': "",

                 'tap_angle': self.angle,
                 'min_tap_angle': self.angle_min,
                 'max_tap_angle': self.angle_max,
                 'id_tap_angle_table': "",

                 'control_mode': control_modes[self.control_mode],

                 # 'min_tap_position': self.tap_changer.min_tap,
                 # 'max_tap_position': self.tap_changer.max_tap,
                 # 'tap_inc_reg_down': self.tap_changer.inc_reg_down,
                 # 'tap_inc_reg_up': self.tap_changer.inc_reg_up,
                 # 'virtual_tap_from': tap_f,
                 # 'virtual_tap_to': tap_t,
                 # 'bus_to_regulated': self.bus_to_regulated,

                 'vset': self.vset,
                 'pset': self.Pset,

                 'base_temperature': self.temp_base,
                 'operational_temperature': self.temp_oper,
                 'alpha': self.alpha
                 }

        elif version == 3:
            d = {'id': self.idtag,
                 'type': 'transformer',
                 'phases': 'ps',
                 'name': self.name,
                 'name_code': self.code,
                 'bus_from': self.bus_from.idtag,
                 'bus_to': self.bus_to.idtag,
                 'active': self.active,
                 'rate': self.rate,
                 'contingency_factor1': self.contingency_factor,
                 'contingency_factor2': self.contingency_factor,
                 'contingency_factor3': self.contingency_factor,

                 'Vnomf': v_from,
                 'Vnomt': v_to,
                 'hv': self.HV,
                 'lv': self.LV,
                 'r': self.R,
                 'x': self.X,
                 'g': self.G,
                 'b': self.B,
                 'tap_module': self.tap_module,
                 'min_tap_module': self.tap_module_min,
                 'max_tap_module': self.tap_module_max,
                 'id_tap_module_table': "",

                 'tap_angle': self.angle,
                 'min_tap_angle': self.angle_min,
                 'max_tap_angle': self.angle_max,
                 'id_tap_angle_table': "",

                 'control_mode': control_modes[self.control_mode],

                 # 'min_tap_position': self.tap_changer.min_tap,
                 # 'max_tap_position': self.tap_changer.max_tap,
                 # 'tap_inc_reg_down': self.tap_changer.inc_reg_down,
                 # 'tap_inc_reg_up': self.tap_changer.inc_reg_up,
                 # 'virtual_tap_from': tap_f,
                 # 'virtual_tap_to': tap_t,
                 # 'bus_to_regulated': self.bus_to_regulated,

                 'vset': self.vset,
                 'pset': self.Pset,

                 'base_temperature': self.temp_base,
                 'operational_temperature': self.temp_oper,
                 'alpha': self.alpha,

                 'overload_cost': self.Cost,
                 'capex': self.capex,
                 'opex': self.opex,
                 'build_status': str(self.build_status.value).lower(),
                 }
        else:
            d = dict()

        return d

    def get_profiles_dict(self, version=3):
        """

        :return:
        """
        if self.active_prof is not None:
            active_prof = self.active_prof.tolist()
            rate_prof = self.rate_prof.tolist()
        else:
            active_prof = list()
            rate_prof = list()

        return {'id': self.idtag,
                'active': active_prof,
                'rate': rate_prof}

    def get_units_dict(self, version=3):
        """
        Get units of the values
        """
        return {'rate': 'MW',
                'r': 'p.u.',
                'x': 'p.u.',
                'b': 'p.u.',
                'g': 'p.u.',
                'base_temperature': 'ºC',
                'operational_temperature': 'ºC',
                'alpha': '1/ºC'}

    def plot_profiles(self, time_series=None, my_index=0, show_fig=True):
        """
        Plot the time series results of this object
        :param time_series: TimeSeries Instance
        :param my_index: index of this object in the simulation
        :param show_fig: Show the figure?
        """

        if time_series is not None:
            fig = plt.figure(figsize=(12, 8))

            ax_1 = fig.add_subplot(211)
            ax_2 = fig.add_subplot(212, sharex=ax_1)

            x = time_series.results.time_array

            # loading
            y = time_series.results.loading.real * 100.0
            df = pd.DataFrame(data=y[:, my_index], index=x, columns=[self.name])
            ax_1.set_title('Loading', fontsize=14)
            ax_1.set_ylabel('Loading [%]', fontsize=11)
            df.plot(ax=ax_1)

            # losses
            y = np.abs(time_series.results.losses)
            df = pd.DataFrame(data=y[:, my_index], index=x, columns=[self.name])
            ax_2.set_title('Losses', fontsize=14)
            ax_2.set_ylabel('Losses [MVA]', fontsize=11)
            df.plot(ax=ax_2)

            plt.legend()
            fig.suptitle(self.name, fontsize=20)

        if show_fig:
            plt.show()

    def get_coordinates(self):
        """
        Get the branch defining coordinates
        """
        return [self.bus_from.get_coordinates(), self.bus_to.get_coordinates()]

    def delete_virtual_taps(self):
        """
        Set the HV and LV parameters such that any virtual tap is null
        """
        self.HV = max(self.bus_from.Vnom, self.bus_to.Vnom)
        self.LV = min(self.bus_from.Vnom, self.bus_to.Vnom)

    def fix_inconsistencies(self, logger: Logger, maximum_difference=0.1):
        """
        Fix the inconsistencies
        :param logger:
        :param maximum_difference: proportion to be under or above (i.e. Transformer HV=41.9, bus HV=45 41.9/45 = 0.93 -> 0.9 <= 0.93 <= 1.1, so its ok
        :return:
        """
        errors = False
        HV = max(self.bus_from.Vnom, self.bus_to.Vnom)
        LV = min(self.bus_from.Vnom, self.bus_to.Vnom)

        if self.LV > self.HV:
            logger.add_warning("HV > LV", self.name, self.HV, HV)
            self.HV, self.LV = self.LV, self.HV
            errors = True

        rHV = self.HV / HV
        rLV = self.LV / LV
        LB = 1 - maximum_difference
        UB = 1 + maximum_difference
        if not (LB <= rHV <= UB):
            logger.add_warning("Corrected transformer HV", self.name, self.HV, HV)
            self.HV = HV
            errors = True

        if not (LB <= rLV <= UB):
            logger.add_warning("Corrected transformer LV", self.name, self.LV, LV)
            self.LV = LV
            errors = True

        if self.R < 0.0:
            logger.add_warning("Corrected transformer R<0", self.name, self.R, -self.R)
            self.R = -self.R
            errors = True

        return errors
