# GridCal
# Copyright (C) 2022 Santiago Peñate Vera
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
from matplotlib import pyplot as plt
from GridCal.Engine.Core.Devices.editable_device import EditableDevice, DeviceType, GCProp
from GridCal.Engine.Core.Devices.enumerations import BuildStatus


class StaticGenerator(EditableDevice):
    """
    Arguments:

        **name** (str, "StaticGen"): Name of the static generator

        **P** (float, 0.0): Active power in MW

        **Q** (float, 0.0): Reactive power in MVAr

        **P_prof** (DataFrame, None): Pandas DataFrame with the active power profile in MW

        **Q_prof** (DataFrame, None): Pandas DataFrame with the reactive power profile in MVAr

        **active** (bool, True): Is the static generator active?

        **mttf** (float, 0.0): Mean time to failure in hours

        **mttr** (float, 0.0): Mean time to recovery in hours

    """

    def __init__(self, name='StaticGen', idtag=None, code='', P=0.0, Q=0.0, P_prof=None, Q_prof=None, active=True,
                 mttf=0.0, mttr=0.0, cost=1200.0,
                 capex=0, opex=0, build_status: BuildStatus = BuildStatus.Commissioned):

        EditableDevice.__init__(self,
                                name=name,
                                idtag=idtag,
                                code=code,
                                active=active,
                                device_type=DeviceType.StaticGeneratorDevice)

        self.bus = None

        self.active_prof = None

        self.mttf = mttf

        self.mttr = mttr

        # Power (MW + jMVAr)
        self.P = P
        self.Q = Q

        # power profile for this load
        self.P_prof = P_prof
        self.Q_prof = Q_prof

        self.Cost = cost
        self.Cost_prof = None

        self.capex = capex

        self.opex = opex

        self.build_status = build_status

        self.register(key='bus', units='', tpe=DeviceType.BusDevice, definition='', editable=False)
        self.register(key='active', units='', tpe=bool, definition='', profile_name='active_prof')
        self.register(key='P', units='MW', tpe=float, definition='Active power', profile_name='P_prof')
        self.register(key='Q', units='MVAr', tpe=float, definition='Reactive power', profile_name='Q_prof')
        self.register(key='mttf', units='h', tpe=float, definition='Mean time to failure')
        self.register(key='mttr', units='h', tpe=float, definition='Mean time to recovery')
        self.register(key='Cost', units='e/MWh', tpe=float, definition='Cost of not served energy. Used in OPF.',
                      profile_name='Cost_prof')
        self.register(key='capex', units='e/MW', tpe=float,
                      definition='Cost of investment. Used in expansion planning.')
        self.register(key='opex', units='e/MWh', tpe=float, definition='Cost of operation. Used in expansion planning.')
        self.register(key='build_status', units='', tpe=BuildStatus,
                      definition='Branch build status. Used in expansion planning.')

    def copy(self):
        """
        Deep copy of this object
        :return:
        """
        return StaticGenerator(name=self.name,
                               P=self.P,
                               Q=self.Q,
                               P_prof=self.P_prof,
                               Q_prof=self.Q_prof,
                               mttf=self.mttf,
                               mttr=self.mttr)

    def get_properties_dict(self, version=3):
        """
        Get json dictionary
        :return:
        """

        data = {'id': self.idtag,
                'phases': 'ps',
                'name': self.name,
                'name_code': self.code,
                'bus': self.bus.idtag,
                'active': self.active,
                'P': self.P,
                'Q': self.Q,
                'capex': self.capex,
                'opex': self.opex,
                'build_status': str(self.build_status.value).lower(),
                'technology': "",
                'Cost': self.Cost
                }

        if self.active_prof is not None:
            data['active_profile'] = self.active_prof.tolist()
            data['P_prof'] = self.P_prof.tolist()
            data['Q_prof'] = self.Q_prof.tolist()
            data['Cost_prof'] = self.Cost_prof.tolist()

        return data

    def get_profiles_dict(self, version=3):
        """

        :return:
        """

        if self.active_prof is not None:
            active_profile = self.active_prof.tolist()
            P_prof = self.P_prof.tolist()
            Q_prof = self.Q_prof.tolist()
        else:
            active_profile = list()
            P_prof = list()
            Q_prof = list()

        return {'id': self.idtag,
                'active': active_profile,
                'p': P_prof,
                'q': Q_prof}

    def get_units_dict(self, version=3):
        """
        Get units of the values
        """
        return {'p': 'MW',
                'q': 'MVAr'}

    def plot_profiles(self, time=None, show_fig=True):
        """
        Plot the time series results of this object
        :param time: array of time values
        :param show_fig: Show the figure?
        """

        if time is not None:
            fig = plt.figure(figsize=(12, 8))

            ax_1 = fig.add_subplot(211)
            ax_2 = fig.add_subplot(212, sharex=ax_1)

            # P
            y = self.P_prof
            df = pd.DataFrame(data=y, index=time, columns=[self.name])
            ax_1.set_title('Active power', fontsize=14)
            ax_1.set_ylabel('MW', fontsize=11)
            df.plot(ax=ax_1)

            # Q
            y = self.Q_prof
            df = pd.DataFrame(data=y, index=time, columns=[self.name])
            ax_2.set_title('Reactive power', fontsize=14)
            ax_2.set_ylabel('MVAr', fontsize=11)
            df.plot(ax=ax_2)

            plt.legend()
            fig.suptitle(self.name, fontsize=20)

            if show_fig:
                plt.show()