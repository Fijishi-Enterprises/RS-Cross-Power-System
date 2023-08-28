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


from GridCal.Engine.Core.Devices.editable_device import EditableDevice, DeviceType, GCProp
from GridCal.Engine.Core.Devices.Aggregation.investments_group import InvestmentsGroup


class Investment(EditableDevice):
    """
    The Contingency object is the container of all the

    Arguments:

        **name** (str, "Contingency"): Name of the contingency

        **type** (float, 10.0): Nominal voltage in kV

        **loads** (list, list()): List of contingency elements

    """

    def __init__(self, idtag=None, device_idtag=None, name="Investment", code='', CAPEX=0.0, OPEX=0.0,
                 group: InvestmentsGroup = None, comment: str = ""):
        """
        Contingency
        :param idtag: String. Element unique identifier
        :param name: String. Contingency name
        :param code: String. Contingency code name
        :param CAPEX: Float. Capital expenditures
        :param OPEX: Float. Operating expenditures
        :param group: ContingencyGroup. Contingency group
        :param comment: Comment
        """

        EditableDevice.__init__(
            self,
            idtag=idtag,
            code=code,
            active=True,
            name=name,
            device_type=DeviceType.InvestmentDevice
        )

        # Contingency type
        self.device_idtag = device_idtag
        self.CAPEX = CAPEX
        self.OPEX = OPEX
        self._group: InvestmentsGroup = group
        self.comment = comment

        self.register(key='device_idtag', units='', tpe=str, definition='Unique ID')
        self.register(key='CAPEX', units='M€', tpe=float,
                      definition='Capital expenditures. This is the initial investment.')
        self.register(key='OPEX', units='M€', tpe=float,
                      definition='Operation expenditures. Maintenance costs among other recurrent costs.')
        self.register(key='group', units='', tpe=DeviceType.InvestmentsGroupDevice, definition='Investment group')
        self.register(key='comment', units='', tpe=str, definition='Comments')

    @property
    def group(self):
        return self._group

    @group.setter
    def group(self, val: InvestmentsGroup):
        self._group = val

    @property
    def category(self):
        """
        Display the group category
        :return:
        """
        return self.group.category

    @category.setter
    def category(self, val):
        # self.group.category = val
        pass

    def get_properties_dict(self, version=3):
        """
        Get json dictionary
        :return:
        """

        return {
            'id': self.idtag,
            'name': self.name,
            'name_code': self.code,
            'group': self._group,
            'device_idtag': self.device_idtag,
            'CAPEX': self.CAPEX,
            'OPEX': self.OPEX,
            'comment': self.comment
        }
