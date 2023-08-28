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
from GridCal.Engine.Core.Devices.Aggregation.contingency_group import ContingencyGroup


class Contingency(EditableDevice):
    """
    The Contingency object is the container of all the

    Arguments:

        **name** (str, "Contingency"): Name of the contingency

        **type** (float, 10.0): Nominal voltage in kV

        **loads** (list, list()): List of contingency elements

    """

    def __init__(self, idtag=None, device_idtag='', name="Contingency", code='', prop='active', value=0.0,
                 group: ContingencyGroup=None):
        """
        Contingency
        :param idtag: String. Element unique identifier
        :param name: String. Contingency name
        :param code: String. Contingency code name
        :param prop: String. Property to modify when contingency is triggered out
        :param value: Float. Property value to apply when contingency happens
        :param group: ContingencyGroup. Contingency group
        """

        EditableDevice.__init__(
            self,
            idtag=idtag,
            code=code,
            active=True,
            name=name,
            device_type=DeviceType.ContingencyDevice
        )

        # Contingency type
        self.device_idtag = device_idtag
        self._prop = prop
        self._value = value
        self._group: ContingencyGroup = group

        self.register(key='idtag', units='', tpe=str, definition='Unique ID', editable=False)
        self.register(key='device_idtag', units='', tpe=str, definition='Unique ID', editable=False)
        self.register(key='prop', units='', tpe=str, definition='Name of the object property to change')
        self.register(key='value', units='', tpe=float, definition='Property value')
        self.register(key='group', units='', tpe=DeviceType.ContingencyGroupDevice, definition='Contingency group')

    @property
    def name(self) -> str:
        return self._name

    @name.setter
    def name(self, val: str):
        self._name = val
        if self.graphic_obj is not None:
            self.graphic_obj.set_label(self._name)

    @property
    def prop(self) -> str:
        """
        Property to modify when contingency is triggered out
        :return:
        """
        return self._prop

    @prop.setter
    def prop(self, val: str):
        if val in ['active']:
            self._prop = val

    @property
    def value(self) -> float:
        """
        Property value to apply when contingency happens
        :return:
        """
        return self._value

    @value.setter
    def value(self, val: float):
        self._value = val

    @property
    def group(self) -> ContingencyGroup:
        """
        Contingency group
        :return:
        """
        return self._group

    @group.setter
    def group(self, val: ContingencyGroup):
        self._group = val

    @property
    def category(self):
        """

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
            'group': self._group.idtag,
            'prop': self.prop,
            'value': self.value
        }
