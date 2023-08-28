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


from typing import Union
from GridCal.Engine.Core.Devices.editable_device import EditableDevice, DeviceType


class GenericAreaGroup(EditableDevice):

    def __init__(self, name='', code='', idtag: Union[str, None] = None, device_type=DeviceType.GenericArea, latitude=0.0, longitude=0.0):
        """

        :param name:
        :param idtag:
        :param device_type:
        :param latitude:
        :param longitude:
        """
        EditableDevice.__init__(self,
                                name=name,
                                code=code,
                                idtag=idtag,
                                active=True,
                                device_type=device_type)

        self.latitude = latitude
        self.longitude = longitude

        self.register(key='longitude', units='deg', tpe=float, definition='longitude of the bus.', profile_name='',
                      editable=False)
        self.register(key='latitude', units='deg', tpe=float, definition='latitude of the bus.', profile_name='',
                      editable=False)

    def get_properties_dict(self, version=3):

        data = {'id': self.idtag,
                'name': self.name,
                'code': self.code
                }
        return data

    def get_profiles_dict(self, version=3):
        data = {'id': self.idtag}
        return data

    def get_units_dict(self, version=3):
        data = {}
        return data


class Area(GenericAreaGroup):

    def __init__(self, name='Area', idtag=None, code='', latitude=0.0, longitude=0.0):
        """

        :param name:
        :param idtag:
        :param latitude:
        :param longitude:
        """
        GenericAreaGroup.__init__(self,
                                  name=name,
                                  idtag=idtag,
                                  code=code,
                                  device_type=DeviceType.AreaDevice,
                                  latitude=latitude,
                                  longitude=longitude)
