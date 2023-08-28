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
from GridCal.Engine.IO.base.units import UnitMultiplier, UnitSymbol, Unit
from GridCal.Engine.IO.psse.devices.psse_object import PSSeObject
from GridCal.Engine.basic_structures import Logger
import GridCal.Engine.Core.Devices as dev


class PSSeSubstation(PSSeObject):

    def __init__(self):
        PSSeObject.__init__(self, "Substation")

        self.IS: int = 0
        self.NAME: str = ""
        self.LATI: float = 0.0
        self.LONG: float = 0.0
        self.SGR: float = 0.0

        self.register_property(property_name="IS",
                               rawx_key='isub',
                               class_type=int,
                               description="Substation number ",
                               min_value=1,
                               max_value=99999)

        self.register_property(property_name="NAME",
                               rawx_key='name',
                               class_type=str,
                               description="ubstation name, NAME may be up to forty characters and may contain "
                                           "any combi-nation of blanks, uppercase letters, numbers and special "
                                           "characters. NAME must beenclosed in single or double quotes if it "
                                           "contains any blanks or special characters.",
                               max_chars=40)

        self.register_property(property_name="LATI",
                               rawx_key='lati',
                               class_type=float,
                               description="Substation latitude in degrees",
                               min_value=-90,
                               max_chars=90,
                               unit=Unit.get_deg())

        self.register_property(property_name="LONG",
                               rawx_key='long',
                               class_type=float,
                               description="Substation longitude in degrees..",
                               min_value=-180,
                               max_chars=180,
                               unit=Unit.get_deg())

        self.register_property(property_name="SGR",
                               rawx_key='sgr',
                               class_type=float,
                               description="Substation grounding DC resistance in ohms.",
                               unit=Unit.get_ohm())

    def parse(self, data, version, logger: Logger):
        """

        :param data:
        :param version:
        :param logger:
        """

        if version >= 35:

            self.IS, self.NAME, self.LATI, self.LONG, self.SGR = data[0]

        else:
            logger.add_warning('Areas not defined for version', str(version))

    def get_raw_line(self, version):

        if version >= 35:
            return self.format_raw_line([self.IS, self.NAME, self.LATI, self.LONG, self.SGR])
        else:
            raise Exception('Areas not defined for version', str(version))

    def get_id(self) -> str:
        return str(self.IS)
