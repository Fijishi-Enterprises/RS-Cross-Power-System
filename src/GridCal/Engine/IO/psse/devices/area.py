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


class PSSeArea(PSSeObject):

    def __init__(self):
        PSSeObject.__init__(self, "Area")

        self.I: int = -1
        self.ISW: int = 0
        self.PDES: float = 0.0
        self.PTOL: float = 0.0
        self.ARNAME: str = ''

        self.register_property(property_name="I",
                               rawx_key='iarea',
                               class_type=int,
                               description="Area number",
                               min_value=1,
                               max_value=9999)

        self.register_property(property_name="ISW",
                               rawx_key='isw',
                               class_type=int,
                               description="Bus number, or extended bus name enclosed in single quotes "
                                           "(refer to Extended Bus Names), of the area slack bus for area "
                                           "interchange control. The bus must be a generator (Type 2) bus in "
                                           "the specified area. Any area containing a system swing bus (Type 3) "
                                           "must have either that swing bus or a bus number of zero specified "
                                           "for its area slack bus number. Any area with an area slack bus number "
                                           "of zero is considered a floating area by the area interchange control "
                                           "option of the power flow solution activities.")

        self.register_property(property_name="PDES",
                               rawx_key='pdes',
                               class_type=float,
                               description="Desired net interchange leaving the area (export); entered in MW. "
                                           "PDES must be specified such that is consistent with the area interchange "
                                           "definition implied by the area interchange control code (tie lines only, "
                                           "or tie lines and loads) to be specified during power flow solutions "
                                           "(refer to Section 6.3.20, Automatic Adjustments and Area Interchange "
                                           "Control).")

        self.register_property(property_name="PTOL",
                               rawx_key='ptol',
                               class_type=float,
                               description="Interchange tolerance bandwidth; entered in MW.")

        self.register_property(property_name="ARNAME",
                               rawx_key='arname',
                               class_type=str,
                               description="Alphanumeric identifier assigned to area I. "
                                           "ARNAME may be up to twelve characters and may contain any combination "
                                           "of blanks, uppercase letters, numbers and special characters. "
                                           "ARNAME must be enclosed in single or double quotes if it contains "
                                           "any blanks or special characters.",
                               max_chars=12)

    def parse(self, data, version, logger: Logger):
        """

        :param data:
        :param version:
        :param logger:
        """

        self.I = -1

        self.ARNAME = ''

        if version >= 29:
            # I, ISW, PDES, PTOL, 'ARNAME'
            self.I, self.ISW, self.PDES, self.PTOL, self.ARNAME = data[0]

            self.ARNAME = self.ARNAME.replace("'", "").strip()
        else:
            logger.add_warning('Areas not defined for version', str(version))

    def get_raw_line(self, version):

        if version >= 29:
            return self.format_raw_line([self.I, self.ISW, self.PDES, self.PTOL, self.ARNAME])
        else:
            raise Exception('Areas not defined for version', str(version))

    def get_id(self) -> str:
        return str(self.I)

    def get_seed(self):
        return "_CA_{}".format(self.I)
