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


class PSSeNode(PSSeObject):

    def __init__(self):
        PSSeObject.__init__(self, "node")

        self.ISUB: int = 0
        self.NI: int = 0
        self.NAME: str = ''
        self.I: int = 0
        self.STATUS: int = 0
        self.VM: float = 0.0
        self.VA: float = 0.0

        self.register_property(property_name="ISUB",
                               rawx_key='isub',
                               class_type=int,
                               description="Substation number",
                               min_value=1,
                               max_value=99999)

        self.register_property(property_name="NI",
                               rawx_key='inode',
                               class_type=int,
                               description="Node number",
                               min_value=1,
                               max_value=9999)

        self.register_property(property_name="NAME",
                               rawx_key='name',
                               class_type=str,
                               description="Node name.",
                               max_chars=12)

        self.register_property(property_name="I",
                               rawx_key='ibus',
                               class_type=int,
                               description="Bus number",
                               min_value=1,
                               max_value=999997)

        self.register_property(property_name="STATUS",
                               rawx_key="stat",
                               class_type=int,
                               description="Switch status, 1: closed, 0: open")

        self.register_property(property_name="VM",
                               rawx_key="vm",
                               class_type=float,
                               description="Bus voltage magnitude",
                               unit=Unit(UnitMultiplier.none, UnitSymbol.pu),
                               min_value=0.0,
                               max_value=2.0)

        self.register_property(property_name="VA",
                               rawx_key="va",
                               class_type=float,
                               description="Bus voltage angle",
                               unit=Unit(UnitMultiplier.none, UnitSymbol.deg),
                               min_value=0.0,
                               max_value=360.0)

    def parse(self, data, version, logger: Logger):
        """

        :param data:
        :param version:
        :param logger:
        """

        if version >= 35:
            # I, ISW, PDES, PTOL, 'ARNAME'
            self.ISUB, self.NI, self.NAME, self.I, self.STATUS, self.VM, self.VA = data[0]

            self.NAME = self.NAME.replace("'", "").strip()
        else:
            logger.add_warning('Areas not defined for version', str(version))

    def get_raw_line(self, version):

        if version >= 29:
            return self.format_raw_line([self.ISUB, self.NI, self.NAME, self.I, self.STATUS, self.VM, self.VA])
        else:
            raise Exception('Areas not defined for version', str(version))

    def get_id(self) -> str:
        return "{0}_{1}".format(self.ISUB, self.NI)
