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
from GridCal.Engine.IO.cim.cgmes_2_4_15.devices.base import Base
from GridCal.Engine.IO.cim.cgmes_2_4_15.cim_enums import cgmesProfile
from GridCal.Engine.IO.cim.cgmes_2_4_15.devices.topological_node import TopologicalNode
from GridCal.Engine.IO.base.units import UnitMultiplier, UnitSymbol


class SvVoltage(Base):

    def __init__(self, rdfid, tpe, resources=list(), class_replacements=dict()):
        """
        General CIM object container
        :param rdfid: RFID
        :param tpe: type of the object (class)
        """
        Base.__init__(self, rdfid='', tpe=tpe, resources=resources, class_replacements=class_replacements)

        self.TopologicalNode: TopologicalNode = None

        self.v: float = 0.0

        self.angle: float = 0.0

        self.register_property(name='TopologicalNode',
                               class_type=TopologicalNode,
                               multiplier=UnitMultiplier.none,
                               unit=UnitSymbol.none,
                               description="The state voltage associated with the topological node. Default: None",
                               profiles=[cgmesProfile.SV])

        self.register_property(name='v',
                               class_type=float,
                               multiplier=UnitMultiplier.k,
                               unit=UnitSymbol.V,
                               description=" The voltage magnitude of the topological node. Default: 0.0",
                               profiles=[cgmesProfile.SV])

        self.register_property(name='angle',
                               class_type=float,
                               multiplier=UnitMultiplier.none,
                               unit=UnitSymbol.deg,
                               description="The voltage angle of the topological node complex voltage with "
                                           "respect to system reference. Default: 0.0",
                               profiles=[cgmesProfile.SV])
