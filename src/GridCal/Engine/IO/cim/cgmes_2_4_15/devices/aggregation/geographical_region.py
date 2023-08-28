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
from GridCal.Engine.IO.cim.cgmes_2_4_15.cim_enums import cgmesProfile
from GridCal.Engine.IO.cim.cgmes_2_4_15.devices.identified_object import IdentifiedObject
from GridCal.Engine.IO.base.units import UnitMultiplier, UnitSymbol
from GridCal.Engine.IO.cim.cgmes_2_4_15.devices.aggregation.modelling_authority import ModelingAuthority


class GeographicalRegion(IdentifiedObject):
    """

    """
    def __init__(self, rdfid, tpe):
        IdentifiedObject.__init__(self, rdfid, tpe)

        self.ModelingAuthority: Union[ModelingAuthority, None] = None
        self.sourcingActor: str = ''
        self.masUri: str = ''

        self.register_property(name='ModelingAuthority',
                               class_type=ModelingAuthority,
                               multiplier=UnitMultiplier.none,
                               unit=UnitSymbol.none,
                               description="",
                               profiles=[cgmesProfile.EQ])
        self.register_property(name='sourcingActor',
                               class_type=str,
                               multiplier=UnitMultiplier.none,
                               unit=UnitSymbol.none,
                               description="",
                               profiles=[cgmesProfile.EQ])
        self.register_property(name='masUri',
                               class_type=str,
                               multiplier=UnitMultiplier.none,
                               unit=UnitSymbol.none,
                               description="",
                               profiles=[cgmesProfile.EQ])

