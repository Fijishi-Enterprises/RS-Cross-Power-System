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
from GridCal.Engine.IO.cim.cgmes_2_4_15.devices.aggregation.operational_limit_set import OperationalLimitSet
from GridCal.Engine.IO.cim.cgmes_2_4_15.cim_enums import cgmesProfile
from GridCal.Engine.IO.cim.cgmes_2_4_15.devices.equipment_container import EquipmentContainer
from GridCal.Engine.IO.cim.cgmes_2_4_15.devices.injections.power_systems_resource import PowerSystemResource
from GridCal.Engine.IO.base.units import UnitMultiplier, UnitSymbol


class Equipment(PowerSystemResource):

    def __init__(self, rdfid, tpe):
        PowerSystemResource.__init__(self, rdfid, tpe)

        self.aggregate: bool = False
        self.EquipmentContainer: EquipmentContainer = None
        self.OperationalLimitSet: OperationalLimitSet = None

        self.register_property(name='aggregate',
                               class_type=bool,
                               multiplier=UnitMultiplier.none,
                               unit=UnitSymbol.none,
                               description="aggregate",
                               profiles=[cgmesProfile.EQ])

        self.register_property(name='EquipmentContainer',
                               class_type=EquipmentContainer,
                               multiplier=UnitMultiplier.none,
                               unit=UnitSymbol.none,
                               description="EquipmentContainer",
                               profiles=[cgmesProfile.EQ, cgmesProfile.EQ_BD, ])

        self.register_property(name='OperationalLimitSet',
                               class_type=OperationalLimitSet,
                               multiplier=UnitMultiplier.none,
                               unit=UnitSymbol.none,
                               description="OperationalLimitSet",
                               profiles=[cgmesProfile.EQ])
