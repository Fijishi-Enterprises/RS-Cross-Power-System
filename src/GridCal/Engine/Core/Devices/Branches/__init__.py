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


from GridCal.Engine.Core.Devices.Branches.branch import Branch, BranchTemplate
from GridCal.Engine.Core.Devices.Branches.dc_line import DcLine
from GridCal.Engine.Core.Devices.Branches.line import Line
from GridCal.Engine.Core.Devices.Branches.hvdc_line import HvdcLine
from GridCal.Engine.Core.Devices.Branches.switch import Switch
from GridCal.Engine.Core.Devices.Branches.transformer import Transformer2W
from GridCal.Engine.Core.Devices.Branches.transformer3w import Transformer3W
from GridCal.Engine.Core.Devices.Branches.upfc import UPFC
from GridCal.Engine.Core.Devices.Branches.vsc import VSC
from GridCal.Engine.Core.Devices.Branches.winding import Winding
from GridCal.Engine.Core.Devices.Branches.wire import Wire
from GridCal.Engine.Core.Devices.Branches.tap_changer import TapChanger

from GridCal.Engine.Core.Devices.Branches.templates import *
