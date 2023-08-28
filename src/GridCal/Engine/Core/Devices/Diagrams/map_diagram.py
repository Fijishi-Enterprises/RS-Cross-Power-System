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

import uuid
from typing import Dict, Union
from GridCal.Engine.Core.Devices.Diagrams.base_diagram import BaseDiagram
from GridCal.Engine.Core.Devices.Diagrams.map_location import MapLocation


class MapDiagram(BaseDiagram):
    """
    Diagram
    """

    def __init__(self, idtag: Union[None, str] = None, name: str = '') -> None:
        """
        
        :param idtag: UUID
        :param name: Diagram name
        """
        BaseDiagram.__init__(self, idtag=idtag, name=name)

