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
import numpy as np
from typing import List
from PySide6.QtCore import Qt, QPoint, QPointF
from PySide6.QtGui import QPen, QCursor, QColor
from PySide6.QtWidgets import QGraphicsItem, QGraphicsEllipseItem, QGraphicsRectItem

from GridCal.Engine.Core.Devices.Branches.transformer3w import Transformer3W
from GridCal.Gui.GridEditorWidget.generic_graphics import ACTIVE, DEACTIVATED
from GridCal.Gui.GuiFunctions import ObjectsModel
from GridCal.Gui.GridEditorWidget.terminal_item import TerminalItem
from GridCal.Gui.GridEditorWidget.winding_graphics import WindingGraphicItem
from GridCal.Gui.messages import yes_no_question


class Transformer3WGraphicItem(QGraphicsRectItem):
    """
      Represents a block in the diagram
      Has an x and y and width and height
      width and height can only be adjusted with a tip in the lower right corner.

      - in and output ports
      - parameters
      - description
    """

    def __init__(self, diagramScene: "DiagramScene",
                 editor: "GridEditor",
                 elm: Transformer3W,
                 pos: QPoint = None,
                 parent=None,
                 index=0):
        """

        :param diagramScene: DiagramScene object
        :param editor: GridEditor object
        :param elm: Transformer3W object
        :param pos: position
        :param parent:
        :param index:
        """
        QGraphicsRectItem.__init__(self, parent=parent)
        self.n_windings = 3
        self.min_w = 180.0
        self.min_h = 20.0
        self.offset = 10
        self.h = 70
        self.w = 80
        self.setRect(0.0, 0.0, self.w, self.h)

        self.api_object: Transformer3W = elm
        self.diagramScene = diagramScene  # this is the parent that hosts the pointer to the circuit
        self.editor = editor

        # color
        self.pen_width = 4
        if self.api_object is not None:
            if self.api_object.active:
                self.color = ACTIVE['color']
                self.style = ACTIVE['style']
            else:
                self.color = DEACTIVATED['color']
                self.style = DEACTIVATED['style']
        else:
            self.color = ACTIVE['color']
            self.style = ACTIVE['style']

        self.setPen(QPen(Qt.transparent, self.pen_width, self.style))
        self.setBrush(Qt.transparent)
        self.setFlag(QGraphicsItem.GraphicsItemFlag.ItemIsMovable)
        self.setFlag(QGraphicsItem.GraphicsItemFlag.ItemIsSelectable)
        self.setCursor(QCursor(Qt.PointingHandCursor))

        # index
        self.index = index

        if pos is not None:
            self.setPos(pos)

        # windings
        diameter = self.w * 0.5
        r = diameter / 2
        angle_0 = -90
        d_angle = 360 / self.n_windings
        angles_deg = [angle_0 + d_angle * i for i in range(self.n_windings)]
        angles = np.deg2rad(angles_deg)
        x = r * np.cos(angles) + self.w / 4
        y = r * np.sin(angles) + self.w / 4
        xt = diameter * np.cos(angles) + diameter
        yt = diameter * np.sin(angles) + diameter

        self.winding_circles: List[QGraphicsEllipseItem] = list()
        self.terminals: List[TerminalItem] = list()
        self.connection_lines: List[WindingGraphicItem | None] = list()

        pen = QPen(self.color, self.pen_width, self.style)

        for i in range(self.n_windings):
            # create objects
            winding_circle = QGraphicsEllipseItem(parent=self)
            winding_circle.setRect(0.0, 0.0, diameter, diameter)
            winding_circle.setPos(x[i], y[i])

            terminal = TerminalItem("t", parent=self, editor=self.editor)
            terminal.setPos(xt[i], yt[i])
            terminal.setRotation(angles_deg[i])

            # set objects style
            winding_circle.setPen(pen)
            terminal.setPen(pen)

            self.winding_circles.append(winding_circle)
            self.terminals.append(terminal)
            self.connection_lines.append(None)

        # set the graphical objects appropriately
        self.api_object.winding1.graphic_obj = self.winding_circles[0]

        self.big_marker = None

        # other actions
        self.set_winding_tool_tips()

    def recolour_mode(self):
        """
        Change the colour according to the system theme
        """
        if self.api_object is not None:
            if self.api_object.active:
                self.color = ACTIVE['color']
                self.style = ACTIVE['style']
            else:
                self.color = DEACTIVATED['color']
                self.style = DEACTIVATED['style']
        else:
            self.color = ACTIVE['color']
            self.style = ACTIVE['style']

        pen = QPen(self.color, self.pen_width, self.style)
        for i in range(self.n_windings):
            self.winding_circles[i].setPen(pen)
            self.terminals[i].setPen(pen)
            self.connection_lines[i].recolour_mode()

    def set_winding_tool_tips(self):
        """
        Set
        :return:
        """
        if self.api_object is not None:
            self.winding_circles[0].setToolTip("Winding 1: {0} kV".format(self.api_object.V1))
            self.winding_circles[1].setToolTip("Winding 2: {0} kV".format(self.api_object.V2))
            self.winding_circles[2].setToolTip("Winding 3: {0} kV".format(self.api_object.V3))
        pass

    def set_label(self, val: str):
        """
        Set the label content
        :param val:
        :return:
        """
        # this function is just for compatibility
        pass

    def mouseMoveEvent(self, event: 'QGraphicsSceneMouseEvent'):
        """
        On mouse move of this object...
        Args:
            event: QGraphicsSceneMouseEvent inherited
        """
        super().mouseMoveEvent(event)

        self.diagramScene.parent_.set_position(device=self.api_object,
                                               x=self.pos().x(),
                                               y=self.pos().y(),
                                               w=self.w,
                                               h=self.h,
                                               r=self.rotation())

        self.update_conn()

    def mousePressEvent(self, QGraphicsSceneMouseEvent):
        """
        mouse press: display the editor
        :param QGraphicsSceneMouseEvent:
        :return:
        """
        if self.api_object is not None:
            mdl = ObjectsModel([self.api_object], self.api_object.editable_headers,
                               parent=self.diagramScene.parent().object_editor_table,
                               editable=True, transposed=True,
                               non_editable_attributes=self.api_object.non_editable_attributes)

            self.diagramScene.parent().object_editor_table.setModel(mdl)

    def add_big_marker(self, color=Qt.red, tool_tip_text=""):
        """
        Add a big marker to the bus
        :param color: Qt Color ot the marker
        :param tool_tip_text: tool tip text to display
        :return:
        """
        if self.big_marker is None:
            self.big_marker = QGraphicsEllipseItem(0, 0, 180, 180, parent=self)
            self.big_marker.setBrush(color)
            self.big_marker.setOpacity(0.5)
            self.big_marker.setToolTip(tool_tip_text)

    def delete_big_marker(self):
        """
        Delete the big marker
        """
        if self.big_marker is not None:
            self.diagramScene.removeItem(self.big_marker)
            self.big_marker = None

    def change_size(self, w, h):
        """

        :param w:
        :param h:
        """
        # Keep for compatibility

        self.diagramScene.parent_.set_position(idtag=self.api_object.idtag,
                                               x=self.pos().x(),
                                               y=self.pos().y(),
                                               w=w,
                                               h=h,
                                               r=self.rotation())

    def set_position(self, x: float, y: float) -> None:
        """
        Set the bus x, y position
        :param x: x in pixels
        :param y: y in pixels
        """
        x = 0 if np.isnan(x) else int(x)
        y = 0 if np.isnan(y) else int(y)
        self.setPos(QPoint(int(x), int(y)))

    def update_conn(self) -> None:
        """
        Update the object
        """

        # Arrange line positions
        for terminal in self.terminals:
            x = int(2 * terminal.w / 3)
            y = int(2 * terminal.h / 3)
            m = QPointF(x, y)
            terminal.process_callbacks(self.pos() + terminal.pos() - m)

    def get_connection_winding(self, from_port: TerminalItem, to_port: TerminalItem):
        """
        Find the winding between the terminals
        :param from_port: "from" terminal [TerminalItem]
        :param to_port: "to" terminal [TerminalItem]
        """
        for i, t in enumerate(self.terminals):
            if t in [from_port, to_port]:
                return i

        raise Exception("Unknown winding")

    def set_connection(self, i: int, bus, conn: WindingGraphicItem):
        """
        Create the connection with a bus
        :param i: winding index 0-2
        :param bus: Bus object to connect to
        :param conn: Connection graphical object [LineGraphicItem]
        """
        if i == 0:
            self.api_object.bus1 = bus
            self.api_object.V1 = bus.Vnom
            self.connection_lines[0] = conn
            self.terminals[0].setZValue(-1)
            self.api_object.winding1.graphic_obj = conn
            conn.api_object = self.api_object.winding1

        elif i == 1:
            self.api_object.bus2 = bus
            self.api_object.V2 = bus.Vnom
            self.connection_lines[1] = conn
            self.terminals[1].setZValue(-1)
            self.api_object.winding2.graphic_obj = conn
            conn.api_object = self.api_object.winding2

        elif i == 2:
            self.api_object.bus3 = bus
            self.api_object.V3 = bus.Vnom
            self.connection_lines[2] = conn
            self.terminals[2].setZValue(-1)
            self.api_object.winding3.graphic_obj = conn
            conn.api_object = self.api_object.winding3

        # update the connection placement
        # from_port.update()
        # to_port.update()
        self.update_conn()
        self.mousePressEvent(None)

    def arrange_children(self) -> None:
        """
        this function is necessary because this graphic item behaves like a bus,
        but the the function itself does nothing
        """
        pass

    def set_tile_color(self, brush: QColor):
        """
        Set the voltage colour
        :param brush: QColor object
        """
        for w in self.winding_circles:
            w.setPen(QPen(brush, self.pen_width, self.style))

    def set_winding_color(self, i, color: QColor):
        """
        Set a winding (loading) colour
        :param i: winding index 0-2
        :param color: QColor
        """
        self.winding_circles[i].setPen(QPen(color, self.pen_width, self.style))
        self.terminals[i].setPen(QPen(color, self.pen_width, self.style))

    def delete_all_connections(self):
        """
        Delete all bus connections
        """
        for t in self.terminals:
            t.remove_all_connections()

        for c in self.connection_lines:
            self.diagramScene.removeItem(c)

    def remove(self, ask=True):
        """
        Remove this element
        @return:
        """
        if ask:
            ok = yes_no_question('Are you sure that you want to remove this bus',
                                 'Remove bus')
        else:
            ok = True

        if ok:
            self.delete_all_connections()
            self.diagramScene.removeItem(self)
            self.diagramScene.circuit.delete_transformer3w(self.api_object)
