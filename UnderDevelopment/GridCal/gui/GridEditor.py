from PyQt5.QtCore import *
from PyQt5.QtGui import *
from PyQt5.QtWidgets import *
from grid.CircuitOO import *
from gui.GuiFunctions import *


class LineUpdateMixin(object):
    def __init__(self, parent):
        super(LineUpdateMixin, self).__init__(parent)
        self.setFlag(QGraphicsItem.ItemSendsScenePositionChanges)

    def itemChange(self, change, value):
        if change == QGraphicsItem.ItemScenePositionHasChanged:
            self.parentItem().update_line(value)
        return super(LineUpdateMixin, self).itemChange(change, value)


class Polygon(LineUpdateMixin, QGraphicsPolygonItem):
    pass


class Square(LineUpdateMixin, QGraphicsRectItem):
    pass


class Circle(LineUpdateMixin, QGraphicsEllipseItem):
    pass


class QLine(LineUpdateMixin, QGraphicsLineItem):
    pass


class GeneralItem(object):

    def editParameters(self):
        pd = ParameterDialog(self.window())
        pd.exec_()

    def contextMenuEvent(self, event):
        menu = QMenu()
        pa = menu.addAction('Parameters')
        pa.triggered.connect(self.editParameters)

        ra1 = menu.addAction('Rotate +90')
        ra1.triggered.connect(self.rotate_clockwise)
        ra2 = menu.addAction('Rotate -90')
        ra2.triggered.connect(self.rotate_counterclockwise)

        ra3 = menu.addAction('Delete all the connections')
        ra3.triggered.connect(self.delete_all_connections)

        da = menu.addAction('Delete')
        da.triggered.connect(self.remove_)

        menu.exec_(event.screenPos())

    def rotate_clockwise(self):
        self.rotate(90)

    def rotate_counterclockwise(self):
        self.rotate(-90)

    def delete_all_connections(self):
        for term in self.terminals:
            term.remove_all_connections()
        # for term in self.lower_terminals:
        #     term.remove_all_connections()

    def remove_(self, delete_in_API=True):
        """

        @param delete_in_API:
        @return:
        """
        self.delete_all_connections()


class BranchGraphicItem(QGraphicsLineItem):
    """
    - fromPort
    - toPort
    """
    def __init__(self, fromPort, toPort, diagramScene, width=5, branch: Branch=None):
        """

        @param fromPort:
        @param toPort:
        @param diagramScene:
        """
        QGraphicsLineItem.__init__(self, None)

        self.api_object = branch
        if self.api_object is not None:
            if self.api_object.is_enabled:
                self.style = Qt.SolidLine
                self.color = Qt.black
            else:
                self.style = Qt.DashLine
                self.color = Qt.gray
        else:
            self.style = Qt.SolidLine
            self.color = Qt.black
        self.width = width
        self.pen_width = width
        self.setPen(QPen(Qt.black, self.width, self.style))
        self.setFlag(self.ItemIsSelectable, True)
        self.setCursor(QCursor(Qt.PointingHandCursor))

        self.pos1 = None
        self.pos2 = None
        self.fromPort = None
        self.toPort = None
        self.diagramScene = diagramScene

        if fromPort:
            self.setFromPort(fromPort)

        if toPort:
            self.setToPort(toPort)

        # Create arrow item:
        # self.line_object = LineItem(self)
        self.diagramScene.addItem(self)

        if fromPort and toPort:
            self.redraw()

    def contextMenuEvent(self, event):
        """
        Show context menu
        @param event:
        @return:
        """
        menu = QMenu()

        ra1 = menu.addAction('Properties')
        ra1.triggered.connect(self.editParameters)

        pe = menu.addAction('Enable/Disable')
        pe.triggered.connect(self.enable_disable_toggle)

        menu.addSeparator()

        ra2 = menu.addAction('Delete')
        ra2.triggered.connect(self.remove)

        menu.exec_(event.screenPos())

    def editParameters(self):
        """
        Display parameters editor for the Bus
        :return:
        """
        dialogue = QDialog(parent=self.diagramScene.parent())
        dialogue.setWindowTitle(self.api_object.name)
        layout = QVBoxLayout()
        grid = QTableView()
        layout.addWidget(grid)
        dialogue.setLayout(layout)

        mdl = ObjectsModel([self.api_object], self.api_object.edit_headers, self.api_object.edit_types,
                           parent=grid, editable=True, transposed=True, non_editable_indices=[1, 2])

        grid.setModel(mdl)
        dialogue.show()

    def mousePressEvent(self, QGraphicsSceneMouseEvent):
        """
        mouse press: display the editor
        :param QGraphicsSceneMouseEvent:
        :return:
        """
        mdl = ObjectsModel([self.api_object], self.api_object.edit_headers, self.api_object.edit_types,
                           parent=self.diagramScene.parent().object_editor_table, editable=True, transposed=True,
                           non_editable_indices=[1, 2])

        self.diagramScene.parent().object_editor_table.setModel(mdl)

    def remove(self):
        """
        Remove this object in the diagram and the API
        @return:
        """
        self.diagramScene.circuit.delete_branch(self.api_object)
        self.diagramScene.removeItem(self)

    def remove_(self):
        """
        Remove this object in the diagram
        @return:
        """
        self.diagramScene.removeItem(self)

    def enable_disable_toggle(self):
        """

        @return:
        """
        if self.api_object.is_enabled:
            self.set_enable(False)
        else:
            self.set_enable(True)

    def set_enable(self, val=True):
        """
        Set the enable value, graphically and in the API
        @param val:
        @return:
        """
        self.api_object.is_enabled = val
        if self.api_object is not None:
            if self.api_object.is_enabled:
                self.style = Qt.SolidLine
                self.color = QtCore.Qt.black
            else:
                self.style = Qt.DashLine
                self.color = QtCore.Qt.gray
        else:
            self.style = Qt.SolidLine
            self.color = QtCore.Qt.black
        self.setPen(QPen(self.color, self.width, self.style))

    def setFromPort(self, fromPort):
        """
        Set the From terminal in a connection
        @param fromPort:
        @return:
        """
        self.fromPort = fromPort
        if self.fromPort:
            self.pos1 = fromPort.scenePos()
            self.fromPort.posCallbacks.append(self.setBeginPos)
            self.fromPort.setZValue(0)

    def setToPort(self, toPort):
        """
        Set the To terminal in a connection
        @param toPort:
        @return:
        """
        self.toPort = toPort
        if self.toPort:
            self.pos2 = toPort.scenePos()
            self.toPort.posCallbacks.append(self.setEndPos)
            self.toPort.setZValue(0)

    def setEndPos(self, endpos):
        """
        Set the starting position
        @param endpos:
        @return:
        """
        self.pos2 = endpos
        self.redraw()

    def setBeginPos(self, pos1):
        """
        Set the starting position
        @param pos1:
        @return:
        """
        self.pos1 = pos1
        self.redraw()

    def redraw(self):
        """
        Redraw the line with the given positions
        @return:
        """
        self.setLine(QLineF(self.pos1, self.pos2))
        self.setZValue(0)


class ParameterDialog(QDialog):

    def __init__(self, parent=None):
        super(ParameterDialog, self).__init__(parent)
        self.button = QPushButton('Ok', self)
        l = QVBoxLayout(self)
        l.addWidget(self.button)
        self.button.clicked.connect(self.OK)

    def OK(self):
        self.close()


class TerminalItem(QGraphicsEllipseItem):
    """
    Represents a connection point to a subsystem
    """

    def __init__(self, name, editor=None, parent=None, h=10, w=10):
        """

        @param name:
        @param editor:
        @param parent:
        """
        QGraphicsEllipseItem.__init__(self, QRectF(-6, -6, h, w), parent)
        self.setCursor(QCursor(QtCore.Qt.CrossCursor))

        # Properties:
        self.setBrush(QBrush(Qt.white))

        # terminal parent object
        self.parent = parent

        self.hosting_connections = list()

        self.editor = editor

        # Name:
        self.name = name
        self.posCallbacks = []
        self.setFlag(self.ItemSendsScenePositionChanges, True)

    def itemChange(self, change, value):
        """

        @param change:
        @param value:
        @return:
        """
        if change == self.ItemScenePositionHasChanged:
            for cb in self.posCallbacks:
                cb(value)
            return value
        return super(TerminalItem, self).itemChange(change, value)

    def mousePressEvent(self, event):
        """
        Start a connection
        Args:
            event:

        Returns:

        """
        self.editor.startConnection(self)
        self.hosting_connections.append(self.editor.startedConnection)

    def remove_all_connections(self):
        """
        Removes all the terminal connections
        Returns:

        """
        n = len(self.hosting_connections)
        for i in range(n-1, -1, -1):
            self.hosting_connections[i].remove_()
            self.hosting_connections.pop(i)


class HandleItem(QGraphicsEllipseItem):
    """
    A handle that can be moved by the mouse: Element to resize the boxes
    """
    def __init__(self, parent=None):
        """

        @param parent:
        """
        QGraphicsEllipseItem.__init__(self, QRectF(-4, -4, 8, 8), parent)
        # super(HandleItem, self).__init__(QRectF(-4, -4, 8, 8), parent)
        self.posChangeCallbacks = []
        self.setBrush(QBrush(Qt.red))
        self.setFlag(self.ItemIsMovable, True)
        self.setFlag(self.ItemSendsScenePositionChanges, True)
        self.setCursor(QCursor(Qt.SizeFDiagCursor))

    def itemChange(self, change, value):
        """

        @param change:
        @param value:
        @return:
        """
        if change == self.ItemPositionChange:
            x, y = value.x(), value.y()
            # TODO: make this a signal?
            # This cannot be a signal because this is not a QObject
            for cb in self.posChangeCallbacks:
                res = cb(x, y)
                if res:
                    x, y = res
                    value = QPointF(x, y)
            return value

        # Call superclass method:
        return super(HandleItem, self).itemChange(change, value)


class LoadGraphicItem(QGraphicsItemGroup):

    def __init__(self, parent, api_obj, diagramScene):
        """

        :param parent:
        :param api_obj:
        """
        # QGraphicsPolygonItem.__init__(self, parent=parent)
        # QGraphicsItemGroup.__init__(self, parent=parent)
        super(LoadGraphicItem, self).__init__(parent)

        self.w = 20.0
        self.h = 20.0

        self.parent = parent

        self.api_object = api_obj

        self.diagramScene = diagramScene

        # Properties of the container:
        # self.setBrush(QtGui.QBrush(QtCore.Qt.black))
        self.setFlags(self.ItemIsSelectable | self.ItemIsMovable)
        self.setCursor(QCursor(QtCore.Qt.PointingHandCursor))
        # self.installSceneEventFilter(self)

        triangle = Polygon(self)
        triangle.setPolygon(QPolygonF([QPointF(0, 0), QPointF(self.w, 0), QPointF(self.w/2, self.h)]))
        triangle.setPen(QPen(Qt.red, 2))

        # line = QGraphicsLineItem(QLineF(QPointF(self.w/2, -10), QPointF(self.w/2, 0)))
        # line .setPen(QPen(Qt.red, 2))
        # triangle.setPos(10, 30)
        # self.
        self.addToGroup(triangle)
        # self.addToGroup(line)

        # line to tie this object with the original bus (the parent)
        self.nexus = QGraphicsLineItem()
        parent.scene().addItem(self.nexus)
        self.update_line(self.pos())

    def update_line(self, pos):
        parent = self.parentItem()
        rect = parent.rect()
        self.nexus.setLine(
            pos.x() + self.w/2, pos.y() + 0,
            parent.x() + rect.width() / 2,
            parent.y() + rect.height(),
        )

    def contextMenuEvent(self, event):
        """
        Display context menu
        @param event:
        @return:
        """
        menu = QMenu()

        da = menu.addAction('Delete')
        da.triggered.connect(self.remove)

        menu.exec_(event.screenPos())

    def remove(self):
        """
        Remove this element
        @return:
        """
        self.diagramScene.removeItem(self.nexus)
        self.diagramScene.removeItem(self)
        self.api_object.bus.loads.remove(self.api_object)

    def mousePressEvent(self, QGraphicsSceneMouseEvent):
        """
        mouse press: display the editor
        :param QGraphicsSceneMouseEvent:
        :return:
        """
        mdl = ObjectsModel([self.api_object], self.api_object.edit_headers, self.api_object.edit_types,
                           parent=self.diagramScene.parent().object_editor_table, editable=True, transposed=True)
        self.diagramScene.parent().object_editor_table.setModel(mdl)


class ShuntGraphicItem(QGraphicsItemGroup):

    def __init__(self, parent, api_obj, diagramScene):
        """

        :param parent:
        :param api_obj:
        """
        # QGraphicsPolygonItem.__init__(self, parent=parent)
        # QGraphicsItemGroup.__init__(self, parent=parent)
        super(ShuntGraphicItem, self).__init__(parent)

        self.w = 10.0
        self.h = 30.0

        self.parent = parent

        self.api_object = api_obj

        self.diagramScene = diagramScene

        pen = QPen(Qt.red, 2)

        # Properties of the container:
        # self.setBrush(QtGui.QBrush(QtCore.Qt.black))
        self.setFlags(self.ItemIsSelectable | self.ItemIsMovable)
        self.setCursor(QCursor(QtCore.Qt.PointingHandCursor))
        # self.installSceneEventFilter(self)

        # line to tie this object with the original bus (the parent)
        self.nexus = QGraphicsLineItem()
        parent.scene().addItem(self.nexus)

        lines = list()
        lines.append(QLineF(QPointF(self.w/2, 0), QPointF(self.w/2, self.h*0.4)))
        lines.append(QLineF(QPointF(0, self.h*0.4), QPointF(self.w, self.h*0.4)))
        lines.append(QLineF(QPointF(0, self.h*0.6), QPointF(self.w, self.h*0.6)))
        lines.append(QLineF(QPointF(self.w/2, self.h*0.6), QPointF(self.w/2, self.h)))
        for l in lines:
            l1 = QLine(self)
            l1.setLine(l)
            l1.setPen(pen)
            self.addToGroup(l1)

        # line = QLine(self)
        # line.setLine(QLineF(QPointF(self.w/2, 0), QPointF(self.w/2, self.h*0.4)))
        # line.setPen(pen)
        # self.addToGroup(line)

        self.setPos(self.parent.x(), self.parent.y() + 100)
        self.update_line(self.pos())

    def update_line(self, pos):
        parent = self.parentItem()
        rect = parent.rect()
        self.nexus.setLine(
            pos.x() + self.w/2,
            pos.y() + 0,
            parent.x() + rect.width() / 2,
            parent.y() + rect.height(),
        )

    def contextMenuEvent(self, event):
        """
        Display context menu
        @param event:
        @return:
        """
        menu = QMenu()

        da = menu.addAction('Delete')
        da.triggered.connect(self.remove)

        menu.exec_(event.screenPos())

    def remove(self):
        """
        Remove this element
        @return:
        """
        self.diagramScene.removeItem(self.nexus)
        self.diagramScene.removeItem(self)
        self.api_object.bus.shunts.remove(self.api_object)

    def mousePressEvent(self, QGraphicsSceneMouseEvent):
        """
        mouse press: display the editor
        :param QGraphicsSceneMouseEvent:
        :return:
        """
        mdl = ObjectsModel([self.api_object], self.api_object.edit_headers, self.api_object.edit_types,
                           parent=self.diagramScene.parent().object_editor_table, editable=True, transposed=True)
        self.diagramScene.parent().object_editor_table.setModel(mdl)


class ControlledGeneratorGraphicItem(QGraphicsItemGroup):

    def __init__(self, parent, api_obj, diagramScene):
        """

        :param parent:
        :param api_obj:
        """
        # QGraphicsPolygonItem.__init__(self, parent=parent)
        # QGraphicsItemGroup.__init__(self, parent=parent)

        super(ControlledGeneratorGraphicItem, self).__init__(parent)

        # self.w = 60.0
        # self.h = 60.0

        self.parent = parent

        self.api_object = api_obj

        self.diagramScene = diagramScene

        color = Qt.red
        pen = QPen(color, 2)

        self.w = 40
        self.h = 40

        # Properties of the container:
        # self.setBrush(QtGui.QBrush(QtCore.Qt.black))
        self.setFlags(self.ItemIsSelectable | self.ItemIsMovable)
        self.setCursor(QCursor(QtCore.Qt.PointingHandCursor))
        # self.installSceneEventFilter(self)

        # line to tie this object with the original bus (the parent)
        self.nexus = QGraphicsLineItem()
        parent.scene().addItem(self.nexus)

        # l1 = QGraphicsLineItem(QLineF(QPointF(self.w/2, 0), QPointF(self.w/2, -10)))
        # l1.setPen(pen)
        # self.addToGroup(l1)

        circle = Circle(parent)
        circle.setRect(0, 0, self.h, self.w)
        circle.setPen(pen)
        self.addToGroup(circle)

        label = QGraphicsTextItem('G', parent=circle)
        label.setDefaultTextColor(color)
        label.setPos(self.h/4, self.w/4)

        self.setPos(self.parent.x(), self.parent.y() + 100)
        self.update_line(self.pos())

    def update_line(self, pos):
        parent = self.parentItem()
        rect = parent.rect()
        self.nexus.setLine(
            pos.x() + self.w/2, pos.y() + 0,
            parent.x() + rect.width() / 2,
            parent.y() + rect.height(),
        )

    def contextMenuEvent(self, event):
        """
        Display context menu
        @param event:
        @return:
        """
        menu = QMenu()

        da = menu.addAction('Delete')
        da.triggered.connect(self.remove)

        menu.exec_(event.screenPos())

    def remove(self):
        """
        Remove this element
        @return:
        """
        self.diagramScene.removeItem(self.nexus)
        self.diagramScene.removeItem(self)
        self.api_object.bus.controlled_generators.remove(self.api_object)

    def mousePressEvent(self, QGraphicsSceneMouseEvent):
        """
        mouse press: display the editor
        :param QGraphicsSceneMouseEvent:
        :return:
        """
        mdl = ObjectsModel([self.api_object], self.api_object.edit_headers, self.api_object.edit_types,
                           parent=self.diagramScene.parent().object_editor_table, editable=True, transposed=True)
        self.diagramScene.parent().object_editor_table.setModel(mdl)


class StaticGeneratorGraphicItem(QGraphicsItemGroup):

    def __init__(self, parent, api_obj, diagramScene):
        """

        :param parent:
        :param api_obj:
        """
        # QGraphicsPolygonItem.__init__(self, parent=parent)
        # QGraphicsItemGroup.__init__(self, parent=parent)
        super(StaticGeneratorGraphicItem, self).__init__(parent)

        self.parent = parent

        self.api_object = api_obj

        self.diagramScene = diagramScene

        color = Qt.red
        pen = QPen(color, 2)

        self.w = 40
        self.h = 40

        # Properties of the container:
        # self.setBrush(QtGui.QBrush(QtCore.Qt.black))
        self.setFlags(self.ItemIsSelectable | self.ItemIsMovable)
        self.setCursor(QCursor(QtCore.Qt.PointingHandCursor))

        # line to tie this object with the original bus (the parent)
        self.nexus = QGraphicsLineItem()
        parent.scene().addItem(self.nexus)

        # l1 = QGraphicsLineItem(QLineF(QPointF(self.w/2, 0), QPointF(self.w/2, -10)))
        # l1.setPen(pen)
        # self.addToGroup(l1)

        square = Square(parent)
        square.setRect(0, 0, self.h, self.w)
        square.setPen(pen)
        self.addToGroup(square)

        label = QGraphicsTextItem('S', parent=square)
        label.setDefaultTextColor(color)
        label.setPos(self.h/4, self.w/4)

        self.setPos(self.parent.x(), self.parent.y() + 100)
        self.update_line(self.pos())

    def update_line(self, pos):
        parent = self.parentItem()
        rect = parent.rect()
        self.nexus.setLine(
            pos.x() + self.w/2, pos.y() + 0,
            parent.x() + rect.width() / 2,
            parent.y() + rect.height(),
        )

    def contextMenuEvent(self, event):
        """
        Display context menu
        @param event:
        @return:
        """
        menu = QMenu()

        da = menu.addAction('Delete')
        da.triggered.connect(self.remove)

        menu.exec_(event.screenPos())

    def remove(self):
        """
        Remove this element
        @return:
        """
        self.diagramScene.removeItem(self.nexus)
        self.diagramScene.removeItem(self)
        self.api_object.bus.static_generators.remove(self.api_object)

    def mousePressEvent(self, QGraphicsSceneMouseEvent):
        """
        mouse press: display the editor
        :param QGraphicsSceneMouseEvent:
        :return:
        """
        mdl = ObjectsModel([self.api_object], self.api_object.edit_headers, self.api_object.edit_types,
                           parent=self.diagramScene.parent().object_editor_table, editable=True, transposed=True)
        self.diagramScene.parent().object_editor_table.setModel(mdl)


class BatteryGraphicItem(QGraphicsItemGroup):

    def __init__(self, parent, api_obj, diagramScene):
        """

        :param parent:
        :param api_obj:
        """
        # QGraphicsPolygonItem.__init__(self, parent=parent)
        # QGraphicsItemGroup.__init__(self, parent=parent)
        super(BatteryGraphicItem, self).__init__(parent)

        self.parent = parent

        self.api_object = api_obj

        self.diagramScene = diagramScene

        color = Qt.red
        pen = QPen(color, 2)

        self.w = 40
        self.h = 40

        # Properties of the container:
        # self.setBrush(QtGui.QBrush(QtCore.Qt.black))
        self.setFlags(self.ItemIsSelectable | self.ItemIsMovable)
        self.setCursor(QCursor(QtCore.Qt.PointingHandCursor))

        # line to tie this object with the original bus (the parent)
        self.nexus = QGraphicsLineItem()
        parent.scene().addItem(self.nexus)

        # l1 = QGraphicsLineItem(QLineF(QPointF(self.w/2, 0), QPointF(self.w/2, -10)))
        # l1.setPen(pen)
        # self.addToGroup(l1)

        square = Square(self)
        square.setRect(0, 0, self.h, self.w)
        square.setPen(pen)
        self.addToGroup(square)

        label = QGraphicsTextItem('B', parent=square)
        label.setDefaultTextColor(color)
        label.setPos(self.h/4, self.w/4)

        self.setPos(self.parent.x(), self.parent.y() + 100)
        self.update_line(self.pos())

    def update_line(self, pos):
        parent = self.parentItem()
        rect = parent.rect()
        self.nexus.setLine(
            pos.x() + self.w/2, pos.y() + 0,
            parent.x() + rect.width() / 2,
            parent.y() + rect.height(),
        )

    def contextMenuEvent(self, event):
        """
        Display context menu
        @param event:
        @return:
        """
        menu = QMenu()

        da = menu.addAction('Delete')
        da.triggered.connect(self.remove)

        menu.exec_(event.screenPos())

    def remove(self):
        """
        Remove this element
        @return:
        """
        self.diagramScene.removeItem(self.nexus)
        self.diagramScene.removeItem(self)
        self.api_object.bus.batteries.remove(self.api_object)

    def mousePressEvent(self, QGraphicsSceneMouseEvent):
        """
        mouse press: display the editor
        :param QGraphicsSceneMouseEvent:
        :return:
        """
        mdl = ObjectsModel([self.api_object], self.api_object.edit_headers, self.api_object.edit_types,
                           parent=self.diagramScene.parent().object_editor_table, editable=True, transposed=True)
        self.diagramScene.parent().object_editor_table.setModel(mdl)


class BusGraphicItem(QGraphicsRectItem, GeneralItem):
    """
      Represents a block in the diagram
      Has an x and y and width and height
      width and height can only be adjusted with a tip in the lower right corner.

      - in and output ports
      - parameters
      - description
    """
    def __init__(self, diagramScene, name='Untitled', parent=None, index=0, editor=None,
                 bus: Bus=None, pos: QPoint=None):
        """

        @param diagramScene:
        @param name:
        @param parent:
        @param index:
        @param editor:
        """
        # QGraphicsRectItem.__init__(self, parent=parent)
        # GeneralItem.__init__(self)
        super(BusGraphicItem, self).__init__(parent)

        self.w = 60.0
        self.h = 60.0

        self.api_object = bus

        self.diagramScene = diagramScene

        self.editor = editor

        self.graphic_children = list()

        # Properties of the rectangle:
        self.setPen(QPen(QtCore.Qt.black, 2))
        self.setBrush(QBrush(QtCore.Qt.black))
        self.setFlags(self.ItemIsSelectable | self.ItemIsMovable)
        self.setCursor(QCursor(QtCore.Qt.PointingHandCursor))

        # index
        self.index = index

        if pos is not None:
            self.setPos(pos)

        # Label:
        self.label = QGraphicsTextItem(bus.name, self)
        self.label.setDefaultTextColor(QtCore.Qt.white)

        # Create corner for resize:
        self.sizer = HandleItem(self)
        self.sizer.setPos(self.w, self.h)
        self.sizer.posChangeCallbacks.append(self.change_size)  # Connect the callback

        self.sizer.setFlag(self.sizer.ItemIsSelectable, True)

        # connection terminals the block:
        self.upper_terminals = []
        self.upper_terminals.append(TerminalItem('n', parent=self, editor=self.editor))  # , h=self.h))
        self.lower_terminals = []
        self.lower_terminals.append(TerminalItem('s', parent=self, editor=self.editor))  # , h=self.h))
        self.right_terminals = []
        self.right_terminals.append(TerminalItem('e', parent=self, editor=self.editor))  # , w=self.w))
        self.left_terminals = []
        self.left_terminals.append(TerminalItem('w', parent=self, editor=self.editor))  # , w=self.w))

        self.terminals = self.upper_terminals + self.lower_terminals + self.right_terminals + self.left_terminals

        # Update size:
        self.change_size(self.w, self.h)

    def change_size(self, w, h):
        """
        Resize block function
        @param w:
        @param h:
        @return:
        """
        # Limit the block size to the minimum size:
        if h < self.h:
            h = self.h
        if w < self.w:
            w = self.w
        self.setRect(0.0, 0.0, w, h)

        offset = 10

        # center label:
        rect = self.label.boundingRect()
        lw, lh = rect.width(), rect.height()
        lx = (w - lw) / 2
        ly = (h - lh) / 2
        self.label.setPos(lx, ly)

        # upper
        n = len(self.upper_terminals)
        y0 = -offset/2
        dx = w / (n+1)
        x0 = dx
        for term in self.upper_terminals:
            term.setPos(x0, y0)
            # term.setPos(x0 - w / 2 + offset / 2, y0)
            x0 += dx

        # lower
        n = len(self.lower_terminals)
        y0 = h + offset
        dx = w / (n+1)
        x0 = dx
        for term in self.lower_terminals:
            term.setPos(x0, y0)
            # term.setPos(x0 - w / 2 + offset / 2, y0)
            x0 += dx

        # right
        n = len(self.right_terminals)
        x0 = w + offset
        dy = h / (n+1)
        y0 = dy
        for term in self.right_terminals:
            term.setPos(x0, y0)
            # term.setPos(x0, y0 - h / 2 + offset / 2)
            y0 += dy

        # left
        n = len(self.left_terminals)
        x0 = - offset
        dy = h / (n+1)
        y0 = dy
        for term in self.left_terminals:
            term.setPos(x0, y0)
            # term.setPos(x0, y0 - h / 2 + offset / 2)
            y0 += dy

        return w, h

    def arrange_children(self):
        """
        This function sorts the load and generators icons
        Returns:
            Nothing
        """
        y0 = self.h + 40
        x = 0
        # print(x, y0)
        for elm in self.graphic_children:
            elm.setPos(x, y0)
            x += elm.w + 10

    def create_children_icons(self):
        """
        Create the icons of the elements that are attached to the API bus object
        Returns:
            Nothing
        """
        for elm in self.api_object.loads:
            self.add_load(elm)

        for elm in self.api_object.static_generators:
            self.add_static_generator(elm)

        for elm in self.api_object.controlled_generators:
            self.add_controlled_generator(elm)

        for elm in self.api_object.shunts:
            self.add_shunt(elm)

        for elm in self.api_object.batteries:
            self.add_battery(elm)

        self.arrange_children()

    def contextMenuEvent(self, event):
        """
        Display context menu
        @param event:
        @return:
        """
        menu = QMenu()
        # pa = menu.addAction('Parameters')
        # pa.triggered.connect(self.editParameters)

        pe = menu.addAction('Enable/Disable')
        pe.triggered.connect(self.enable_disable_toggle)

        pl = menu.addAction('Plot profiles')
        pl.triggered.connect(self.plot_profiles)

        ra1 = menu.addAction('Rotate +90')
        ra1.triggered.connect(self.rotate_clockwise)
        ra2 = menu.addAction('Rotate -90')
        ra2.triggered.connect(self.rotate_counterclockwise)

        menu.addSeparator()

        ra3 = menu.addAction('Delete all the connections')
        ra3.triggered.connect(self.delete_all_connections)

        da = menu.addAction('Delete')
        da.triggered.connect(self.remove)

        menu.addSeparator()

        al = menu.addAction('Add load')
        al.triggered.connect(self.add_load)

        ash = menu.addAction('Add shunt')
        ash.triggered.connect(self.add_shunt)

        acg = menu.addAction('Add controlled generator')
        acg.triggered.connect(self.add_controlled_generator)

        asg = menu.addAction('Add static generator')
        asg.triggered.connect(self.add_static_generator)

        ab = menu.addAction('Add battery')
        ab.triggered.connect(self.add_battery)

        menu.addSeparator()

        arr = menu.addAction('Arrange')
        arr.triggered.connect(self.arrange_children)

        menu.exec_(event.screenPos())

    def remove(self):
        """
        Remove this element
        @return:
        """
        self.delete_all_connections()
        self.diagramScene.removeItem(self)
        self.diagramScene.circuit.delete_bus(self.api_object)

    def enable_disable_toggle(self):
        """
        Toggle bus element state
        @return:
        """
        self.api_object.is_enabled = not self.api_object.is_enabled
        print('Enabled:', self.api_object.is_enabled)

        if self.api_object.is_enabled:
            self.setBrush(QBrush(QtCore.Qt.black))

            for term in self.terminals:
                for host in term.hosting_connections:
                    host.set_enable(val=True)
        else:
            self.setBrush(QBrush(QtCore.Qt.gray))

            for term in self.terminals:
                for host in term.hosting_connections:
                    host.set_enable(val=False)

    def plot_profiles(self):
        """

        @return:
        """
        # t = self.diagramScene.circuit.master_time_array
        # self.api_object.plot_profiles(time_idx=t)
        self.api_object.plot_profiles()

    def editParameters(self):
        """
        Display parameters editor for the Bus
        :return:
        """
        dialogue = QDialog(parent=self.diagramScene.parent())
        dialogue.setWindowTitle(self.api_object.name)
        layout = QVBoxLayout()
        grid = QTableView()
        layout.addWidget(grid)
        dialogue.setLayout(layout)

        mdl = ObjectsModel([self.api_object], self.api_object.edit_headers, self.api_object.edit_types,
                           parent=grid, editable=True, transposed=True)

        grid.setModel(mdl)
        dialogue.show()

    def mousePressEvent(self, QGraphicsSceneMouseEvent):
        """
        mouse press: display the editor
        :param QGraphicsSceneMouseEvent:
        :return:
        """
        mdl = ObjectsModel([self.api_object], self.api_object.edit_headers, self.api_object.edit_types,
                           parent=self.diagramScene.parent().object_editor_table, editable=True, transposed=True)
        self.diagramScene.parent().object_editor_table.setModel(mdl)

    def add_load(self, api_obj=None):
        """

        Returns:

        """
        if api_obj is None or type(api_obj) is bool:
            api_obj = Load()
            api_obj.bus = self.api_object
            self.api_object.loads.append(api_obj)

        _grph = LoadGraphicItem(self, api_obj, self.diagramScene)
        api_obj.graphic_obj = _grph
        self.graphic_children.append(_grph)
        self.arrange_children()

    def add_shunt(self, api_obj=None):
        """

        Returns:

        """
        if api_obj is None or type(api_obj) is bool:
            api_obj = Shunt()
            api_obj.bus = self.api_object
            self.api_object.shunts.append(api_obj)

        _grph = ShuntGraphicItem(self, api_obj, self.diagramScene)
        api_obj.graphic_obj = _grph
        self.graphic_children.append(_grph)
        self.arrange_children()

    def add_controlled_generator(self, api_obj=None):
        """

        Returns:

        """
        if api_obj is None or type(api_obj) is bool:
            api_obj = ControlledGenerator()
            api_obj.bus = self.api_object
            self.api_object.controlled_generators.append(api_obj)

        _grph = ControlledGeneratorGraphicItem(self, api_obj, self.diagramScene)
        api_obj.graphic_obj = _grph
        self.graphic_children.append(_grph)
        self.arrange_children()

    def add_static_generator(self, api_obj=None):
        """

        Returns:

        """
        if api_obj is None or type(api_obj) is bool:
            api_obj = StaticGenerator()
            api_obj.bus = self.api_object
            self.api_object.static_generators.append(api_obj)

        _grph = StaticGeneratorGraphicItem(self, api_obj, self.diagramScene)
        api_obj.graphic_obj = _grph
        self.graphic_children.append(_grph)
        self.arrange_children()

    def add_battery(self, api_obj=None):
        """

        Returns:

        """
        if api_obj is None or type(api_obj) is bool:
            api_obj = Battery()
            api_obj.bus = self.api_object
            self.api_object.batteries.append(api_obj)

        _grph = BatteryGraphicItem(self, api_obj, self.diagramScene)
        api_obj.graphic_obj = _grph
        self.graphic_children.append(_grph)
        self.arrange_children()


class EditorGraphicsView(QGraphicsView):
    """
    Editor where the diagram is displayed
    """
    def __init__(self, scene, parent=None, editor=None):
        """

        @param scene:
        @param parent:
        @param editor:
        """
        QGraphicsView.__init__(self, scene, parent)

        # self.setBackgroundBrush(QColor(0,66,255,180))
        self.setDragMode(QGraphicsView.RubberBandDrag)
        self.setRubberBandSelectionMode(Qt.IntersectsItemShape)
        self.setMouseTracking(True)
        self.setInteractive(True)
        self.scene_ = scene
        self.setRenderHints(QPainter.Antialiasing | QPainter.SmoothPixmapTransform)
        self.editor = editor
        self.last_n = 1
        self.setAlignment(Qt.AlignCenter)

    def dragEnterEvent(self, event):
        """

        @param event:
        @return:
        """
        if event.mimeData().hasFormat('component/name'):
            event.accept()

    def dragMoveEvent(self, event):
        """
        Move element
        @param event:
        @return:
        """
        if event.mimeData().hasFormat('component/name'):
            event.accept()

    def dropEvent(self, event):
        """
        Create an element
        @param event:
        @return:
        """
        if event.mimeData().hasFormat('component/name'):
            objtype = event.mimeData().data('component/name')
            # name = str(objtype)

            print(str(event.mimeData().data('component/name')))

            elm = None
            data = QByteArray()
            stream = QDataStream(data, QIODevice.WriteOnly)
            stream.writeQString('Bus')
            if objtype == data:
                name = 'Bus ' + str(self.last_n)
                self.last_n += 1
                obj = Bus(name=name)
                elm = BusGraphicItem(diagramScene=self.scene(), name=name, editor=self.editor, bus=obj)
                obj.graphic_obj = elm
                self.scene_.circuit.add_bus(obj)  # weird but only way to have graphical-API communication

            if elm is not None:
                elm.setPos(self.mapToScene(event.pos()))
                self.scene_.addItem(elm)
                # self.scene_.circuit.add_bus(obj) # weird but only way to have graphical-API communication
                print('Block created')

    def wheelEvent(self, event):
        """
        Zoom
        @param event:
        @return:
        """
        self.setTransformationAnchor(QGraphicsView.AnchorUnderMouse)

        # Scale the view / do the zoom
        scale_factor = 1.15
        # print(event.angleDelta().x(), event.angleDelta().y(), event.angleDelta().manhattanLength() )
        if event.angleDelta().y() > 0:
            # Zoom in
            self.scale(scale_factor, scale_factor)

        else:
            # Zooming out
            self.scale(1.0 / scale_factor, 1.0 / scale_factor)

    def add_bus(self, bus: Bus, explode_factor=1.0):
        elm = BusGraphicItem(diagramScene=self.scene(), name=bus.name, editor=self.editor, bus=bus)
        elm.setPos(self.mapToScene(QPoint(bus.x * explode_factor, bus.y * explode_factor)))
        self.scene_.addItem(elm)
        return elm


class LibraryModel(QStandardItemModel):
    """
    Items model to host the draggable icons
    """
    def __init__(self, parent=None):
        """

        @param parent:
        """
        QStandardItemModel.__init__(self, parent)

    def mimeTypes(self):
        """

        @return:
        """
        return ['component/name']

    def mimeData(self, idxs):
        """

        @param idxs:
        @return:
        """
        mimedata = QMimeData()
        for idx in idxs:
            if idx.isValid():
                txt = self.data(idx, Qt.DisplayRole)

                data = QByteArray()
                stream = QDataStream(data, QIODevice.WriteOnly)
                stream.writeQString(txt)

                mimedata.setData('component/name', data)
        return mimedata


class DiagramScene(QGraphicsScene):

    def __init__(self, parent=None, circuit: MultiCircuit=None):
        """

        @param parent:
        """
        super(DiagramScene, self).__init__(parent)
        self.parent_ = parent
        self.circuit = circuit

    def mouseMoveEvent(self, mouseEvent):
        """

        @param mouseEvent:
        @return:
        """
        self.parent_.sceneMouseMoveEvent(mouseEvent)
        super(DiagramScene, self).mouseMoveEvent(mouseEvent)

    def mouseReleaseEvent(self, mouseEvent):
        """

        @param mouseEvent:
        @return:
        """
        self.parent_.sceneMouseReleaseEvent(mouseEvent)
        super(DiagramScene, self).mouseReleaseEvent(mouseEvent)


class ObjectFactory(object):

    def get_box(self):
        """

        @return:
        """
        pixmap = QPixmap(40, 40)
        pixmap.fill()
        painter = QPainter(pixmap)
        painter.fillRect(0, 0, 40, 40, Qt.black)
        # painter.setBrush(Qt.red)
        # painter.drawEllipse(36, 2, 20, 20)
        # painter.setBrush(Qt.yellow)
        # painter.drawEllipse(20, 20, 20, 20)
        painter.end()

        return QIcon(pixmap)

    def get_circle(self):
        """

        @return:
        """
        pixmap = QPixmap(40, 40)
        pixmap.fill()
        painter = QPainter(pixmap)
        # painter.fillRect(10, 10, 80, 80, Qt.black)
        painter.setBrush(Qt.red)
        painter.drawEllipse(0, 0, 40, 40)
        # painter.setBrush(Qt.yellow)
        # painter.drawEllipse(20, 20, 20, 20)
        painter.end()

        return QIcon(pixmap)