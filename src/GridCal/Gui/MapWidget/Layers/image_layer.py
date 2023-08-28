from typing import Callable, List, Union, Tuple, Dict
from GridCal.Gui.MapWidget.Layers.place import Place
from GridCal.Gui.MapWidget.Layers.layer_types import LayerType
from PySide6.QtGui import QPixmap


class ImageData:
    """
    PolylineDataPoint
    """

    def __init__(self,
                 lon: float,
                 lat: float,
                 pmap: QPixmap,
                 w: int,
                 h: int,
                 placement: Place,
                 offset_x: float,
                 offset_y: float,
                 radius: int,
                 colour: Tuple[int, int, int, int],
                 udata: Dict):
        """

        :param lon: longitude
        :param lat: latitude
        :param pmap: QPixmap
        :param w: width (px)
        :param h: height (px)
        :param placement: placement type
        :param offset_x: X offset in pixels
        :param offset_y: Y offset in pixels
        :param radius: Radius  (px)
        :param colour: colour tuple (R, G, B, A)
        :param udata: dictionary of extra properties
        """
        self.lon = lon
        self.lat = lat
        self.pmap = pmap
        self.w = w
        self.h = h
        self.placement = placement
        self.offset_x = offset_x
        self.offset_y = offset_y
        self.radius = radius
        self.colour = colour
        self.udata = udata


class ImageLayer:
    """
    A Layer object.
    """

    DefaultDelta = 50  # default selection delta

    def __init__(self,
                 layer_id: int = 0,
                 painter: Callable = None,
                 data: List[ImageData] = None,
                 map_rel: bool = True,
                 visible: bool = False,
                 show_levels: List[int] = None,
                 selectable: bool = False,
                 name: str = "<no name given>") -> None:
        """
        Initialise the Layer object.

        id           unique layer ID
        painter      render function
        data         the layer data
        map_rel      True if layer is map-relative, else layer-relative
        visible      layer visibility
        show_levels  list of levels at which to auto-show the level
        selectable   True if select operates on this layer, else False
        name         the name of the layer (for debug)
        ltype        a layer 'type' flag
        """

        self.painter = painter  # routine to draw layer
        self.data: List[ImageData] = data  # data that defines the layer
        self.map_rel = map_rel  # True if layer is map relative
        self.visible = visible  # True if layer visible
        self.show_levels = show_levels  # None or list of levels to auto-show
        self.selectable = selectable  # True if we can select on this layer
        self.delta = self.DefaultDelta  # minimum distance for selection
        self.name = name  # name of this layer
        self.type = LayerType.Image  # type of layer
        self.layer_id = layer_id  # ID of this layer

    @staticmethod
    def get_i18n_kw(kwargs, kws, default):
        """Get alternate international keyword value.

        kwargs   dictionary to look for keyword value
        kws      iterable of keyword spelling strings
        default  default value if no keyword found

        Returns the keyword value.
        """

        result = None
        for kw_str in kws[:-1]:
            result = kwargs.get(kw_str, None)
            if result:
                break
        else:
            result = kwargs.get(kws[-1], default)

        return result

    def setData(self, data: List[ImageData]):
        """

        :param data:
        :return:
        """

        self.data = data

    def __str__(self):
        return ('<Polyline Layer: id=%d, name=%s, map_rel=%s, visible=%s>'
                % (self.layer_id, self.name, str(self.map_rel), str(self.visible)))

