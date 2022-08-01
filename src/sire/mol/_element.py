
from ..legacy.Mol import Element


def _color(obj):
    """Return the color of the element, as a RGB triple
       (three integers from 0 to 255). Note that this is
       the color typically used to represent the element
       in a molecular viewer, not the actual color of
       the element.
    """
    def _float_to_int(v):
        i = int(v * 255)
        return max(0, min(i, 255))

    return (_float_to_int(obj.red()),
            _float_to_int(obj.green()),
            _float_to_int(obj.blue()))


def _hex_color(obj):
    """Return the color of the element as a hex string"""
    rgb = obj.color()

    def to_hex(v):
        h = hex(v)[2:]
        if len(h) < 2:
            return "0"+h
        else:
            return h

    return f"0x{to_hex(rgb[0])}{to_hex(rgb[1])}{to_hex(rgb[2])}"


def _color_name(obj):
    """Return the color name of the element. Note that this is
       the color typically used to represent the element
       in a molecular viewer, not the actual color of
       the element.
    """
    colors = {
        1: "white",     # hydrogen
        6: "charcoal",  # carbon
        7: "blue",      # nitrogen
        8: "red",       # oxygen
       16: "yellow"     # sulphur
    }

    try:
        return colors[obj.num_protons()]
    except KeyError:
        return "silver" # default


Element.color = _color
Element.hex_color = _hex_color
Element.color_name = _color_name
