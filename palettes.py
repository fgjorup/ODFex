from PyQt5 import QtGui

class container(object):
    pass

def get_defaults_thm():
    ###################
    #  Theme Details  #
    # From xrdPlanner #
    ###################
    thm = container()
    thm.color_dark = '#404040'                    # [color]  Global dark color
    thm.color_light = '#EEEEEE'                   # [color]  Global light color
    # light mode
    thm.light_conic_label_fill = '#FFFFFF'        # [color]  Contour label fill color
    thm.light_conic_ref_color = '#DCDCDC'         # [color]  Reference contour color
    thm.light_beamstop_color = '#FF000080'        # [color]  Beamstop color
    thm.light_beamstop_edge_color = '#FF0000'     # [color]  Beamstop edge color
    thm.light_det_module_color = '#404040'        # [color]  Detector module border color
    thm.light_det_module_fill = '#404040'         # [color]  Detector module background color
    thm.light_plot_bg_color = '#FFFFFF'           # [color]  Plot background color
    thm.light_unit_label_color = '#808080'        # [color]  Label color
    thm.light_unit_label_fill = '#FFFFFF'         # [color]  Label fill color
    thm.light_slider_border_color = '#808080'     # [color]  Slider frame border color
    thm.light_slider_bg_color = '#AAC0C0C0'       # [color]  Slider frame background color
    thm.light_slider_bg_hover = '#C0C0C0'         # [color]  Slider frame hover color
    thm.light_slider_label_color = '#000000'      # [color]  Slider frame label color
    thm.light_overlay_threshold_color = '#FF0000' # [color]  Map threshold color
    # dark mode
    thm.dark_conic_label_fill = '#000000'         # [color]  Contour label fill color
    thm.dark_conic_ref_color = '#303030'          # [color]  Reference contour color
    thm.dark_beamstop_color = '#FF0000AA'         # [color]  Beamstop color
    thm.dark_beamstop_edge_color = '#FF0000'      # [color]  Beamstop edge color
    thm.dark_det_module_color = '#EEEEEE'         # [color]  Detector module border color
    thm.dark_det_module_fill = '#EEEEEE'          # [color]  Detector module background color
    thm.dark_plot_bg_color = '#000000'            # [color]  Plot background color
    thm.dark_unit_label_color = '#C0C0C0'         # [color]  Label color
    thm.dark_unit_label_fill = '#000000'          # [color]  Label fill color
    thm.dark_slider_border_color = '#202020'      # [color]  Slider frame border color
    thm.dark_slider_bg_color = '#AA303030'        # [color]  Slider frame background color
    thm.dark_slider_bg_hover = '#303030'          # [color]  Slider frame hover color
    thm.dark_slider_label_color = '#C0C0C0'       # [color]  Slider frame label color
    thm.dark_overlay_threshold_color = '#FF0000'  # [color]  Map threshold color

    return thm

def get_palette(use_dark):
    
    thm = get_defaults_thm()

    _color_dark = QtGui.QColor(thm.color_dark)
    _color_light = QtGui.QColor(thm.color_light)

    _red = QtGui.QColor('#FF0000')
    _green = QtGui.QColor('#00FF00')
    _white = QtGui.QColor('#FFFFFF')

    # # set highlight text color depending on the lightness of the colormap
    # _highlight_color = self.cont_cmap.map(0.0, mode='qcolor')
    # if _highlight_color.lightnessF() < 0.5:
    #     _highlight_text = _color_light.lighter(150)
    # else:
    #     _highlight_text = _color_dark.darker(150)
    # define color palette
    if use_dark:
        # palette
        palette = QtGui.QPalette()
        palette.setColor(QtGui.QPalette.ColorRole.Window,          _color_dark)
        palette.setColor(QtGui.QPalette.ColorRole.WindowText,      _white)#_color_light)
        palette.setColor(QtGui.QPalette.ColorRole.Button,          _color_dark)
        palette.setColor(QtGui.QPalette.ColorRole.ButtonText,      _color_light)
        palette.setColor(QtGui.QPalette.ColorRole.Base,            _color_dark)
        palette.setColor(QtGui.QPalette.ColorRole.AlternateBase,   _color_dark.lighter(110))
        palette.setColor(QtGui.QPalette.ColorRole.Highlight,       _color_light.darker(150))#_highlight_color)
        palette.setColor(QtGui.QPalette.ColorRole.HighlightedText, _color_dark.darker(150))#_highlight_text)
        palette.setColor(QtGui.QPalette.ColorRole.Text,            _color_light)
        palette.setColor(QtGui.QPalette.ColorRole.PlaceholderText, _color_light)
        palette.setColor(QtGui.QPalette.ColorRole.BrightText,      _color_light)
        palette.setColor(QtGui.QPalette.ColorRole.ToolTipBase,     _color_dark)
        palette.setColor(QtGui.QPalette.ColorRole.ToolTipText,     _color_light)
    else:
        # palette
        palette = QtGui.QPalette()
        palette.setColor(QtGui.QPalette.ColorRole.Window,          _color_light)
        palette.setColor(QtGui.QPalette.ColorRole.WindowText,      _color_dark)
        palette.setColor(QtGui.QPalette.ColorRole.Button,          _color_light)
        palette.setColor(QtGui.QPalette.ColorRole.ButtonText,      _color_dark)
        palette.setColor(QtGui.QPalette.ColorRole.Base,            _color_light)
        palette.setColor(QtGui.QPalette.ColorRole.AlternateBase,   _color_light.darker(110))
        palette.setColor(QtGui.QPalette.ColorRole.Highlight,       _color_dark)
        palette.setColor(QtGui.QPalette.ColorRole.HighlightedText, _color_light.lighter(150))
        palette.setColor(QtGui.QPalette.ColorRole.Text,            _color_dark)
        palette.setColor(QtGui.QPalette.ColorRole.PlaceholderText, _color_dark)
        palette.setColor(QtGui.QPalette.ColorRole.BrightText,      _color_dark)
        palette.setColor(QtGui.QPalette.ColorRole.ToolTipBase,     _color_light)
        palette.setColor(QtGui.QPalette.ColorRole.ToolTipText,     _color_dark)

    return palette

    # # apply palette to app
    # app = QtWidgets.QApplication.instance()
    # app.setPalette(palette)
    # app.setStyle('Fusion')