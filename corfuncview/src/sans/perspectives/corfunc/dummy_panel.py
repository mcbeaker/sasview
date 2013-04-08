import wx
from wx.lib.scrolledpanel import ScrolledPanel
from sans.guiframe.panel_base import PanelBase

class DummyPanel(ScrolledPanel, PanelBase):
    """
    This panel describes proportion of invariants 
    """
    ## Internal nickname for the window, used by the AUI manager
    window_name = "CorFunc"
    ## Name to appear on the window title bar
    window_caption = "CorFunc"
    CENTER_PANE = True
    
    def __init__(self, parent, id=-1, plots=None, 
                 standalone=False, **kwargs):
        """
        """
        kwargs["size"] = (100, 100)
        kwargs["style"] = wx.FULL_REPAINT_ON_RESIZE
        ScrolledPanel.__init__(self, parent, id=id, **kwargs)
        PanelBase.__init__(self, parent)
        self.SetupScrolling()
