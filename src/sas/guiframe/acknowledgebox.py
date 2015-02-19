'''
Created on Feb 18, 2015

@author: jkrzywon
'''

__id__ = "$Id: acknoweldgebox.py 2015-18-02 jkrzywon $"
__revision__ = "$Revision: 1193 $"

import wx
import wx.lib.hyperlink
import random
import os.path
import os
try:
    # Try to find a local config
    import imp
    path = os.getcwd()
    if(os.path.isfile("%s/%s.py" % (path, 'local_config'))) or \
      (os.path.isfile("%s/%s.pyc" % (path, 'local_config'))):
        fObj, path, descr = imp.find_module('local_config', [path])
        config = imp.load_module('local_config', fObj, path, descr)  
    else:
        # Try simply importing local_config
        import local_config as config
except:
    # Didn't find local config, load the default 
    import config
    

class DialogAcknowledge(wx.Dialog):
    """
    "Acknowledgement" Dialog Box
    
    Shows the current method for acknowledging SasView in
    scholarly publications.
    
    """
    
    def __init__(self, *args, **kwds):
        
        kwds["style"] = wx.DEFAULT_DIALOG_STYLE
        wx.Dialog.__init__(self, *args, **kwds)
        
        self.ack = wx.TextCtrl(self, style=wx.TE_LEFT|wx.TE_MULTILINE|wx.TE_BESTWRAP|wx.TE_READONLY|wx.TE_NO_VSCROLL)
        self.ack.SetBackgroundColour(wx.NullColour)
        self.ack.SetValue(config._acknowledgement_publications)
        self.ack.SetMinSize((-1,60))
        self.preamble = wx.StaticText(self, -1, config._acknowledgement_preamble)
        self.static_line = wx.StaticLine(self, 0)
        
        self.__set_properties()
        self.__do_layout()
        
    def __set_properties(self):
        """
        """
        # begin wxGlade: DialogAbout.__set_properties
        self.SetTitle("Acknowledging SasView")
        self.SetSize((500, 225))
        # end wxGlade

    def __do_layout(self):
        """
        """
        # begin wxGlade: DialogAbout.__do_layout
        sizer_main = wx.BoxSizer(wx.VERTICAL)
        sizer_titles = wx.BoxSizer(wx.VERTICAL)
        sizer_titles.Add(self.preamble, 0, wx.ALL|wx.EXPAND, 0)
        sizer_titles.Add(self.static_line, 0, wx.EXPAND, 0)
        sizer_titles.Add(self.ack, 0, wx.ALL|wx.EXPAND, 0)
        sizer_main.Add(sizer_titles, -1, wx.ALL|wx.EXPAND, 5)
        self.SetAutoLayout(True)
        self.SetSizer(sizer_main)
        self.Layout()
        self.Centre()
        # end wxGlade
        

##### testing code ############################################################
class MyApp(wx.App):
    """
    """
    def OnInit(self):
        """
        """
        wx.InitAllImageHandlers()
        dialog = DialogAcknowledge(None, -1, "")
        self.SetTopWindow(dialog)
        dialog.ShowModal()
        dialog.Destroy()
        return 1

# end of class MyApp

if __name__ == "__main__":
    app = MyApp(0)
    app.MainLoop()