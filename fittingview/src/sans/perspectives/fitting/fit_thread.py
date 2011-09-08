
import sys
from data_util.calcthread import CalcThread

def map_getattr(classInstance, classFunc, *args):
    """
    Take an instance of a class and a function name as a string.
    Execute class.function and return result
    """
    return  getattr(classInstance,classFunc)(*args)

def map_apply(arguments):
    return apply(arguments[0], arguments[1:])

class FitThread(CalcThread):
    """Thread performing the fit """
    
    def __init__(self, 
                  fn,
                  page_id,
                   handler,
                  pars=None,
                 completefn = None,
                 updatefn   = None,
                 yieldtime  = 0.01,
                 worktime   = 0.01,
                 ftol       = None):
        CalcThread.__init__(self,completefn,
                 updatefn,
                 yieldtime,
                 worktime)
        self.handler = handler
        self.fitter = fn
        self.pars = pars
        self.page_id = page_id
        self.starttime = 0
        self.updatefn = updatefn
        #Relative error desired in the sum of squares.
        self.ftol = ftol
   
    def isquit(self):
        """
        :raise KeyboardInterrupt: when the thread is interrupted
        
        """
        try:
            CalcThread.isquit(self)
        except KeyboardInterrupt:
            raise KeyboardInterrupt
       
    def compute(self):
        """
        Perform a fit 
        """
        msg = ""
        try:
            list_handler = []
            list_curr_thread = [] 
            list_ftol = []
            list_map_get_attr = []
            list_fit_function = []
            list_q = []
            for i in range(len(self.fitter)):
                list_handler.append(None)
                list_q.append(None)
                list_curr_thread.append(None)
                list_ftol.append(self.ftol)
                list_fit_function.append('fit')
                list_map_get_attr.append(map_getattr)
            from multiprocessing import Pool
            inputs = zip(list_map_get_attr,self.fitter, list_fit_function,
                         list_handler, list_q, list_curr_thread,list_ftol)
            result =  Pool(1).map(func=map_apply, 
                               iterable=inputs)
            #self.handler.starting_fit()
            self.complete(result= result,
                          page_id=self.page_id,
                          pars = self.pars)
           
        except KeyboardInterrupt, msg:
            # Thread was interrupted, just proceed and re-raise.
            # Real code should not print, but this is an example...
            #print "keyboard exception"
            #Stop on exception during fitting. Todo: need to put 
            #some mssg and reset progress bar.
            raise
            #if self.handler is not None:
            #    self.handler.error(msg=msg)
        except:
            if self.handler is not None:
                self.handler.error(msg=str(sys.exc_value))
           
        
    