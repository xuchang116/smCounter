import sys
import datetime
import logging

# stream handler object that redirects to Python logging facility
class RedirectToLogger(object):
   def __init__(self):
      self.logger = logging.getLogger()
 
   def write(self, buf):
      for line in buf.rstrip().splitlines():
         self.logger.debug(line.rstrip())

   def flush(self):
      pass
 
# initialize this logging module
def init(logToConsole, logFileName):

   # format a log file name
   now = datetime.datetime.now()
   timestamp = now.strftime("_%Y.%m.%d_%H.%M.%S")

   # set up logging to a log file
   logDateFormat = "%Y-%m-%d %H:%M:%S"
   logFormat = "%(asctime)s.%(msecs)03d %(message)s"
   logging.basicConfig(level=logging.DEBUG,
                       format=logFormat,
                       datefmt=logDateFormat,
                       filename=logFileName,
                       filemode='w')

   # create a stream redirection object for stdout, so print() goes to logger
   rtl = RedirectToLogger()
   sys.stdout = rtl
   
   # if not logging to console, redirect stderr to the log file
   if not logToConsole:
      sys.stderr = rtl
      
   # if logging to console route: stdout==> logger ==> stderr & logfile
   else:
      console = logging.StreamHandler()
      console.setLevel(logging.DEBUG)
      formatter = logging.Formatter(logFormat,logDateFormat)
      console.setFormatter(formatter)
      logging.getLogger().addHandler(console)
