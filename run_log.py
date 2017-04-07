import sys
import datetime
import logging

# global hack
sysStdOut = sys.stdout
sysStdErr = sys.stderr

#----------------------------------------------------------------------
# stream handler object that redirects to Python logging facility
#----------------------------------------------------------------------
class RedirectToLogger(object):
   def __init__(self):
      self.logger = logging.getLogger()
 
   def write(self, buf):
      for line in buf.rstrip().splitlines():
         self.logger.debug(line.rstrip())
         
   def flush(self):
      pass
 
#----------------------------------------------------------------------
# initialize the logging
#----------------------------------------------------------------------
def init(logFilePrefix):

   # format a log file name
   now = datetime.datetime.now()
   timestamp = now.strftime("_%Y.%m.%d_%H.%M.%S")
   logFileName = logFilePrefix + ".run-log" + timestamp + ".txt"

   # set up logging to a log file
   logDateFormat = "%Y-%m-%d %H:%M:%S"
   logFormat = "%(asctime)s.%(msecs)03d %(message)s"
   logging.basicConfig(level=logging.DEBUG,
                       format=logFormat,
                       datefmt=logDateFormat,
                       filename=logFileName,
                       filemode="w")

   # create a stream redirection object for stdout, so print() goes to logger
   rtl = RedirectToLogger()
   sys.stdout = rtl
   sys.stderr = rtl

#----------------------------------------------------------------------
# close log file and restore stdout and stderr
#----------------------------------------------------------------------
def close():
   sys.stdout = sysStdOut
   sys.stdout = sysStdErr  
   logger = logging.getLogger()
   handler = logger.handlers[0]
   handler.stream.close()
   logger.removeHandler(handler)
