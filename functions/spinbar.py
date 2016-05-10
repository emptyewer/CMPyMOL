import threading
import time

class SpinCursor(threading.Thread):
    """ A console spin cursor class """
    def __init__(self, parent, msg='Please Wait...', speed=3, ):
        self.parent = parent
        self.flag = False
        # Any message to print first ?
        self.msg = msg
        # Complete printed string
        self.string = ''
        # Speed is given as number of spins a second
        # Use it to calculate spin wait time
        self.waittime = 1.0/float(speed * 4)
        self.spinchars = ('-', '\\ ', '| ', '/ ')
        threading.Thread.__init__(self, None, None, "Spin Thread")
        
    def spin(self):
        """ Perform a single spin """
        for x in self.spinchars:
            self.string = self.msg + "\t" + x + "\r"
            self.parent.status_bar.showMessage(self.string.encode('utf-8'))
            time.sleep(self.waittime)

    def run(self):
        while (not self.flag):
            self.spin()
        self.parent.status_bar.showMessage(" " * (len(self.string) + 1))
        
    def stop(self):
        self.parent.status_bar.showMessage('')
        self.flag = True

# class StatusSpinCursor(threading.Thread):
#     """ A console spin cursor class """
#     def __init__(self, status, msg='Please Wait...', speed=3):
#         self.out = sys.stdout
#         self.flag = False
#         self.status = status
#         # Any message to print first ?
#         self.msg = msg
#         # Complete printed string
#         self.string = ''
#         # Speed is given as number of spins a second
#         # Use it to calculate spin wait time
#         self.waittime = 1.0/float(speed * 4)
#         self.spinchars = ('-', '\\ ', '| ', '/ ')
#         threading.Thread.__init__(self, None, None, "Spin Thread")
#
#     def spin(self):
#         """ Perform a single spin """
#         for x in self.spinchars:
#             self.string = self.msg + "\t" + x + "\r"
#             self.status.showMessage(self.string.encode('utf-8'))
#             time.sleep(self.waittime)
#
#     def run(self):
#         while (not self.flag):
#             self.spin()
#         self.out.write(" " * (len(self.string) + 1))
#
#     def stop(self):
#         self.flag = True

