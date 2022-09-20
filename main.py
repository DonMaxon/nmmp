from PyQt5 import uic
from PyQt5.QtWidgets import QApplication
import os
import solving as s

Form, Window = uic.loadUiType('interface.ui')
app = QApplication([])
window = Window()
form = Form()
form.setupUi(window)
window.show()
app.exec_()