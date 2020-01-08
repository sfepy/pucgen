#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
PUCGEN - GUI for Periodic Unit Cell generator

"""

import sys
import os
from optparse import OptionParser
from inspect import getargspec
from copy import deepcopy
from ast import literal_eval

from PyQt5.QtWidgets import QApplication, QMainWindow, QWidget,\
     QHBoxLayout, QVBoxLayout, QTabWidget, QLineEdit,\
     QLabel, QPushButton, QFrame, QFileDialog, QMessageBox,\
     QComboBox, QListWidget, QDialog, QDialogButtonBox
from PyQt5.QtGui import QFont, QPixmap, QPalette, QColor

from pucgen import PUC, pucgen_classes
from vtk_viewer import VTKViewer

params_dict = {}


def MessageBox(s, type='warning'):
    title, qmb = {
        'warning': ('Warning', QMessageBox.Warning),
        'error': ('Error', QMessageBox.Critical)
    }[type]

    msg = QMessageBox()
    msg.setIcon(qmb)
    msg.setText(s)
    msg.setWindowTitle(title)
    msg.setStandardButtons(QMessageBox.Ok)
    msg.exec_()


class SetParamsDialog(QDialog):
    def __init__(self, parent, name, params, err_lines=[]):
        super(SetParamsDialog, self).__init__(parent=parent)
        self.setWindowTitle('Edit %s parameters' % name)

        self.vbox = QVBoxLayout()
        self.params_dict = {}
        params_tr_dict = {
            'mat_id': 'material id',
            'central_point': 'central point',
            'el_size': 'element size',
        }
        for ii, (nm, val) in enumerate(params):
            hbox = QHBoxLayout()
            hbox.addWidget(QLabel('%s:' % params_tr_dict.get(nm, nm)))
            le = QLineEdit(val)
            if ii in err_lines:
                le.setStyleSheet('background-color: red')
            self.params_dict[nm] = le  
            hbox.addStretch(1)
            hbox.addWidget(le)
            self.vbox.addLayout(hbox)

        btns = QDialogButtonBox.Ok | QDialogButtonBox.Cancel
        btnbox = QDialogButtonBox(btns)
        btnbox.accepted.connect(self.accept)
        btnbox.rejected.connect(self.reject)
        self.vbox.addWidget(btnbox)
        self.setLayout(self.vbox)


class MainWindow(QMainWindow):

    def __init__(self):
        QMainWindow.__init__(self)

        self.base_cell_mat_id = None
        self.next_component_id = 1 
        self.components = []

        self.initUI()

    def init_GeneratorTab(self):

        self.gen_classes = pucgen_classes
        self.class_args = {}
        clslist = []

        for cls in self.gen_classes:
            args, _, _, defaults = getargspec(cls.__init__)

            clsargs = zip(args[1:], defaults)
            self.class_args[cls.__name__] = clsargs
            clslist.append(cls.__name__)

        vbox = QVBoxLayout()
        vbox.setSpacing(10)

        vbox1 = QVBoxLayout()
        btn_new = QPushButton('New', self)
        btn_new.clicked.connect(self.new_component)
        self.cmb_new = QComboBox()
        self.cmb_new.addItems(clslist[1:])
        btn_edit = QPushButton('Edit', self)
        btn_edit.clicked.connect(self.edit_component)
        btn_delete = QPushButton('Delete', self)
        btn_delete.clicked.connect(self.delete_component)
        btn_active = QPushButton('Activate/Deactivate', self)
        btn_active.clicked.connect(self.change_active)
        vbox1.addWidget(btn_new)
        vbox1.addWidget(self.cmb_new)
        hr = QFrame()
        hr.setFrameShape(QFrame.HLine)
        vbox1.addWidget(hr)
        vbox1.addWidget(btn_edit)
        vbox1.addWidget(btn_delete)
        vbox1.addWidget(btn_active)
        vbox1.addStretch(1)

        vbox2 = QVBoxLayout()
        self.listbox = QListWidget()
        self.listbox_update(selected=0)
        vbox2.addWidget(self.listbox)
        # vbox2.addStretch(1)

        hbox = QHBoxLayout()
        hbox.addStretch(1)
        hbox.addLayout(vbox1)
        hbox.addStretch(1)
        hbox.addLayout(vbox2)
        hbox.addStretch(1)

        vbox.addStretch(1)
        vbox.addLayout(hbox)
        vbox.addStretch(1)

        self.new_component(base_cell=True)

        return vbox

    def listbox_update(self, selected=0):
        self.listbox.clear()
        for ii, (cls, pars, act) in enumerate(self.components):
            flag = ' ' if act else '#'
            self.listbox.addItem('%s%s (%d)' % (flag, cls.__name__,
                                                pars['mat_id']))
            if not act:
                self.listbox.item(ii).setForeground(QColor('lightGray'))

        if selected is not None:
            self.listbox.setCurrentRow(selected)

    def new_component(self, **kwargs):
        base_cell = kwargs.get('base_cell', False)
        idx = 0 if base_cell else self.cmb_new.currentIndex() + 1
        cls = self.gen_classes[idx]
        kwargs = {k: deepcopy(v) for k, v in self.class_args[cls.__name__]}
        kwargs['mat_id'] = self.next_component_id

        self.components.append([cls, kwargs, True])

        self.next_component_id += 1
        self.listbox_update(selected=len(self.components) - 1)

    def edit_component(self):
        idx = self.listbox.currentRow()
        cls, pars, act = self.components[idx]
        clsargs = self.class_args[cls.__name__]

        status = 0
        err_lines = []
        pars0 = deepcopy(pars)
        scpars = []
        for k, _ in clsargs:
            val = pars[k]
            sval = ', '.join(map(str, val)) if hasattr(val, '__len__')\
                else str(val)
            scpars.append([k, sval])

        while not(status):
            dlg = SetParamsDialog(self, cls.__name__, scpars,
                                  err_lines)
            if dlg.exec_():
                status = 1
                ncpars = dlg.params_dict
                err_lines = []
                for ii, (k, _) in enumerate(clsargs):
                    val = ncpars[k].text()
                    scpars[ii][1] = val
                    if val.isalpha():
                        pars[k] = val
                    else:
                        try:
                            pars[k] = literal_eval(val)
                        except SyntaxError:
                            status = 0
                            err_lines.append(ii)

            else:
                status = -1

        if status > 0:
            self.listbox_update(selected=idx)
        else:
            for k, v in pars0.items():
                pars[k] = v
            # WarningBox('Parameters of %s not changed!' % clsname)

    def delete_component(self):
        idx = self.listbox.currentRow()
        if idx == 0:
            MessageBox('Can not delete Base Cell!')
        else:
            self.components.pop(idx)
            self.listbox_update(selected=idx - 1)

    def change_active(self):
        idx = self.listbox.currentRow()
        if idx == 0:
            MessageBox('Can not deactivate Base Cell!')
        else:
            self.components[idx][2] = not(self.components[idx][2])
            self.listbox_update(selected=idx)

    def save_puc(self, **kwargs):
        fname = kwargs.get('fname', None)
        if fname is None:
            fname, _ = QFileDialog.getSaveFileName(self, 'Save PUC file',
                                                   filter='Files (*.puc)')

        if len(fname) > 0:
            PUC.save_puc(fname, self.components)

    def load_puc(self, **kwargs):
        fname = kwargs.get('fname', None)
        if fname is None:
            fname, _ = QFileDialog.getOpenFileName(self, 'Load PUC file',
                                                   filter='Files (*.puc)')

        if len(fname) > 0:
            self.components = PUC.load_puc(fname)
            self.listbox_update()

    def initUI(self):
        import os.path as op
        path_to_script = op.dirname(os.path.abspath(__file__))

        cw = QWidget()
        self.setCentralWidget(cw)
        vbox = QVBoxLayout()
        vbox.setSpacing(10)

        vbox1 = QVBoxLayout()
        font_info = QFont()
        font_info.setItalic(True)
        font_info.setPixelSize(12)
        info = QLabel('Version: 0.1')
        info.setFont(font_info)
        vbox1.addStretch(1)
        vbox1.addWidget(info)

        vbox2 = QVBoxLayout()
        pucgen_logo = QLabel()
        logopath = os.path.join(path_to_script, 'pucgen_logo.png')
        logo = QPixmap(logopath)
        pucgen_logo.setPixmap(logo)
        vbox2.addWidget(pucgen_logo)

        hbox = QHBoxLayout()
        hbox.addStretch(1)
        hbox.addLayout(vbox2)
        hbox.addLayout(vbox1)
        hbox.addStretch(1)
        vbox.addLayout(hbox)

        tabs = QTabWidget()
        tab1 = QWidget()

        tab1.setLayout(self.init_GeneratorTab())
        tabs.addTab(tab1, 'Generator')
        vbox.addWidget(tabs)

        # save, load, generate, quit
        hbox = QHBoxLayout()
        hbox.addStretch(1)
        btn_save = QPushButton('Save', self)
        btn_save.clicked.connect(self.save_puc)
        btn_load = QPushButton('Load', self)
        btn_load.clicked.connect(self.load_puc)
        btn_generate = QPushButton('Generate', self)
        btn_generate.clicked.connect(self.generate)
        btn_quit = QPushButton('Quit', self)
        btn_quit.clicked.connect(self.quit)
        hbox.addWidget(btn_save)
        hbox.addWidget(btn_load)
        hbox.addWidget(btn_generate)
        hbox.addWidget(btn_quit)
        hbox.addStretch(1)

        vbox.addLayout(hbox)

        cw.setLayout(vbox)
        self.setWindowTitle('PUCGen')
        self.show()

    def quit(self, event):
        self.close()

    def generate(self, **kwargs):
        fname = kwargs.get('fname', None)
        if fname is None:
            fname, _ = QFileDialog.getSaveFileName(self, 'Save VTK file',
                                                   filter='Files (*.vtk)')

        if len(fname) > 0:
            self.puc = PUC(cell_mat_id=None)
            for cls, pars, act in self.components:
                if act:
                    self.puc.add(cls(**pars))

            failed = False
            try:
                self.puc(fname)
            except:
                MessageBox('Gmsh error!', 'error')
                failed = True

            if not failed:
                viewer = VTKViewer(self, fname,
                                   self.components[0][1].get('mat_id'))
                viewer.exec_()

def main():
    app = QApplication(sys.argv)
    mw = MainWindow()
    sys.exit(app.exec_())

if __name__ == "__main__":
    main()