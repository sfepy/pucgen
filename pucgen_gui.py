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
from functools import reduce

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

def OverwriteBox(filename):
    msg = QMessageBox()
    msg.setIcon(QMessageBox.Warning)
    msg.setText('The file "%s" already exists. Do you wish to overwrite it?' % filename)
    msg.setWindowTitle('Overwrite File?')
    msg.setStandardButtons(QMessageBox.Ok | QMessageBox.Cancel)
    btn_ok = msg.button(QMessageBox.Ok)
    btn_ok.setText('Overwrite')
    msg.exec_()

    if msg.clickedButton() == btn_ok:
        return True
    else:
        return False

def select_vtk_file(parent, s, le, default=''):
    dlg = QFileDialog(parent, s, default, 'VTK(*.vtk)')
    dlg.setLabelText(QFileDialog.Accept, 'Choose')
    if dlg.exec_():
        le.setText(dlg.selectedFiles()[0])
        le.setStyleSheet('background-color: white')

check_edit_types = {
    'np': ['mat_id', 'radius', 'length', 'thickness',
          'es_dmin', 'es_dmax', 'es_in', 'size', 'size_x', 'el_size'],
    'n3p': ['radius', 'size', 'grid'],
    'n3': ['central_point', 'direction'],
    's': ['direction'],
    'o': ['filename_in', 'filename_out'],
}

def check_edits(edits):
    pars = {}
    ok_ = True

    for k in edits.keys():
        ok = True
        val = edits[k].text()
        if k in check_edit_types['o'] and len(val) > 0:
            pars[k] = val
        elif k in check_edit_types['s'] and val.isalpha():
            pars[k] = val
        else:
            try:
                pars[k] = literal_eval(val)
            except:
                ok = False
   
            if ok:
                if k in check_edit_types['np']\
                    and isinstance(pars[k], (int, float))\
                    and pars[k] > 0:
                    pass
                elif k in check_edit_types['n3p']\
                    and isinstance(pars[k], tuple) and len(pars[k]) == 3\
                    and reduce(lambda a, b: a and b,
                               map(lambda x: x>=0, pars[k])):
                    pass
                elif  k in check_edit_types['n3']\
                    and isinstance(pars[k], tuple) and len(pars[k]) == 3:
                    pass
                else:
                    ok = False

        if ok:
            edits[k].setStyleSheet('background-color: white')
        else:
            edits[k].setStyleSheet('background-color: red')

        ok_ = ok_ and ok

    return pars, ok_

class SetParamsDialog(QDialog):
    def __init__(self, parent, name, params):
        super(SetParamsDialog, self).__init__(parent=parent)

        self.edits = {}
        params_tr_dict = {
            'mat_id': 'material id',
            'central_point': 'central point',
            'el_size': 'element size',
        }

        self.setWindowTitle('Edit %s parameters' % name)
        self.vbox = QVBoxLayout()
        for nm, val in params:
            hbox = QHBoxLayout()
            hbox.addWidget(QLabel('%s:' % params_tr_dict.get(nm, nm)))
            le = QLineEdit(val)
            self.edits[nm] = le  
            hbox.addStretch(1)
            hbox.addWidget(le)
            self.vbox.addLayout(hbox)

        btns = QDialogButtonBox.Ok | QDialogButtonBox.Cancel
        btnbox = QDialogButtonBox(btns)
        btnbox.accepted.connect(self.check_values)
        btnbox.rejected.connect(self.reject)
        self.vbox.addWidget(btnbox)
        self.setLayout(self.vbox)

    def check_values(self):
        self.params_dict, ok = check_edits(self.edits)
        if ok:
            return self.accept()

class MainWindow(QMainWindow):

    def __init__(self):
        QMainWindow.__init__(self)

        self.base_cell_mat_id = None
        self.next_component_id = 1 
        self.components = []

        self.initUI()

    def init_RepeaterTab(self):
        vbox = QVBoxLayout()
        vbox.setSpacing(10)

        vbox = QVBoxLayout()

        hbox1 = QHBoxLayout()
        hbox1.addWidget(QLabel('Input VTK file:'))
        hbox1.addStretch(1)
        vbox.addLayout(hbox1)

        hbox2 = QHBoxLayout()
        self.repeter_in_file = QLineEdit('')
        hbox2.addWidget(self.repeter_in_file)
        btn_select_in_file = QPushButton('...', self)
        btn_select_in_file.clicked.connect(self.repeater_select_in_file)
        width = btn_select_in_file.fontMetrics().boundingRect('...').width() + 20
        btn_select_in_file.setMaximumWidth(width)
        hbox2.addWidget(btn_select_in_file)
        vbox.addLayout(hbox2)

        hbox1b = QHBoxLayout()
        hbox1b.addWidget(QLabel('Output VTK file:'))
        hbox1b.addStretch(1)
        vbox.addLayout(hbox1b)

        hbox2b = QHBoxLayout()
        self.repeter_out_file = QLineEdit('')
        hbox2b.addWidget(self.repeter_out_file)
        btn_select_out_file = QPushButton('...', self)
        btn_select_out_file.clicked.connect(self.repeater_select_out_file)
        width = btn_select_out_file.fontMetrics().boundingRect('...').width() + 20
        btn_select_out_file.setMaximumWidth(width)
        hbox2b.addWidget(btn_select_out_file)
        vbox.addLayout(hbox2b)

        hbox3 = QHBoxLayout()
        hbox3.addWidget(QLabel('Grid:'))
        grid = QLineEdit('1, 1, 1')
        hbox3.addWidget(grid)
        hbox3.addStretch(1)
        vbox.addLayout(hbox3)

        hbox4 = QHBoxLayout()
        hbox4.addWidget(QLabel('Size_x:'))
        size_x = QLineEdit('1')
        hbox4.addWidget(size_x)
        hbox4.addStretch(1)
        vbox.addLayout(hbox4)

        hbox5 = QHBoxLayout()
        hbox5.addStretch(1)
        btn_gen_grid = QPushButton('Generate grid', self)
        btn_gen_grid.clicked.connect(self.repeater_generate_grid)
        hbox5.addWidget(btn_gen_grid)
        hbox5.addStretch(1)

        vbox.addStretch(1)
        vbox.addLayout(hbox5)
        vbox.addStretch(1)

        self.edits_to_check = {'grid': grid, 'size_x': size_x,
                               'filename_in': self.repeter_in_file,
                               'filename_out': self.repeter_out_file}

        return vbox

    def repeater_generate_grid(self):
        from gen_mesh_utils import repeater

        pars, ok = check_edits(self.edits_to_check)

        if ok:
            out_file = pars['filename_out']
            if os.path.exists(out_file):
                ok = ok & OverwriteBox(os.path.split(out_file)[1])

        if ok:
            repeater(pars['filename_in'],
                     out_file,
                     pars['grid'], pars['size_x'])
            viewer = VTKViewer(self, out_file, mat_id=None)
            viewer.exec_()

    def repeater_select_in_file(self):
        select_vtk_file(self, 'Input VTK file', self.repeter_in_file)

    def repeater_select_out_file(self):
        select_vtk_file(self, 'Output VTK file', self.repeter_out_file)

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

        hbox = QHBoxLayout()
        hbox.addStretch(1)
        hbox.addLayout(vbox1)
        hbox.addStretch(1)
        hbox.addLayout(vbox2)
        hbox.addStretch(1)

        vbox.addStretch(1)
        vbox.addLayout(hbox)
        vbox.addStretch(1)

        hbox = QHBoxLayout()
        hbox.addWidget(QLabel('Output VTK file:'))
        hbox.addStretch(1)
        vbox.addLayout(hbox)

        hbox = QHBoxLayout()
        self.generator_out_file = QLineEdit('')
        hbox.addWidget(self.generator_out_file)
        btn_select_in_file = QPushButton('...', self)
        btn_select_in_file.clicked.connect(self.generator_select_out_file)
        width = btn_select_in_file.fontMetrics().boundingRect('...').width() + 20
        btn_select_in_file.setMaximumWidth(width)
        hbox.addWidget(btn_select_in_file)
        vbox.addLayout(hbox)

        hbox = QHBoxLayout()
        hbox.addStretch(1)
        btn_save = QPushButton('Save', self)
        btn_save.clicked.connect(self.save_puc)
        hbox.addWidget(btn_save)
        btn_load = QPushButton('Load', self)
        btn_load.clicked.connect(self.load_puc)
        hbox.addWidget(btn_load)
        btn_generate = QPushButton('Generate', self)
        btn_generate.clicked.connect(self.generate)
        hbox.addWidget(btn_generate)
        hbox.addStretch(1)

        vbox.addLayout(hbox)

        self.new_component(base_cell=True)

        return vbox

    def generator_select_out_file(self):
        select_vtk_file(self, 'Output VTK file', self.generator_out_file)

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

        scpars = []
        for k, _ in clsargs:
            val = pars[k]
            sval = ', '.join(map(str, val)) if hasattr(val, '__len__')\
                else str(val)
            scpars.append([k, sval])

        dlg = SetParamsDialog(self, cls.__name__, scpars)
        if dlg.exec_():
            for k in pars.keys():
                 pars[k] = dlg.params_dict[k]

            self.listbox_update(selected=idx)

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

        tab3 = QWidget()
        tab3.setLayout(self.init_RepeaterTab())
        tabs.addTab(tab3, 'Repeater')

        vbox.addWidget(tabs)

        # save, load, generate, quit
        hbox = QHBoxLayout()
        hbox.addStretch(1)
        btn_quit = QPushButton('Quit', self)
        btn_quit.clicked.connect(self.quit)
        hbox.addWidget(btn_quit)
        hbox.addStretch(1)

        vbox.addLayout(hbox)

        cw.setLayout(vbox)

        self.setWindowTitle('PUCGen')
        self.show()

    def quit(self, event):
        self.close()

    def generate(self):
        pars, ok = check_edits({'filename_out': self.generator_out_file})

        if ok:
            out_file = pars['filename_out']
            if os.path.exists(out_file):
                ok = ok & OverwriteBox(os.path.split(out_file)[1])

        if ok and len(out_file) > 0:
            self.puc = PUC(cell_mat_id=None)
            for cls, pars, act in self.components:
                if act:
                    self.puc.add(cls(**pars))

            failed = False
            try:
                self.puc(out_file)
            except:
                MessageBox('Gmsh error!', 'error')
                failed = True

            if not failed:
                viewer = VTKViewer(self, out_file,
                                   self.components[0][1].get('mat_id'))
                viewer.exec_()

def main():
    app = QApplication(sys.argv)
    mw = MainWindow()
    sys.exit(app.exec_())

if __name__ == "__main__":
    main()