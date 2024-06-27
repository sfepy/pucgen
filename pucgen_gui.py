#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
PUCGEN - GUI for Periodic Unit Cell generator

"""

import sys
import os
from inspect import getargspec
from copy import deepcopy
from ast import literal_eval
from functools import reduce

from PyQt6.QtWidgets import (QApplication, QMainWindow, QWidget,
     QHBoxLayout, QVBoxLayout, QTabWidget, QLineEdit,
     QLabel, QPushButton, QFrame, QFileDialog, QMessageBox,
     QComboBox, QListWidget, QDialog, QDialogButtonBox)
from PyQt6.QtGui import QFont, QPixmap, QPalette, QColor

from pucgen import PUC, pucgen_classes, version
from vtk_viewer import VTKViewer

params_dict = {}


def MessageBox(s, type='warning'):
    title, qmb = {
        'warning': ('Warning', QMessageBox.Icon.Warning),
        'error': ('Error', QMessageBox.Icon.Critical)
    }[type]

    msg = QMessageBox()
    msg.setIcon(qmb)
    msg.setText(s)
    msg.setWindowTitle(title)
    msg.setStandardButtons(QMessageBox.StandardButton.Ok)
    msg.exec()


def select_vtk_open(parent, s, default=''):
    return QFileDialog.getOpenFileName(parent, s, default, 'VTK(*.vtk)')[0]


def select_vtk_save(parent, s, default=''):
    return QFileDialog.getSaveFileName(parent, s, default, 'VTK(*.vtk)')[0]


check_edit_types = {
    'np': ['dimension', 'mat_id', 'el_size', 'size_x'],
    'n3p': ['grid'],
    'n3': ['central_point', 'direction'],
    'nn': ['dimension'],
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
                if (k in check_edit_types['np']
                    and isinstance(pars[k], (int, float)) and pars[k] > 0):
                    pass
                elif (k in check_edit_types['n3p']
                      and isinstance(pars[k], tuple) and len(pars[k]) == 3
                      and reduce(lambda a, b: a and b,
                                 map(lambda x: x>=0, pars[k]))):
                    pass
                elif (k in check_edit_types['n3']
                      and isinstance(pars[k], tuple) and len(pars[k]) == 3):
                    pass
                elif (k in check_edit_types['nn']
                      and isinstance(pars[k], tuple)
                      and reduce(lambda a, b: a and b,
                        map(lambda x: isinstance(x, (int, float)), pars[k]))):
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

        self.setWindowTitle(f'Edit {name} parameters')
        self.vbox = QVBoxLayout()
        for nm, val in params:
            hbox = QHBoxLayout()
            hbox.addWidget(QLabel(f'{params_tr_dict.get(nm, nm)}:'))
            le = QLineEdit(val)
            self.edits[nm] = le
            hbox.addStretch(1)
            hbox.addWidget(le)
            self.vbox.addLayout(hbox)

        btns = QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel
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
        self.puc = None

        self.initUI()

    def init_RepeaterTab(self):
        vbox = QVBoxLayout()
        vbox.setSpacing(10)

        vbox = QVBoxLayout()

        hbox1 = QHBoxLayout()
        hbox1.addStretch(1)
        btn_select_in_file = QPushButton('Open VTK', self)
        btn_select_in_file.clicked.connect(self.repeater_open_in_file)
        hbox1.addWidget(btn_select_in_file)
        hbox1.addStretch(1)
        vbox.addLayout(hbox1)

        hbox2a = QHBoxLayout()
        self.mesh_info_line1 = QLabel('')
        hbox2a.addWidget(self.mesh_info_line1)
        vbox.addLayout(hbox2a)

        hbox2b = QHBoxLayout()
        self.mesh_info_line2 = QLabel('')
        hbox2b.addWidget(self.mesh_info_line2)
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

        vbox.addStretch(1)

        hbox5 = QHBoxLayout()
        hbox5.addStretch(1)
        btn_gen_grid = QPushButton('Generate grid', self)
        btn_gen_grid.clicked.connect(self.repeater_generate_grid)
        hbox5.addWidget(btn_gen_grid)
        hbox5.addStretch(1)
        vbox.addLayout(hbox5)

        hbox5a = QHBoxLayout()
        self.rep_info_line = QLabel('')
        hbox5a.addWidget(self.rep_info_line)
        vbox.addLayout(hbox5a)
        vbox.addStretch(1)

        self.edits_to_check = {'grid': grid, 'size_x': size_x}

        return vbox

    def repeater_generate_grid(self):
        from gen_mesh_utils import repeat_cell

        fname = select_vtk_save(self, 'Output VTK file')

        pars, ok = check_edits(self.edits_to_check)

        if ok and (len(fname) > 0):
            repeat_cell(self.repeater_mesh, fname,
                        pars['grid'], pars['size_x'])
            self.rep_info_line.setText(f'generated {fname}')
            # viewer = VTKViewer(self, fname, mat_id=None)
            # viewer.exec()

    def repeater_open_in_file(self):
        fname = select_vtk_open(self, 'Input VTK file')

        if len(fname) > 0:
            from gen_mesh_utils import meshio_read
            mesh = meshio_read(fname)
            self.mesh_info_line1.setText(fname)
            text = f'points: {mesh.points.shape[0]}, '
            text += ', '.join([f'{cg.type}: {cg.data.shape[0]}'
                               for cg in mesh.cells])
            self.mesh_info_line2.setText(text)
            self.repeater_mesh = mesh

    def init_GeneratorTab(self):

        self.gen_classes = pucgen_classes
        self.class_args = {}
        clslist = []

        for cls in self.gen_classes:
            args, _, _, defaults = getargspec(cls.__init__)

            clsargs = tuple(zip(args[1:], defaults))
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
        hr.setFrameShape(QFrame.Shape.HLine)
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

        hbox = QHBoxLayout()
        self.gen_info_line = QLabel('')
        hbox.addWidget(self.gen_info_line)
        vbox.addLayout(hbox)

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
        cls, pars, _ = self.components[idx]
        clsargs = self.class_args[cls.__name__]

        scpars = []
        for k, _ in clsargs:
            val = pars[k]
            sval = ', '.join(map(str, val)) if hasattr(val, '__len__')\
                else str(val)
            scpars.append([k, sval])

        dlg = SetParamsDialog(self, cls.__name__, scpars)
        if dlg.exec():
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
            self.gen_info_line.setText(f'saved to {fname}')

    def load_puc(self, **kwargs):
        fname = kwargs.get('fname', None)
        if fname is None:
            fname, _ = QFileDialog.getOpenFileName(self, 'Load PUC file',
                                                   filter='Files (*.puc)')

        if len(fname) > 0:
            self.components = PUC.load_puc(fname)
            self.listbox_update()
            self.gen_info_line.setText(f'loaded {fname}')

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
        info = QLabel(f'Version: {version}')
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
        fname = select_vtk_save(self, 'Output VTK file')

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

            self.gen_info_line.setText(f'generated {fname}')

            # if not failed:
            #     viewer = VTKViewer(self, out_file,
            #                        self.components[0][1].get('mat_id'))
            #     viewer.exec()

def main():
    app = QApplication(sys.argv)
    mw = MainWindow()
    mw.show()
    sys.exit(app.exec())

if __name__ == "__main__":
    main()
