import sys
import numpy as nm
import vtk
from PyQt5.QtWidgets import QVBoxLayout, QSpacerItem, QFrame,\
    QSizePolicy, QDialog, QPushButton, QCheckBox, QSlider, QApplication
from PyQt5.QtCore import Qt, QSize

import vtk.qt
vtk.qt.QVTKRWIBase = 'QGLWidget'

from vtk.qt.QVTKRenderWindowInteractor import QVTKRenderWindowInteractor


class VTKViewer(QDialog):
    @staticmethod
    def view_vtk(filename, mat_id=1):
        reader = vtk.vtkUnstructuredGridReader()
        reader.SetFileName(filename)
        reader.Update()

        cell_data = reader.GetOutput().GetCellData()
        array_names = [cell_data.GetArrayName(k)
            for k in range(cell_data.GetNumberOfArrays())]
        if len(array_names) >= 1:
            ca, cb = cell_data.GetArray(0).GetRange()
        else:
            ca, cb = 1, 1
        lut = vtk.vtkLookupTable()
        lut.SetAlphaRange(1, 1)
        lut.SetTableRange(ca, cb)
        mapper = vtk.vtkDataSetMapper()
        mapper.SetScalarModeToUseCellFieldData()
        mapper.SetColorModeToMapScalars()
        mapper.ScalarVisibilityOn()
        if len(array_names) >= 1:
            mapper.SelectColorArray(array_names[0])
        mapper.SetLookupTable(lut)
        mapper.SetScalarRange(ca,cb)
        actor = vtk.vtkActor()
        actor.SetMapper(mapper)
        actor.GetProperty().EdgeVisibilityOn()
        actor.GetProperty().SetOpacity(1.0)

        # Create the Renderer
        renderer = vtk.vtkRenderer()
        renderer.AddActor(actor)

        if mat_id is not None:
            mat_id_field = cell_data.GetArray(0)
            ids = vtk.vtkIdTypeArray()
            ids.SetNumberOfComponents(1)
            for ii in range(mat_id_field.GetNumberOfValues()):
                if mat_id_field.GetValue(ii) == mat_id:
                    ids.InsertNextValue(ii)

            selectionNode = vtk.vtkSelectionNode()
            selectionNode.SetFieldType(0) # CELLS
            selectionNode.SetContentType(4) # 4 INDICES
            selectionNode.SetSelectionList(ids)
            selection = vtk.vtkSelection()
            selection.AddNode(selectionNode)

            extractSelection = vtk.vtkExtractSelection()
            extractSelection.SetInputConnection(0, reader.GetOutputPort())
            extractSelection.SetInputData(1, selection)
            extractSelection.Update()

            part0 = vtk.vtkUnstructuredGrid()
            part0.ShallowCopy(extractSelection.GetOutput())

            selectionNode.GetProperties().Set(vtk.vtkSelectionNode.INVERSE(), 1) #invert the selection
            extractSelection.Update()

            part1 = vtk.vtkUnstructuredGrid()
            part1.ShallowCopy(extractSelection.GetOutput())

            mapper.SetInputData(part1)

            mmapper = vtk.vtkDataSetMapper()
            mmapper.SetInputData(part0)
            mmapper.SetScalarModeToUseCellData()
            mmapper.SetLookupTable(lut)
            mmapper.SetScalarRange(ca,cb)
            mactor = vtk.vtkActor()
            mactor.SetMapper(mmapper)
            mactor.GetProperty().EdgeVisibilityOn()
            mactor.GetProperty().SetColor(0, 0, 0)
            mactor.GetProperty().SetOpacity(1)
            mactor.SetVisibility(0)

            outline = vtk.vtkOutlineFilter()
            outline.SetInputConnection(reader.GetOutputPort())
            omapper = vtk.vtkPolyDataMapper()
            omapper.SetInputConnection(outline.GetOutputPort())
            oactor = vtk.vtkActor()
            oactor.SetMapper(omapper)
            oactor.GetProperty().SetColor(0, 0, 0)

            renderer.AddActor(mactor)
            renderer.AddActor(oactor)
        else:
            mapper.SetInputConnection(reader.GetOutputPort())
            mactor = None

        renderer.SetBackground(1, 1, 1) # Set background to white

        return renderer, mactor

    def __init__(self, parent, filename, mat_id=1):
        super(VTKViewer, self).__init__(parent=parent)

        self.setWindowTitle('VTKViewer: %s' % filename)

        vbox = QVBoxLayout()
        frame = QFrame()
        vbox.addWidget(frame)

        vtkWidget = QVTKRenderWindowInteractor(frame)
        vbox.addWidget(vtkWidget)

        ren, self.obj = self.view_vtk(filename, mat_id=mat_id)

        ren_win = vtkWidget.GetRenderWindow()
        ren_win.AddRenderer(ren)
        self.iren = ren_win.GetInteractor()
        vtkWidget.Initialize()
        vtkWidget.Start()

        if self.obj is not None:
            self.chbox = QCheckBox('matrix visibility')
            self.chbox.setChecked(False)
            self.chbox.stateChanged.connect(self.change_chbox_value)
            vbox.addWidget(self.chbox)
            self.slider = QSlider(Qt.Horizontal)
            # self.slider.setMinimum(0)
            # self.slider.setMaximum(1)
            self.slider.setValue(self.slider.maximum())
            self.slider.valueChanged.connect(self.change_slider_value)
            vbox.addWidget(self.slider)
            self.slider.setEnabled(False)

        vbox.addItem(QSpacerItem(0, 40, QSizePolicy.Minimum, QSizePolicy.Expanding))
        btn_quit = QPushButton('Quit Viewer', self)
        btn_quit.clicked.connect(self.close)

        vbox.addWidget(btn_quit)
        self.setLayout(vbox)

        self.toshot = self

    def change_chbox_value(self):
        if self.chbox.isChecked():
            self.obj.SetVisibility(True)
            self.slider.setEnabled(True)
        else:
            self.obj.SetVisibility(False)
            self.slider.setEnabled(False)
        self.iren.Render()

    def change_slider_value(self):
        val = (self.slider.value() - self.slider.minimum())\
            / float(self.slider.maximum())
        self.obj.GetProperty().SetOpacity(val)
        self.iren.Render()


def main():
    if len(sys.argv) > 1:
        app = QApplication(sys.argv)
        viewer = VTKViewer(None, sys.argv[1], mat_id=None)
        viewer.show()
        sys.exit(app.exec_())
    else:
        print("Usage: python vtk_viewer <filename.vtk>")


if __name__ == "__main__":
    main()