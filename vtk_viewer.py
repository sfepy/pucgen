import sys
import numpy as nm
import vtk
from PyQt5.QtWidgets import QVBoxLayout, QHBoxLayout, QSpacerItem,\
    QSizePolicy, QDialog, QPushButton, QCheckBox, QSlider, QApplication
from PyQt5.QtCore import Qt
from vtk.util import numpy_support
import vtk.qt
vtk.qt.QVTKRWIBase = 'QGLWidget'

from vtk.qt.QVTKRenderWindowInteractor import QVTKRenderWindowInteractor


class VTKViewer(QDialog):
    @staticmethod
    def view_vtk(filename, mat_id=1):
        reader = vtk.vtkUnstructuredGridReader()
        reader.SetFileName(filename)
        reader.Update()


        ca, cb = reader.GetOutput().GetCellData().GetScalars().GetRange()
        lut = vtk.vtkLookupTable()
        # lut.SetHueRange(0.667, 0)
        lut.SetAlphaRange(1, 1)
        lut.SetTableRange(ca, cb)
        mapper = vtk.vtkDataSetMapper()
        mapper.SetScalarModeToUseCellData()
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
            cell_data = reader.GetOutput().GetCellData()
            cell_data.SetActiveScalars('mat_id')
            mat_id_field = numpy_support.vtk_to_numpy(cell_data.GetScalars())

            ids = numpy_support.numpy_to_vtk(nm.where(mat_id_field == mat_id)[0])

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

        self.vbox = QVBoxLayout()
        vtkWidget = QVTKRenderWindowInteractor()
        vtkWidget.Initialize()
        vtkWidget.Start()
        ren, self.obj = self.view_vtk(filename, mat_id=mat_id)
        self.ren_win = vtkWidget.GetRenderWindow()
        self.ren_win.AddRenderer(ren)
        self.vbox.addWidget(vtkWidget)
        if self.obj is not None:
            self.chbox = QCheckBox('matrix visibility')
            self.chbox.setChecked(False)
            self.chbox.stateChanged.connect(self.change_chbox_value)
            self.vbox.addWidget(self.chbox)
            self.slider = QSlider(Qt.Horizontal)
            # self.slider.setMinimum(0)
            # self.slider.setMaximum(1)
            self.slider.setValue(self.slider.maximum())
            self.slider.valueChanged.connect(self.change_slider_value)
            self.vbox.addWidget(self.slider)
        self.vbox.addItem(QSpacerItem(0, 40, QSizePolicy.Minimum, QSizePolicy.Expanding))
        btn_quit = QPushButton('Quit Viewer', self)
        btn_quit.clicked.connect(self.close)
        self.vbox.addWidget(btn_quit)
        self.setLayout(self.vbox)

        self.toshot = self
        vtkWidget.show()

    def change_chbox_value(self):
        if self.chbox.isChecked():
            self.obj.SetVisibility(True)
        else:
            self.obj.SetVisibility(False)
        self.ren_win.GetInteractor().Render()

    def change_slider_value(self):
        val = (self.slider.value() - self.slider.minimum())\
            / float(self.slider.maximum())
        self.obj.GetProperty().SetOpacity(val)
        self.ren_win.GetInteractor().Render()


def main():
    if len(sys.argv) > 1:
        app = QApplication(sys.argv)
        viewer = VTKViewer(None, sys.argv[1], mat_id=1)
        # viewer = VTKViewer(None, sys.argv[1], mat_id=None)
        viewer.show()
        sys.exit(app.exec_())
    else:
        print("Usage: python vtk_viewer <filename.vtk>")


if __name__ == "__main__":
    main()