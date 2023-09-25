# 导入pyside2所必须的模块
from PySide2.QtWidgets import QApplication, QMainWindow, \
    QPushButton, QPlainTextEdit, QMessageBox, QLineEdit, QCheckBox, \
    QFileDialog, QDoubleSpinBox, QTableWidgetItem
from PySide2.QtWidgets import QApplication, QGraphicsView, QGraphicsScene
from pathlib import Path
from PySide2.QtUiTools import QUiLoader
from PySide2.QtCore import QFile, QTimer
import PySide2.QtCore as QtCore
# 从ansys.fluent.core导入pyfluent包
import ansys.fluent.core as pyfluent
#导入ansys.fluent.visualizatuon 可视化后处理包
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from ansys.fluent.visualization import set_config
#导入后处理包中的Graphics方法用来绘制云图
from ansys.fluent.visualization.pyvista import Graphics
set_config(blocking=True,set_view_on_display = 'isometric')
# 导入sys系统模块
import sys
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
#导入csv模块
import csv
#导入json模块
import json
import time

#定义登录窗口
class Login:
    def __init__(self):
        qfile_stats = QFile("login1.ui")
        qfile_stats.open(QFile.ReadOnly)
        qfile_stats.close()
        self.main_window = QMainWindow()
        self.ui = QUiLoader().load(qfile_stats)
        self.ui.setWindowFlags(QtCore.Qt.FramelessWindowHint)
        self.ui.setAttribute(QtCore.Qt.WA_TranslucentBackground)



#定义控制台输出的类ConsoleOutput
class ConsoleOutput:
    def __init__(self, text_edit):
        self.text_edit = text_edit

    def write(self, message):
        cursor = self.text_edit.textCursor()
        cursor.movePosition(cursor.End)
        cursor.insertText(message)
        self.text_edit.setTextCursor(cursor)
        self.text_edit.ensureCursorVisible()

    def flush(self):
        pass

# 定义操作界面UI窗口的Stats窗口
class Stats(QMainWindow):
    def __init__(self):
        super().__init__()
        qfile_stats = QFile(r"untitled2.ui")
        qfile_stats.open(QFile.ReadOnly)
        qfile_stats.close()
        self.main_window = QMainWindow()
        self.ui = QUiLoader().load(qfile_stats)
        self.ui.file1.clicked.connect(self.choose_mesh_file_path) #点击file1控件选择网格文件
        sys.stdout = ConsoleOutput(self.ui.o_put)
        sys.stderr = ConsoleOutput(self.ui.o_put)#控制台输出转移至UI文件之中
        self.ui.showgui.clicked.connect(self.showGui)
        self.ui.nogui.clicked.connect(self.showGui) #点击showgui或者nogui控件选择是否显示Fluent的GUI界面
        self.ui.file3.clicked.connect(self.choose_case_file_path) #点击file3控件选择保存case文件的保存路径
        self.ui.file2.clicked.connect(self.choose_result_file_path)#点击flie2控件选择保存result文件的保存路径
        self.ui.x_axis.valueChanged.connect(self.gravity)
        self.ui.y_axis.valueChanged.connect(self.gravity)
        self.ui.z_axis.valueChanged.connect(self.gravity) #设置重力大小
        self.ui.energy_equ.clicked.connect(self.energy) # 点击energy_equ控件打开能量方程
        self.ui.viscous.activated.connect(self.viscous_q) #点击viscous选择特定的能量方程
        self.ui.Id_boundry.clicked.connect(self.identificationBroundary) # 点击Id_boundary控件识别边界条件
        self.ui.csv_file.clicked.connect(self.read_csv_file) #点击csv——file控件读取边界条件网格
        self.ui.initialization.clicked.connect(self.initializationCal) #点击initiallization控件进行计算初始化
        self.ui.number.valueChanged.connect(self.sloveNum) #点击number控件进行求解步数的设定
        self.ui.calulate.clicked.connect(self.cal) #点击cal控件进行数值计算
        self.ui.methods.activated.connect(self.viscous_q)
        self.ui.show_mesh.clicked.connect(self.show_mesh) #点击showmesh控件显示网格
        self.ui.tunnel_in.clicked.connect(self.tunnel_in_rplot)
        self.ui.tunnel_out.clicked.connect(self.tunnel_out_rplot)
        self.ui.polution.clicked.connect(self.mid_tunnel_co)
        self.ui.velocity.clicked.connect(self.mid_tunnel_ve)
        self.ui.casefile.clicked.connect(self.case_setting)
        self.ui.meshfloder.clicked.connect(self.mesh_floder)
        self.ui.cal.clicked.connect(self.calulate)

    # 定义显示gui窗口的函数
    def showGui(self):
        global session
        if self.ui.showgui.isChecked() == True and self.ui.nogui.isChecked() == False:
            # 使用求解器模式打开
            session = pyfluent.launch_fluent(precision="double", processor_count=8, show_gui=True,
                                             mode="solver")
            session.tui.file.read_case(case_file_name="ISO.cas.h5")
            print( "标准CASE文件读取完成！！！")
        elif self.ui.showgui.isChecked() == False and self.ui.nogui.isChecked() == True:
            session = pyfluent.launch_fluent(precision="double", processor_count=8, show_gui=False,
                                             mode="solver")
            session.tui.file.read_case(case_file_name=r"ISO.cas.h5")
            print("标准CASE文件读取完成！！！" )

    #定义选择网格文件函数
    def choose_mesh_file_path(self):
        file_dialog = QFileDialog()
        file_dialog.setNameFilter("Mesh Files (*.msh)")
        file_dialog.setDefaultSuffix('msh')
        file_directory, _ = file_dialog.getOpenFileName(self, "选择.msh文件", "", "Mesh Files (*.msh)")
        if file_directory:
            self.ui.mesh_file_path.setText(file_directory)  # 返回文件地址路径
            session.tui.file.replace_mesh(case_file_name=file_directory)  # 替换网格数据
            session.tui.mesh.check()  # 检查网格
            session.tui.define.units('length', 'm')  # 定义单位长度
            session.tui.define.operating_conditions.operating_pressure('101325')  # 定义工作压力
            print("网格文件读取完成！！！")

    # 定义求解重力加速度方向以及位置
    def gravity(self):
        self.ui.x_axis.setMaximum(10.00)
        self.ui.x_axis.setMinimum(-10.00)
        self.ui.y_axis.setMaximum(10.00)
        self.ui.y_axis.setMinimum(-10.00)
        self.ui.z_axis.setMaximum(10.00)
        self.ui.z_axis.setMinimum(-10.00)
        self.ui.x_axis.setWrapping(True)
        self.ui.x_axis.setSingleStep(0.01)
        self.ui.y_axis.setWrapping(True)
        self.ui.y_axis.setSingleStep(0.01)
        self.ui.z_axis.setWrapping(True)
        self.ui.z_axis.setSingleStep(0.01)
        x_axis = str(self.ui.x_axis.value())
        y_axis = str(self.ui.y_axis.value())
        z_axis = str(self.ui.z_axis.value())
        axis = [x_axis,y_axis,z_axis]
        return axis

    #定义保存计算结果文件的函数
    def choose_result_file_path(self):
        file_dialog = QFileDialog()
        file_dialog.setNameFilter("DAT Files (*.dat.h5)")
        file_dialog.setDefaultSuffix('dat.h5')
        file_directory, _ = file_dialog.getSaveFileName(self, "保存文件", "", "dat Files (*.dat.h5)")
        if file_directory:
            session.tui.file.write_data(file_directory)
            print("结果文件保存完成！！！")
            self.ui.case_file_path.setText(file_directory)

    #定义保存case文件的函数
    def choose_case_file_path(self):
        file_dialog = QFileDialog()
        file_dialog.setNameFilter("Case Files (*.cas)")
        file_dialog.setDefaultSuffix('cas')
        file_directory, _ = file_dialog.getSaveFileName(self, "保存文件", "", "Case Files (*.cas)")
        if file_directory:
            session.tui.file.write_case(case_file_name=file_directory)
            print("case文件保存完成！！！")
            self.ui.case_file_path.setText(file_directory)



    #定义激活能量方程的函数
    def energy(self):
       if self.ui.energy_equ.isChecked() == True:
           session.tui.define.models.energy('yes', ',', ',', ',', ',')  # 激活能量方程
           axis = self.gravity()
           session.tui.define.operating_conditions.gravity('yes', axis[0],axis[1],axis[2])  # 定义坐标轴
           print("\033[91m" + "能量方程已激活！！！" + "\033[0m")

    #定义湍流方程的函数
    def viscous_q(self):
        if self.ui.viscous.currentText() == 'k-kl-w':
            session.tui.define.models.viscous.k_kl_w('yes')
            print("湍流方程已设置！！！")
        elif self.ui.viscous.currentText() == 'ke_realizable':
            session.tui.define.models.viscous.ke_realizable('yes')
            print("湍流方程已设置！！！")
        elif self.ui.viscous.currentText() == 'kw_bsl':
            session.tui.define.models.viscous.kw_bsl('yes')
            print("湍流方程已设置！！！")
        elif self.ui.viscous.currentText() == 'kw_geko':
            session.tui.define.models.viscous.kw_geko('yes')
            print("湍流方程已设置！！！")
        elif self.ui.viscous.currentText() == 'kw_sst':
            session.tui.define.models.viscous.kw_sst('yes')
            print("湍流方程已设置！！！")
        elif self.ui.viscous.currentText() == 'kw_standard':
            session.tui.define.models.viscous.kw_standard('yes')
            print("湍流方程已设置！！！")

    #识别流体边界条件
    def identificationBroundary(self):
        session.tui.define.boundary_conditions.list_zones()

    #定义设置求解步数的边界条件
    def boundaryNumber(self):
        self.ui.boundry_number.setMaximum(10.00)
        self.ui.boundry_number.setMinimum(1)
        self.ui.boundry_number.setWrapping(True)
        rownumber = self.ui.boundry_number.value()
        self.ui.bounary.setRowCount(rownumber)

    #定义读取边界条件表格的函数
    def read_csv_file(self):
        file_dialog = QFileDialog()
        file_dialog.setNameFilter("CSV Files (*.csv)")
        file_dialog.setDefaultSuffix('csv')
        file_directory, _ = file_dialog.getOpenFileName(self, "选择.csv文件", "", "CSV Files (*.csv)")
        if file_directory:
            self.ui.csv_file_path.setText(file_directory)
            with open(file_directory) as file:
                reader = csv.reader(file)
                header_row = next(reader)  # 读取表头行
                print(header_row)
                data = []  # 存储非空行数据
                for row in reader:
                    if not all(cell.strip() == '' for cell in row):
                        data.append(row)
            self.ui.bounary.setHorizontalHeaderLabels(header_row)
            self.ui.bounary.setRowCount(len(data))
            self.ui.bounary.setColumnCount(len(data[0]))
            for i, row in enumerate(data):
                for j, cell in enumerate(row):
                    self.ui.bounary.setItem(i, j, QTableWidgetItem(cell))
            for d in data:
                if d[0] == "velocity-inlet":
                    self.velocityInlet(d)
                elif d[0] == "pressure-inlet":
                    self.pressureInlet(d)
                elif d[0] == "pressure-outlet":
                    self.pressureOutlet(d)
                elif d[0] == "mass-flow-inlet":
                    self.massFlowInlet(d)

    #定义velocity——inlet的边界条件
    def velocityInlet(self,setting):
        session.tui.define.boundary_conditions.zone_type(setting[1], setting[0])
        session.setup.boundary_conditions.velocity_inlet[setting[1]].vmag = {'option': 'constant or expression',
                                                                           'constant': float(setting[2])}
        session.setup.boundary_conditions.velocity_inlet[setting[1]].ke_spec = 'Intensity and Hydraulic Diameter'
        session.setup.boundary_conditions.velocity_inlet[setting[1]].turb_hydraulic_diam = float(setting[3])
        session.setup.boundary_conditions.velocity_inlet[setting[1]].turb_intensity = float(setting[4])
        session.setup.boundary_conditions.velocity_inlet[setting[1]].species_in_mole_fractions = False
        session.setup.boundary_conditions.velocity_inlet[setting[1]].mf = json.loads(setting[5])
        print(session.setup.boundary_conditions.velocity_inlet[setting[1]].mf())
        session.setup.boundary_conditions.velocity_inlet[setting[1]].t = {'option': 'constant or expression',
                                                                        'constant':float(setting[6])}
        print(session.setup.boundary_conditions.velocity_inlet[setting[1]].t())
        print(f"{setting[1]}边界条件已设置完成")


    #定义pressure-inlet边界条件
    def pressureInlet(self, setting):
        session.tui.define.boundary_conditions.zone_type(setting[1], setting[0])
        session.setup.boundary_conditions.pressure_inlet[setting[1]].p0 = {'option': 'constant or expression',
                                                                             'constant': float(setting[2])}
        session.setup.boundary_conditions.pressure_inlet[setting[1]].ke_spec = 'Intensity and Hydraulic Diameter'
        session.setup.boundary_conditions.pressure_inlet[setting[1]].turb_hydraulic_diam =float (setting[3])
        session.setup.boundary_conditions.pressure_inlet[setting[1]].turb_intensity = float(setting[4])
        session.setup.boundary_conditions.pressure_inlet[setting[1]].species_in_mole_fractions = False
        session.setup.boundary_conditions.pressure_inlet[setting[1]].mf = json.loads(setting[5])
        print(session.setup.boundary_conditions.pressure_inlet[setting[1]].mf())
        session.setup.boundary_conditions.pressure_inlet[setting[1]].t0= {'option': 'constant or expression',
                                                                          'constant': float(setting[6])}
        print(session.setup.boundary_conditions.pressure_inlet[setting[1]].t0())
        if setting[7] == '隧道入口':
            session.tui.solve.report_definitions.add(
                "in",
                "surface-massflowrate",
                "surface-names",
                setting[1],
                "()",
                "quit",
            )
            session.tui.solve.report_files.add('in')
            session.tui.solve.report_files.edit('tunnel-in', "report-defs", ["in"],'file-name','tunnel-in.out')
            session.tui.solve.report_plots.add('tunnel-in')
            session.tui.solve.report_plots.edit('tunnel-in', "report-defs", ["in"])
        print(f"{setting[1]}边界条件已设置完成")

    #定义pressure-outlet的边界条件
    def pressureOutlet(self,setting):
        session.tui.define.boundary_conditions.zone_type(setting[1], setting[0])
        session.setup.boundary_conditions.pressure_outlet[setting[1]].p = {'option': 'constant or expression',
                                                                           'constant': float(setting[2])}
        session.setup.boundary_conditions.pressure_outlet[setting[1]].ke_spec = 'Intensity and Hydraulic Diameter'
        session.setup.boundary_conditions.pressure_outlet[setting[1]].turb_hydraulic_diam = float(setting[3])
        session.setup.boundary_conditions.pressure_outlet[setting[1]].turb_intensity = float(setting[4])
        session.setup.boundary_conditions.pressure_outlet[setting[1]].species_in_mole_fractions = False
        session.setup.boundary_conditions.pressure_outlet[setting[1]].mf = json.loads(setting[5])
        print(session.setup.boundary_conditions.pressure_outlet[setting[1]].species_in_mole_fractions())
        session.setup.boundary_conditions.pressure_outlet[setting[1]].t0 = {'option': 'constant or expression',
                                                                          'constant': float(setting[6])}
        print(session.setup.boundary_conditions.pressure_outlet[setting[1]].t0)
        if setting[7] == '隧道出口':
            session.tui.solve.report_definitions.add(
                "out",
                "surface-massflowrate",
                "surface-names",
                setting[1],
                "()",
                "quit",
            )
            session.tui.solve.report_files.add('out')
            session.tui.solve.report_files.edit('tunnel-out', "report-defs", ["out"],'file-name','tunnel-out.out')
            session.tui.solve.report_plots.add('tunnel-out')
            session.tui.solve.report_plots.edit('tunnel-out', "report-defs", ["out"])
        print(f"{setting[1]}边界条件已设置完成")

    #定义mass——flow-inlet的边界条件
    def massFlowInlet(self, setting):
        session.tui.define.boundary_conditions.zone_type(setting[1], setting[0])
        session.setup.boundary_conditions.mass_flow_inlet[setting[1]].mass_flow = {'option': 'constant or expression',
                                                                           'constant': float(setting[2])}
        session.setup.boundary_conditions.mass_flow_inlet[setting[1]].ke_spec = 'Intensity and Hydraulic Diameter'
        session.setup.boundary_conditions.mass_flow_inlet[setting[1]].turb_hydraulic_diam =float(setting[3])
        session.setup.boundary_conditions.mass_flow_inlet[setting[1]].turb_intensity = float(setting[4])
        session.setup.boundary_conditions.mass_flow_inlet[setting[1]].species_in_mole_fractions = False
        session.setup.boundary_conditions.mass_flow_inlet[setting[1]].mf = json.loads(setting[5])
        print(session.setup.boundary_conditions.mass_flow_inlet[setting[1]].mf())
        session.setup.boundary_conditions.mass_flow_inlet[setting[1]].t0 = {'option': 'constant or expression',
                                                                           'constant': float(setting[6])}
        print(session.setup.boundary_conditions.mass_flow_inlet[setting[1]].t0)
        print(f"{setting[1]}边界条件已设置完成")

    #定义初始化函数的边界条件
    def initializationCal(self):
        if self.ui.initialization.isChecked() == True:
            # 进行初始化
            session.tui.define.boundary_conditions.fluid(
                "fluid",
                "yes",
                "carbon-monoxide-air",
                "no",
                "no",
                "no",
                "no",
                "0",
                "no",
                "0",
                "no",
                "0",
                "no",
                "0",
                "no",
                "0",
                "no",
                "1",
                "no",
                "no",
                "no",
                "no",
                "no",
            )
            session.tui.solve.initialize.hyb_initialization()
            print("计算初始化完成")

    #设置求解方法
    def methods(self):
        if self.ui.methods.currentText() == 'SIMPLE':
            session.solution.methods.p_v_coupling.flow_scheme = 'SIl.MPLE'
            print("求解方法已设置！！！")
        elif self.ui.methods.currentText() == 'PISO':
            session.solution.methods.p_v_coupling.flow_scheme = 'PISO'
            print("求解方法已设置！！！")
        elif self.ui.methods.currentText() == 'SIMPLEC':
            session.solution.methods.p_v_coupling.flow_scheme = 'SIMPLEC'
            print("求解方法已设置！！！")
        elif self.ui.methods.currentText() == 'Coupled':
            session.solution.methods.p_v_coupling.flow_scheme = 'Coupled'
            print("求解方法已设置！！！")

    #定义求解步数的函数
    def sloveNum(self):
        global solve_num
        solve_num = [0]
        if self.ui.number.value() != solve_num[-1]:
            solve_num.append(self.ui.number.value())
        return solve_num[-1]

    #定义计算函数
    def cal(self):
        session.tui.solve.iterate(self.sloveNum())

        # 显示网格
    def show_mesh(self):
        global graphics
        graphics = Graphics(session=session)
        mesh1 = graphics.Meshes['mesh-1']
        mesh1.show_edges = True
        name = mesh1.surfaces_list.allowed_values
        mesh1.surfaces_list = name
        mesh1.display('window-0')


    def tunnel_in_rplot(self):
        x = []
        y = []
        with open(r"in.out", 'r') as f:
            for index, line in enumerate(f):
                if index >= 4:
                    columns = line.strip().split()
                    x.append(float(columns[0]))
                    y.append(float(columns[1]))

        fig = Figure()
        ax = fig.add_subplot(111)
        ax.plot(x, y)
        canvas = FigureCanvas(fig)
        canvas.draw()
        # Set the size of the canvas to match the size of the QGraphicsView
        canvas.setGeometry(self.ui.in_rplot.geometry())
        ax.set_position([0, 0, 1, 1])
        ax.tick_params(axis='both', which='major', pad=10, labelsize=6)
        ax.xaxis.labelpad = 10
        ax.yaxis.labelpad = 10
        fig.tight_layout()
        plt.show()
        # Clear the contents of the QGraphicsView before adding the canvas
        self.ui.in_rplot.setScene(QGraphicsScene())
        # Add the canvas to the QGraphicsView
        self.ui.in_rplot.setScene(QGraphicsScene(self.ui.in_rplot))
        self.ui.in_rplot.scene().addWidget(canvas)
        ax.figure.canvas.draw_idle()


    def tunnel_out_rplot(self):
        x = []
        y = []
        with open(r"out.out", 'r') as f:
            for index, line in enumerate(f):
                if index >= 4:
                    columns = line.strip().split()
                    x.append(float(columns[0]))
                    y.append(float(columns[1]))

        fig = Figure()
        ax = fig.add_subplot(111)
        ax.plot(x, y)
        canvas = FigureCanvas(fig)
        canvas.draw()
        # Set the size of the canvas to match the size of the QGraphicsView
        canvas.setGeometry(self.ui.out_rplot.geometry())
        ax.set_position([0, 0, 1, 1])
        ax.tick_params(axis='both', which='major', pad=10, labelsize=6)
        ax.xaxis.labelpad = 10
        ax.yaxis.labelpad = 10
        fig.tight_layout()
        plt.show()
        self.ui.in_rplot.setScene(QGraphicsScene())
        self.ui.out_rplot.setScene(QGraphicsScene(self.ui.out_rplot))
        self.ui.out_rplot.scene().addWidget(canvas)
        ax.figure.canvas.draw_idle()


    def mid_tunnel_co(self):
        graphics = Graphics(session=session)
        surf_plane = graphics.Surfaces["mid-tunnel"]
        surf_plane.definition.type = "plane-surface"
        plane_surface = surf_plane.definition.plane_surface
        if self.ui.planestyle.currentText() == 'zx-plane':
            surf_plane.definition.plane_surface.creation_method = 'zx-plane'
            plane_surface.zx_plane.y = float(self.ui.planevalue.value())
        elif self.ui.planestyle.currentText() == 'yz-plane':
            surf_plane.definition.plane_surface.creation_method = 'yz-plane'
            plane_surface.yz_plane.x = float(self.ui.planevalue.value())
        elif self.ui.planestyle.currentText() == 'xy-plane':
            surf_plane.definition.plane_surface.creation_method = 'xy-plane'
            plane_surface.xy_plane.z = float(self.ui.planevalue.value())
        co_contour = graphics.Contours["mid-co"]
        co_contour.field = "molef-co"
        co_contour.surfaces_list = ["mid-tunnel"]
        co_contour.display("window-2")

    def mid_tunnel_ve(self):
        graphics = Graphics(session=session)
        surf_plane = graphics.Surfaces["mid-tunnel"]
        surf_plane.definition.type = "plane-surface"
        plane_surface = surf_plane.definition.plane_surface
        if self.ui.planestyle.currentText() == 'zx-plane':
            surf_plane.definition.plane_surface.creation_method = 'zx-plane'
            plane_surface.zx_plane.y = float(self.ui.planevalue.value())
        elif self.ui.planestyle.currentText() == 'yz-plane':
            surf_plane.definition.plane_surface.creation_method = 'yz-plane'
            plane_surface.yz_plane.x = float(self.ui.planevalue.value())
        elif self.ui.planestyle.currentText() == 'xy-plane':
            surf_plane.definition.plane_surface.creation_method = 'xy-plane'
            plane_surface.xy_plane.z = float(self.ui.planevalue.value())
        co_contour = graphics.Contours["mid-ve"]
        co_contour.field = "velocity-magnitude"
        co_contour.surfaces_list = ["mid-tunnel"]
        co_contour.display("window-3")

    def case_setting(self):
        file_dialog = QFileDialog()
        file_dialog.setNameFilter("Cas Files (*.cas.h5)")
        file_dialog.setDefaultSuffix('cas.h5')
        file_directory, _ = file_dialog.getOpenFileName(self, "选择.cas文件", "", "Cas Files (*.cas.h5)")
        if file_directory:
            self.ui.case_setting_path.setText(file_directory)
            try:
                global solver_session
                solver_session = pyfluent.launch_fluent(precision="double", processor_count=8, show_gui=True,
                                                 mode="solver")
                solver_session.tui.file.read_case(case_file_name=file_directory)  # 定义工作压力
                print("case文件设置读取完成！！！")
                self.ui.case_setting_path.setText(file_directory)
            except NameError:
                print('未打开fluent，现在打开fluent GUI界面')
                solver_session = pyfluent.launch_fluent(precision="double", processor_count=8, show_gui=True,
                                                 mode="solver")
                solver_session.tui.file.read_case(case_file_name=file_directory)
                self.ui.case_setting_path.setText(file_directory)

    def mesh_floder(self):
        file_dialog = QFileDialog()
        file_dialog.setFileMode(QFileDialog.Directory)
        file_directory = file_dialog.getExistingDirectory(self, "选择文件夹", "")
        if file_directory:
            self.ui.mesh_floder.setText(file_directory)
            return file_directory

    def calulate(self):
        file_floder = Path(self.mesh_floder())
        mesh_flie_list = file_floder.glob("*.msh*")
        for mesh in mesh_flie_list:
            solver_session.tui.file.replace_mesh(case_file_name=mesh)
            solver_session.tui.solve.iterate(100)
            solver_session.tui.file.write_case(case_file_name=mesh)
        self.ui.Bar1.value = 10
        time.sleep(30)
        self.ui.Bar1.value = 80
    # 打印所选文件夹路径

if __name__ == '__main__':
    app = QApplication(sys.argv)
    login = Login()
    stat = Stats()
    login.ui.show()
    login.ui.enter.clicked.connect(stat.ui.show)
    login.ui.enter.clicked.connect(login.ui.hide)


    app.exec_()