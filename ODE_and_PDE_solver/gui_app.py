#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
gui_app.py
Interfaz PyQt5 para ejecutar latex_solver.py y visualizar resultados.
Asume que latex_solver.py está en el mismo directorio que este archivo.
"""

import sys
import os
import subprocess
import re
from PyQt5 import QtWidgets, QtCore
from PyQt5.QtWidgets import QFileDialog, QMessageBox
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
import pandas as pd
import numpy as np
from matplotlib.animation import FuncAnimation

# Nombre del script solver (se buscará en el mismo directorio que gui_app.py)
SCRIPT_NAME = "latex_solver.py"


class SolverThread(QtCore.QThread):
    """
    Hilo que ejecuta el solver (latex_solver.py) como subprocess y emite stdout en tiempo real.
    """
    line_emitted = QtCore.pyqtSignal(str)
    finished_files = QtCore.pyqtSignal(list)

    def __init__(self, latex_str, script_path, python_exec=None):
        super().__init__()
        self.latex_str = latex_str
        self.script_path = script_path
        self.python_exec = python_exec if python_exec is not None else sys.executable

    def run(self):
        # Construir la lista de argumentos sin usar shell para evitar problemas de quoting
        args = [self.python_exec, self.script_path, "--latex", self.latex_str]
        try:
            proc = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
        except Exception as e:
            self.line_emitted.emit(f"ERROR launching solver: {e}")
            return

        written_files = []
        # leer stdout en tiempo real
        while True:
            line = proc.stdout.readline()
            if line == '' and proc.poll() is not None:
                break
            if line:
                line = line.rstrip()
                self.line_emitted.emit(line)
                # detectar patrones "Wrote <filename>.csv" o "Wrote <name>"
                m = re.search(r'Wrote\s+([^\s,]+\.csv)', line)
                if m:
                    written_files.append(m.group(1))
        proc.wait()
        self.finished_files.emit(written_files)


class MplCanvas(FigureCanvas):
    def __init__(self, parent=None, width=6, height=4, dpi=100):
        fig = Figure(figsize=(width, height), dpi=dpi)
        self.ax = fig.add_subplot(111)
        super().__init__(fig)


class MainWindow(QtWidgets.QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("LaTeX → Solver GUI")
        self.resize(1000, 700)

        # Widgets
        central = QtWidgets.QWidget()
        self.setCentralWidget(central)
        layout = QtWidgets.QVBoxLayout(central)

        top_h = QtWidgets.QHBoxLayout()
        layout.addLayout(top_h)

        # LaTeX input
        self.latex_edit = QtWidgets.QPlainTextEdit()
        self.latex_edit.setPlaceholderText(r'Introduce la ecuación en LaTeX, p.ej. \frac{d^2 y}{dt^2} + y = 0')
        top_h.addWidget(self.latex_edit, 2)

        right_v = QtWidgets.QVBoxLayout()
        top_h.addLayout(right_v, 1)

        # Buttons
        self.run_btn = QtWidgets.QPushButton("Run solver")
        self.run_btn.clicked.connect(self.on_run)
        right_v.addWidget(self.run_btn)

        self.open_btn = QtWidgets.QPushButton("Open CSV...")
        self.open_btn.clicked.connect(self.on_open_csv)
        right_v.addWidget(self.open_btn)

        self.clear_btn = QtWidgets.QPushButton("Clear log")
        self.clear_btn.clicked.connect(self.on_clear_log)
        right_v.addWidget(self.clear_btn)

        right_v.addStretch()

        # Status / log
        self.log = QtWidgets.QPlainTextEdit()
        self.log.setReadOnly(True)
        layout.addWidget(self.log, 1)

        # Plot canvas
        self.canvas = MplCanvas(self, width=8, height=4, dpi=100)
        layout.addWidget(self.canvas, 3)

        # Animation handle
        self.ani = None

        # default python exec (mismo intérprete que ejecuta la GUI)
        self.python_exec = sys.executable

        # Resolver path: script en mismo directorio que este archivo
        script_dir = os.path.dirname(os.path.abspath(__file__))  # directorio donde vive gui_app.py
        self.script_path = os.path.join(script_dir, SCRIPT_NAME)
        if not os.path.exists(self.script_path):
            self.append_log(f"ADVERTENCIA: no se encontró {SCRIPT_NAME} en {script_dir}")

    def append_log(self, s):
        self.log.appendPlainText(s)

    def on_clear_log(self):
        self.log.clear()

    def on_run(self):
        latex = self.latex_edit.toPlainText().strip()
        if latex == "":
            QMessageBox.warning(self, "Input required", "Introduce la ecuación en LaTeX antes de ejecutar.")
            return
        # comprobar que el script existe
        if not os.path.exists(self.script_path):
            QMessageBox.critical(self, "Script missing", f"No se encontró {self.script_path}. Coloca {SCRIPT_NAME} en la misma carpeta.")
            return

        self.append_log(f"Launching solver with LaTeX: {latex}")
        self.run_btn.setEnabled(False)
        self.thread = SolverThread(latex, self.script_path, python_exec=self.python_exec)
        self.thread.line_emitted.connect(self.append_log)
        self.thread.finished_files.connect(self.on_solver_finished)
        self.thread.start()

    def on_solver_finished(self, files):
        self.append_log("Solver finished. Files written: " + ", ".join(files))
        self.run_btn.setEnabled(True)
        # si hay archivos csv, abrir el primero automáticamente
        for f in files:
            # si f es nombre relativo, convertir a ruta absoluta (está en cwd del proceso: normalmente script_dir)
            potential = f
            if not os.path.isabs(potential):
                potential = os.path.join(os.path.dirname(self.script_path), potential)
            if os.path.exists(potential) and potential.endswith(".csv"):
                self.load_and_animate_csv(potential)
                break

    def on_open_csv(self):
        path, _ = QFileDialog.getOpenFileName(self, "Open CSV", ".", "CSV Files (*.csv)")
        if path:
            self.load_and_animate_csv(path)

    def load_and_animate_csv(self, path):
        self.append_log(f"Loading CSV: {path}")
        try:
            df = pd.read_csv(path)
        except Exception as e:
            QMessageBox.critical(self, "Error reading CSV", str(e))
            return

        cols = df.columns.tolist()
        # heurística: si hay columna 't'
        if 't' in cols:
            times = df['t'].values
            data = df.drop(columns=['t']).values
            col_names = [c for c in cols if c != 't']
        else:
            data = df.values
            times = np.arange(data.shape[0])
            col_names = cols

        # reconstruir datos si vienen pares real/imag
        plot_data = None
        if any("_real_" in c for c in cols):
            # intenta reconstruir magnitud a partir de pares *_real_i, *_imag_i
            real_cols = sorted([c for c in cols if "_real_" in c])
            imag_cols = sorted([c for c in cols if "_imag_" in c])
            if len(real_cols) == len(imag_cols):
                n = len(real_cols)
                plot_data = np.zeros((data.shape[0], n))
                for i in range(n):
                    rcol = real_cols[i]
                    icol = imag_cols[i]
                    plot_data[:, i] = np.abs(df[rcol].values + 1j*df[icol].values)
        if plot_data is None:
            # datos reales directos
            plot_data = np.array(data, dtype=float)

        # prepare plotting
        self.canvas.ax.clear()
        x = np.arange(plot_data.shape[1])
        self.line, = self.canvas.ax.plot(x, plot_data[0, :])
        self.canvas.ax.set_xlabel("space index")
        self.canvas.ax.set_ylabel("value")
        ymin = np.nanmin(plot_data); ymax = np.nanmax(plot_data)
        if ymin == ymax:
            ymin -= 1.0; ymax += 1.0
        self.canvas.ax.set_ylim(ymin, ymax)
        self.canvas.draw()

        # create animation
        if self.ani:
            self.ani.event_source.stop()
            self.ani = None

        def update(i):
            self.line.set_ydata(plot_data[i, :])
            return self.line,

        self.ani = FuncAnimation(self.canvas.figure, update, frames=plot_data.shape[0], interval=50, blit=True)
        self.append_log("Animation started.")

def main():
    app = QtWidgets.QApplication(sys.argv)
    w = MainWindow()
    w.show()
    sys.exit(app.exec_())


if __name__ == "__main__":
    main()
