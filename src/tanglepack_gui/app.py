# app.py — start the Qt app and show your Workbench GUI.
from PySide6.QtWidgets import QApplication
import sys

# If your MainWindow imports TangleWorkbench internally, no args needed:
from .main_window import MainWindow

from tanglepack.TangleWorkbench import TangleWorkbench


def main() -> None:
    app = QApplication(sys.argv)
    w = MainWindow(TangleWorkbench)  # uses your TangleWorkbench under the hood
    w.resize(1200, 900)
    w.show()
    sys.exit(app.exec())


if __name__ == "__main__":
    main()
