# src/tanglepack_gui/main_window.py
from PySide6.QtWidgets import (
    QMainWindow,
    QDockWidget,
    QWidget,
    QVBoxLayout,
    QLineEdit,
    QPushButton,
    QLabel,
    QCheckBox,
)
from PySide6.QtGui import QAction
from PySide6.QtCore import Qt
from PySide6.QtWidgets import QMainWindow, QDockWidget, QMessageBox
from .adapters.map_parser import parse_map_text
from tanglepack_gui.theme import LIGHT, DARK
import numpy as np
import pyqtgraph as pg
from sympy import symbols, sympify, lambdify, Tuple
from .view.canvas import Canvas
from .adapters.workbench_view import manifold_arrays_for_fp

x, y = symbols("x y")


class MainWindow(QMainWindow):
    def __init__(self, workbench_cls):
        super().__init__()
        self.setWindowTitle("Tangle Workbench GUI")
        self.canvas = Canvas(self)
        self.setCentralWidget(self.canvas)

        # simple controls (you'll likely replace with a proper ControlPanel)
        dock = QDockWidget("Controls", self)
        panel = QWidget()
        dock.setWidget(panel)
        self.addDockWidget(Qt.LeftDockWidgetArea, dock)
        # self.addDockWidget(QMainWindow.DockWidgetArea.LeftDockWidgetArea, dock)

        v = QVBoxLayout(panel)
        self.f_edit = QLineEdit("y - 10 + x**2, -x")
        self.inv_edit = QLineEdit("-y, x + 10 - y**2")
        self.x0_edit = QLineEdit("4, -4")
        btn_build = QPushButton("Build Map")
        btn_fp = QPushButton("Find Fixed Point")
        btn_init = QPushButton("Init Both Manifolds")
        btn_grow = QPushButton("Grow ×1")
        btn_ints = QPushButton("Compute Intersections")
        btn_cut = QPushButton("Create Bridges")

        self.fp_label = QLabel("FP: —")
        v.addWidget(self.fp_label)

        # --- Orientation controls ---
        self.orient_label = QLabel("Orientation (approx dirs):")
        self.u_edit = QLineEdit("-1, 0")  # unstable approx dir (ux, uy)
        self.s_edit = QLineEdit("0, 1")  # stable approx dir   (sx, sy)
        self.auto_orient_chk = QCheckBox("Auto-orient before init")
        self.auto_orient_chk.setChecked(True)
        self.apply_orient_btn = QPushButton("Apply Orientation Now")

        v.addWidget(self.orient_label)
        v.addWidget(QLabel("Unstable dir (ux, uy):"))
        v.addWidget(self.u_edit)
        v.addWidget(QLabel("Stable dir (sx, sy):"))
        v.addWidget(self.s_edit)
        v.addWidget(self.auto_orient_chk)
        v.addWidget(self.apply_orient_btn)

        self.apply_orient_btn.clicked.connect(self.on_apply_orientation)

        for w in (
            self.f_edit,
            self.inv_edit,
            self.x0_edit,
            btn_build,
            btn_fp,
            btn_init,
            btn_grow,
            btn_ints,
            btn_cut,
        ):
            v.addWidget(w)

        self.workbench_cls = workbench_cls
        self.workbench = None
        self.fp = None

        btn_build.clicked.connect(self.on_build)
        btn_fp.clicked.connect(self.on_find_fp)
        btn_init.clicked.connect(self.on_init_both)
        btn_grow.clicked.connect(self.on_grow_once)
        btn_ints.clicked.connect(self.on_intersections)
        btn_cut.clicked.connect(self.on_cut)

        self._build_menus()

        self.canvas.clicked.connect(self.on_canvas_click)

    def _build_menus(self):
        menubar = self.menuBar()
        view_menu = menubar.addMenu("&View")

        act_light = QAction("Light mode", self)
        act_light.setCheckable(True)
        act_dark = QAction("Dark (blackboard)", self)
        act_dark.setCheckable(True)

        # exclusive behavior
        def select_light():
            act_dark.setChecked(False)
            act_light.setChecked(True)
            self.canvas.set_theme(LIGHT)

        def select_dark():
            act_light.setChecked(False)
            act_dark.setChecked(True)
            self.canvas.set_theme(DARK)

        act_light.triggered.connect(select_light)
        act_dark.triggered.connect(select_dark)

        # default: start in light (or change to dark if you prefer)
        act_light.setChecked(True)

        view_menu.addAction(act_light)
        view_menu.addAction(act_dark)

    def _parse_vec2(self, text: str) -> np.ndarray:
        try:
            a, b = [float(s.strip()) for s in text.split(",")]
            v = np.array([a, b], dtype=float)
            if not np.isfinite(v).all():
                raise ValueError("Non-finite value.")
            # normalize gently; direction only matters
            n = np.linalg.norm(v)
            return v if n == 0 else v / n
        except Exception:
            raise ValueError(f"Expected 'a, b' (two numbers). Got: {text!r}")

    def on_apply_orientation(self):
        if self.workbench is None or self.fp is None:
            QMessageBox.warning(
                self, "Not ready", "Build the map and find a fixed point first."
            )
            return
        try:
            u_dir = self._parse_vec2(self.u_edit.text())
            s_dir = self._parse_vec2(self.s_edit.text())
            approx_dirs = {"unstable": u_dir, "stable": s_dir}
            self.workbench.orient_eigenvectors(self.fp, approx_dirs)
            self.statusBar().showMessage("Eigenvectors oriented.", 3000)
        except Exception as e:
            QMessageBox.critical(self, "Orientation error", str(e))

    # @staticmethod
    # def _split_vector_expr(txt: str):
    #     """
    #     Accepts 'fx, fy' or '(fx, fy)' or '[fx, fy]'.
    #     Returns a pair of SymPy expressions (fx_expr, fy_expr).
    #     Raises ValueError with a helpful message on parse issues.
    #     """
    #     txt = txt.strip()
    #     if not txt:
    #         raise ValueError("Empty function text.")

    #     # Parenthesize to nudge sympy toward a vector; still handle python tuple fallback.
    #     expr = sympify(f"({txt})", convert_xor=True)

    #     # Normalize to a 2-tuple of sympy expressions
    #     if isinstance(expr, Tuple):  # SymPy Tuple
    #         parts = tuple(expr)
    #     elif isinstance(expr, (tuple, list)):  # Python tuple/list
    #         parts = tuple(expr)
    #     else:
    #         # Single expression given — user forgot the comma
    #         raise ValueError(
    #             "Expected two expressions separated by a comma: e.g. 'x + y**2, 0.3*x - 0.2*y'."
    #         )

    #     if len(parts) != 2:
    #         raise ValueError(f"Expected exactly 2 expressions, got {len(parts)}.")

    #     return parts[0], parts[1]

    # # def parse_map(self, txt):
    # #     fx, fy = sympify(f"({txt})").as_tuple()
    # #     return lambdify((x, y), (fx, fy), "numpy")

    # def parse_map(self, txt: str):
    #     fx_expr, fy_expr = self._split_vector_expr(txt)
    #     # Optional: sanity check for free symbols
    #     free = fx_expr.free_symbols | fy_expr.free_symbols
    #     allowed = {x, y}
    #     extra = free - allowed
    #     if extra:
    #         names = ", ".join(sorted(str(s) for s in extra))
    #         raise ValueError(
    #             f"Unknown symbol(s): {names}. Only 'x' and 'y' are allowed."
    #         )

    #     # Numpy backend is fine; SymPy will map standard funcs properly.
    #     return lambdify((x, y), (fx_expr, fy_expr), modules="numpy")

    # # def on_build(self):
    # #     from tanglepack.TangleWorkbench import TangleWorkbench as WB

    # #     f = self.parse_map(self.f_edit.text())
    # #     finv = self.parse_map(self.inv_edit.text())
    # #     self.workbench = WB(f, finv)  # DynamicalSystem created inside WB
    # #     # ready to go  (WB.__init__)  ⟶  DynamicalSystem(f, finv)
    # #     #                                + FixedPointSolver + ManifoldMachine + Tangle
    # #     # (all inside TangleWorkbench)
    # #     # :contentReference[oaicite:15]{index=15}

    # def on_build(self):
    #     try:
    #         f = self.parse_map(self.f_edit.text())
    #         finv = self.parse_map(self.inv_edit.text())
    #     except Exception as e:
    #         QMessageBox.critical(self, "Parse error", str(e))
    #         return

    #     try:
    #         from tanglepack.TangleWorkbench import TangleWorkbench

    #         self.workbench = TangleWorkbench(f, finv)
    #         # Optionally update a status bar:
    #         self.statusBar().showMessage("System built.", 3000)
    #     except Exception as e:
    #         QMessageBox.critical(
    #             self, "Workbench error", f"Failed to build system:\n{e}"
    #         )
    #         self.workbench = None

    def on_build(self):
        try:
            f = parse_map_text(self.f_edit.text())
            finv = parse_map_text(self.inv_edit.text())
        except Exception as e:
            QMessageBox.critical(self, "Parse error", str(e))
            return
        try:
            self.workbench = self.workbench_cls(f, finv)
            self.statusBar().showMessage("System built.", 3000)
        except Exception as e:
            QMessageBox.critical(
                self, "Workbench error", f"Failed to build system:\n{e}"
            )
            self.workbench = None

    # def on_find_fp(self):
    #     x0, y0 = [float(s.strip()) for s in self.x0_edit.text().split(",")]
    #     self.fp = self.workbench.construct_fixed_point(np.array([x0, y0], dtype=float))
    #     # persists in workbench.fixed_points; returns FixedPoint  :contentReference[oaicite:16]{index=16}

    def on_find_fp(self):
        if self.workbench is None:
            QMessageBox.warning(self, "Not ready", "Build the map first.")
            return

        x0, y0 = [float(s.strip()) for s in self.x0_edit.text().split(",")]
        self.fp = self.workbench.construct_fixed_point(np.array([x0, y0], dtype=float))

        # Extract the point for period-1 (or first point if higher period)
        coords = np.asarray(self.fp.coordinates)
        if coords.ndim == 2 and coords.shape[1] == 2:  # (period, 2)
            pt = coords[0]
        elif coords.ndim == 1 and coords.shape[0] == 2:  # (2,)
            pt = coords
        else:
            # Fallback: flatten to 2
            pt = coords.reshape(-1)[:2]

        # 1) put it in the label
        self.fp_label.setText(f"FP: ({pt[0]:.6f}, {pt[1]:.6f})")

        # 2) draw the point on the canvas (with hover tooltip)
        self.canvas.draw_fixed_point(pt)

        # Optional: center the view a bit around FP
        # (uncomment if you want)
        x, y = float(pt[0]), float(pt[1])
        self.canvas.setXRange(x - 5, x + 5, padding=0)
        self.canvas.setYRange(y - 5, y + 5, padding=0)

    # def on_init_both(self):
    #     self.workbench.initialize_both_manifolds(
    #         self.fp
    #     )  # fills workbench.manifolds  :contentReference[oaicite:17]{index=17}
    #     self._draw_all_manifolds()

    def on_init_both(self):
        if self.workbench is None or self.fp is None:
            QMessageBox.warning(self, "Not ready", "Find a fixed point first.")
            return

        if self.auto_orient_chk.isChecked():
            try:
                u_dir = self._parse_vec2(self.u_edit.text())
                s_dir = self._parse_vec2(self.s_edit.text())
                self.workbench.orient_eigenvectors(
                    self.fp, {"unstable": u_dir, "stable": s_dir}
                )
            except Exception as e:
                QMessageBox.critical(self, "Orientation error", str(e))
                return

        self.workbench.initialize_both_manifolds(self.fp)
        self._ensure_branch_indices()
        self._draw_all_manifolds()

    # def _draw_all_manifolds(self):
    #     # draw each manifold as one polyline; we’ll make them clickable by segment IDs later
    #     for key, M in self.workbench.manifolds.items():
    #         arr = M.get_point_array()
    #         # No per-segment IDs yet; this is the “coarse” view
    #         self.canvas.plot(arr[:, 0], arr[:, 1], pen=None)  # quick sketch

    def _ensure_branch_indices(self):
        if self.workbench is None or self.fp is None:
            return
        for (kfp, kstab, oi, bi), M in list(self.workbench.manifolds.items()):
            if kfp is self.fp and getattr(M, "branch_index", None) is None:
                M.branch_index = bi

    def _draw_all_manifolds(self):
        self.canvas.clear()

        if self.workbench is None or self.fp is None:
            return

        if self.fp is not None:
            coords = np.asarray(self.fp.coordinates)
            pt = (
                coords[0]
                if (coords.ndim == 2 and coords.shape[1] == 2)
                else coords.reshape(-1)[:2]
            )
            self.canvas.draw_fixed_point(pt)

        # iterate manifolds via the helper (now branch-aware)
        from .adapters.workbench_view import manifold_arrays_for_fp

        for (kstab, oi, bi), arr in manifold_arrays_for_fp(self.workbench, self.fp):
            # choose a pen per stability (explicit colors survive theme toggles)
            pen = (
                pg.mkPen("b", width=2)
                if kstab == "unstable"
                else pg.mkPen("r", width=2)
            )
            self.canvas.plot(arr[:, 0], arr[:, 1], pen=pen)

    def on_grow_once(self):
        self.workbench.grow_n_times(
            self.fp, "unstable", 1
        )  # or "stable"  :contentReference[oaicite:18]{index=18}
        self._draw_all_manifolds()

    def on_intersections(self):
        pts = np.array(
            list(self.workbench.compute_intersections(self.fp))
        )  # :contentReference[oaicite:19]{index=19}
        if pts.size:
            self.canvas.scatter_intersections(pts)

    def on_cut(self):
        bridges = self.workbench.create_bridges(
            self.fp
        )  # returns Bridge objects  :contentReference[oaicite:20]{index=20}
        # You can draw bridges as separate polylines via Bridge.plot(...) or traverse nodes.

    def on_canvas_click(self, x, y):
        # Fine picking by segment ID (reuses Tangle’s R-tree & lookup):
        # 1) coarse candidates (data-space bbox)  2) pixel-space refine
        vb = self.canvas.getViewBox()
        sp = vb.mapViewToScene(pg.Point(x, y))
        # Candidates from the R-tree
        tree = (
            self.workbench.Tangle._rtree
        )  # spatial index  :contentReference[oaicite:21]{index=21}
        lu = (
            self.workbench.Tangle._seg_lookup
        )  # id -> _Segment  :contentReference[oaicite:22]{index=22}
        sids = list(tree.intersection((x, y, x, y)))
        pick, dmin = None, 9e9
        for sid in sids:
            seg = lu[sid]
            P = np.array([seg.p0.get_point(), seg.p0_seg1.get_point()])
            # pixel-space distance to line segment
            A = np.array(vb.mapViewToScene(pg.Point(*P[0])).toTuple())
            B = np.array(vb.mapViewToScene(pg.Point(*P[1])).toTuple())
            S = np.array([sp.x(), sp.y()])
            AB = B - A
            AP = S - A
            t = np.clip(AP.dot(AB) / AB.dot(AB), 0, 1)
            proj = A + t * AB
            d = np.linalg.norm(S - proj)
            if d < dmin:
                dmin, pick = d, sid
        if pick is not None and dmin < 6:
            # toggle selection, then restyle
            sel = set(self.canvas.selected)
            if pick in sel:
                sel.remove(pick)
            else:
                sel.add(pick)
            self.canvas.style_selection(sel)
