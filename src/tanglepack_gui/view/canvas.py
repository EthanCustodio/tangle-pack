# src/tanglepack_gui/view/canvas.py
from tanglepack_gui.theme import PlotTheme, LIGHT, DARK, apply_pyqtgraph_theme
from PySide6.QtWidgets import QMainWindow, QToolTip
from PySide6.QtGui import QCursor
from PySide6.QtCore import Signal
import pyqtgraph as pg
import numpy as np


class Canvas(pg.PlotWidget):
    clicked = Signal(float, float)

    def __init__(self, parent=None):
        super().__init__(parent)
        self._theme: PlotTheme = LIGHT  # default
        apply_pyqtgraph_theme(self._theme)
        self.setBackground("w")

        self.showGrid(x=True, y=True, alpha=0.25)
        self.items = {}
        self.selected = set()

        # --- fixed point state ---
        self._fp_item = None
        self._fp_xy = None  # np.array([x, y])

        self.scene().sigMouseClicked.connect(self._on_click)
        self.scene().sigMouseMoved.connect(self._on_move)

    def set_theme(self, theme: PlotTheme):
        """Apply theme to canvas and refresh visuals."""
        self._theme = theme
        apply_pyqtgraph_theme(theme)  # axes, ticks, grid colors
        self.setBackground(theme.background)
        # re-enable grid with desired alpha since foreground changed
        self.showGrid(x=True, y=True, alpha=0.25)

        # restyle fixed point if present
        if self._fp_item is not None:
            self._fp_item.setBrush(pg.mkBrush(theme.fp_brush))
            self._fp_item.setPen(pg.mkPen(theme.fp_pen, width=1))
        self.update()

    def current_theme(self) -> PlotTheme:
        return self._theme

    def _on_click(self, ev):
        if ev.double():
            return
        vb = self.getViewBox()
        if not vb.sceneBoundingRect().contains(ev.scenePos()):
            return
        p = vb.mapSceneToView(ev.scenePos())
        self.clicked.emit(p.x(), p.y())

    def _on_move(self, pos):
        vb = self.getViewBox()
        if not vb.sceneBoundingRect().contains(pos):
            return
        p = vb.mapSceneToView(pos)

        # If we have an FP, check pixel distance to it
        if self._fp_xy is not None:
            fp_scene = vb.mapViewToScene(pg.Point(self._fp_xy[0], self._fp_xy[1]))
            dx = pos.x() - fp_scene.x()
            dy = pos.y() - fp_scene.y()
            dpx = (dx * dx + dy * dy) ** 0.5
            if dpx < 10:  # hover radius in pixels
                QToolTip.showText(
                    QCursor.pos(),
                    f"FP: ({self._fp_xy[0]:.6f}, {self._fp_xy[1]:.6f})",
                    self,
                )
            else:
                # hide tooltip when leaving the hover radius
                QToolTip.hideText()

    # def draw_fixed_point(self, xy, size=9):
    #     xy = np.asarray(xy, dtype=float)
    #     self._fp_xy = xy
    #     if self._fp_item is None:
    #         self._fp_item = pg.ScatterPlotItem(
    #             x=[xy[0]],
    #             y=[xy[1]],
    #             size=size,
    #             symbol="o",
    #             pen=pg.mkPen(width=1),
    #             brush=pg.mkBrush(0, 0, 0),  # black dot
    #         )
    #         self.addItem(self._fp_item)
    #     else:
    #         self._fp_item.setData([xy[0]], [xy[1]])

    def draw_fixed_point(self, xy, size=9):
        xy = np.asarray(xy, dtype=float)
        self._fp_xy = xy
        if self._fp_item is None:
            self._fp_item = pg.ScatterPlotItem(
                x=[xy[0]],
                y=[xy[1]],
                size=size,
                symbol="o",
                pen=pg.mkPen(self._theme.fp_pen, width=1),
                brush=pg.mkBrush(self._theme.fp_brush),
            )
            self.addItem(self._fp_item)
        else:
            self._fp_item.setData([xy[0]], [xy[1]])
            self._fp_item.setBrush(pg.mkBrush(self._theme.fp_brush))
            self._fp_item.setPen(pg.mkPen(self._theme.fp_pen, width=1))

    def upsert_segment(self, sid, xy, selected=False):
        item = self.items.get(sid)
        penw = 5 if selected else 2
        if item is None:
            item = self.plot(xy[:, 0], xy[:, 1], pen=pg.mkPen(width=penw))
            self.items[sid] = item
        else:
            item.setData(xy[:, 0], xy[:, 1])
            item.setPen(pg.mkPen(width=penw))

    def style_selection(self, selected_ids: set[int]):
        self.selected = set(selected_ids)
        for sid, it in self.items.items():
            it.setPen(pg.mkPen(width=5 if sid in self.selected else 2))

    def scatter_intersections(self, pts):
        if pts is None or len(pts) == 0:
            return
        self.plot(pts[:, 0], pts[:, 1], pen=None, symbol="o", symbolSize=6)
