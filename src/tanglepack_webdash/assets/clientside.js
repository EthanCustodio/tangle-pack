// src/tanglepack_webdash/assets/clientside.js
window.dash_clientside = Object.assign({}, window.dash_clientside, {
  fp: {
    _innerPlotDiv: function() {
      const wrap = document.getElementById("plot");
      return wrap && wrap.querySelector(".js-plotly-plot");
    },
    _ensureTrace: function(gd, name) {
      const idx = (gd.data || []).findIndex(t => t && t.name === name);
      if (idx >= 0) return idx;
      Plotly.addTraces(gd, {
        name,
        x: [],
        y: [],
        mode: "markers",
        type: "scattergl",
        marker: { size: 8 },
        hovertemplate: "x=%{x:.3f}<br>y=%{y:.3f}<extra></extra>"
      });
      return gd.data.length - 1;
    },
    _fmt3: v => Number.isFinite(v)
      ? (Math.abs(v) >= 1e4 ? v.toExponential(3) : v.toFixed(3))
      : "--",

    click_to_data: function(n, evt, figure, mode, pts) {
      console.log("[click_to_data]", {n, mode, type: evt && evt.type, ptslen: (pts && pts.length) || 0});
      if (!evt || !figure) return pts || [];

      const gd = window.dash_clientside.fp._innerPlotDiv();
      if (!gd || !gd._fullLayout) return pts || [];

      const full = gd._fullLayout;
      const rect = gd.getBoundingClientRect();
      const px = evt.clientX - rect.left - full.margin.l;
      const py = evt.clientY - rect.top  - full.margin.t;

      const w = (full._size && full._size.w) || (gd.clientWidth  - full.margin.l - full.margin.r);
      const h = (full._size && full._size.h) || (gd.clientHeight - full.margin.t - full.margin.b);
      if (px < 0 || py < 0 || px > w || py > h) return pts || [];

      const x = full.xaxis.p2c(px);
      const y = full.yaxis.p2c(py);

      const out = Array.isArray(pts) ? pts.slice() : [];
      out.push([x, y]);
      return out;
    },

    // B) store -> text list "(x, y)" with 3 decimals
    points_to_label: function(points) {
      if (!Array.isArray(points) || points.length === 0) return "";
      const f = window.dash_clientside.fp._fmt3;
      return points.map(p => `(${f(p[0])}, ${f(p[1])})`).join("\n");
    },

    // C) store -> keep a visible marker layer in sync
    sync_points_trace: function(points) {
      const gd = window.dash_clientside.fp._innerPlotDiv();
      if (!gd) return "";
      const idx = window.dash_clientside.fp._ensureTrace(gd, "Clicked Points");
      const xs = (points || []).map(p => p[0]);
      const ys = (points || []).map(p => p[1]);
      Plotly.restyle(gd, { x: [xs], y: [ys] }, [idx]);
      return ""; // write to hidden sink
    },

    // D) mousemove -> live cursor readout
    move_to_label: function(n, evt, figure) {
      if (!evt || !figure || evt.type !== "mousemove") return window.dash_clientside.no_update;
      const gd = window.dash_clientside.fp._innerPlotDiv();
      if (!gd || !gd._fullLayout) return window.dash_clientside.no_update;
      const full = gd._fullLayout, rect = gd.getBoundingClientRect();
      const px = evt.clientX - rect.left - full.margin.l;
      const py = evt.clientY - rect.top  - full.margin.t;
      const w = (full._size && full._size.w) || (gd.clientWidth  - full.margin.l - full.margin.r);
      const h = (full._size && full._size.h) || (gd.clientHeight - full.margin.t - full.margin.b);
      if (px < 0 || py < 0 || px > w || py > h) return window.dash_clientside.no_update;
      const x = full.xaxis.p2c(px), y = full.yaxis.p2c(py);
      const f = window.dash_clientside.fp._fmt3;
      return `x: ${f(x)} | y: ${f(y)}`;
    },
    
    pull_points: function(n_intervals) {
      // Mirror the front-end array into the Dash Store
      return Array.isArray(window.__fp_points) ? window.__fp_points : [];
    },
    
    set_mode: function(mode) {
      window.__click_mode = mode; // e.g., "none" or "fp_guess" or "select_bridge"
      return ""; // write to hidden div
    },
    
    pull_points_dirty: function(n) {
      if (!Array.isArray(window.__fp_points)) return [];
      if (!window.__fp_dirty) return window.dash_clientside.no_update;
      window.__fp_dirty = false;
      return window.__fp_points;
    },
  },
  
  bridge: {
    // Pull bridge selection from front-end into Dash Store
    pull_selection: function(n) {
      if (typeof window.__selected_bridge_idx === "undefined") return null;
      if (!window.__bridge_dirty) return window.dash_clientside.no_update;
      window.__bridge_dirty = false;
      return window.__selected_bridge_idx;
    },
    
    // Update bridge data in front-end (called when bridges are created/updated)
    push_bridges_data: function(bridges_json) {
      // bridges_json should be an array of {points: [[x,y], ...], ...}
      window.__bridges_data = bridges_json;
      return "";
    },
    
    // Format selection info text
    selection_to_label: function(idx, distance) {
      if (idx === null || idx === undefined) {
        return "No bridge selected";
      }
      const dist_str = (typeof distance === "number") ? ` (distance: ${distance.toFixed(4)})` : "";
      return `Bridge ${idx + 1} selected${dist_str}`;
    },
  }
});