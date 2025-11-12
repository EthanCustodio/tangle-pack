// assets/bridge_selection.js
(function () {
  // front-end source of truth for selected bridge
  if (typeof window.__selected_bridge_idx === "undefined") window.__selected_bridge_idx = null;
  if (typeof window.__bridge_dirty === "undefined") window.__bridge_dirty = false;
  if (typeof window.__bridges_data === "undefined") window.__bridges_data = null;

  function innerPlot() {
    const wrap = document.getElementById("plot");
    return wrap && wrap.querySelector(".js-plotly-plot");
  }

  function findClosestBridge(clickX, clickY, bridges, maxDistance = 0.5) {
    if (!bridges || bridges.length === 0) return null;
    
    let minDist = Infinity;
    let closestIdx = null;
    
    bridges.forEach((bridge, idx) => {
      if (!bridge || !bridge.points || bridge.points.length === 0) return;
      
      // Calculate distance from click to all points on the bridge
      bridge.points.forEach(point => {
        const dx = point[0] - clickX;
        const dy = point[1] - clickY;
        const dist = Math.sqrt(dx * dx + dy * dy);
        
        if (dist < minDist) {
          minDist = dist;
          closestIdx = idx;
        }
      });
    });
    
    if (minDist <= maxDistance) {
      return { idx: closestIdx, distance: minDist };
    }
    return null;
  }

  function onClick(evt) {
    // ✅ gate by dropdown: only in "select_bridge"
    if (window.__click_mode !== "select_bridge") return;

    const gd = innerPlot();
    if (!gd || !gd._fullLayout) return;

    const full = gd._fullLayout;
    const rect = gd.getBoundingClientRect();
    const px = evt.clientX - rect.left - full.margin.l;
    const py = evt.clientY - rect.top  - full.margin.t;

    const w = (full._size && full._size.w) || (gd.clientWidth  - full.margin.l - full.margin.r);
    const h = (full._size && full._size.h) || (gd.clientHeight - full.margin.t - full.margin.b);
    if (px < 0 || py < 0 || px > w || py > h) return;

    const x = full.xaxis.p2c(px);
    const y = full.yaxis.p2c(py);

    // Find closest bridge
    const result = findClosestBridge(x, y, window.__bridges_data);
    
    if (result !== null) {
      window.__selected_bridge_idx = result.idx;
      window.__bridge_dirty = true;
      console.log(`[bridge select] Bridge ${result.idx + 1} selected (distance: ${result.distance.toFixed(4)})`);
    } else {
      window.__selected_bridge_idx = null;
      window.__bridge_dirty = true;
      console.log("[bridge select] No bridge near click");
    }
  }

  function attach() {
    const gd = innerPlot();
    if (!gd) { setTimeout(attach, 120); return; }

    gd.removeEventListener("click", onClick);
    gd.addEventListener("click", onClick);

    console.log("✓ bridge_selection.js attached");
  }

  if (document.readyState === "loading") {
    document.addEventListener("DOMContentLoaded", attach);
  } else {
    attach();
  }
})();