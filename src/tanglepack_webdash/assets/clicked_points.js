// assets/clicked_points.js
(function () {
  const CLEAR_BTN_ID = "btn-clear-guesses";

  // Front-end source of truth
  if (!Array.isArray(window.__fp_points)) window.__fp_points = [];
  if (typeof window.__fp_dirty === "undefined") window.__fp_dirty = false;
  if (typeof window.__click_mode === "undefined") window.__click_mode = "none";

  function innerPlot() {
    const wrap = document.getElementById("plot");
    return wrap && wrap.querySelector(".js-plotly-plot");
  }

  function onClick(evt) {
    // Only capture clicks in "fp_guess" mode
    if (window.__click_mode !== "fp_guess") {
      console.log("[fp click] ignored - mode is:", window.__click_mode);
      return;
    }

    const gd = innerPlot();
    if (!gd || !gd._fullLayout) {
      console.log("[fp click] no plot or layout");
      return;
    }

    const full = gd._fullLayout;
    const rect = gd.getBoundingClientRect();
    const px = evt.clientX - rect.left - full.margin.l;
    const py = evt.clientY - rect.top - full.margin.t;

    const w = (full._size && full._size.w) || (gd.clientWidth - full.margin.l - full.margin.r);
    const h = (full._size && full._size.h) || (gd.clientHeight - full.margin.t - full.margin.b);
    
    if (px < 0 || py < 0 || px > w || py > h) {
      console.log("[fp click] outside plot area");
      return;
    }

    const x = full.xaxis.p2c(px);
    const y = full.yaxis.p2c(py);

    window.__fp_points.push([x, y]);
    window.__fp_dirty = true;
    console.log("[fp click] added point:", x.toFixed(3), y.toFixed(3), "| total:", window.__fp_points.length);
  }

  function onClearClick() {
    console.log("[fp clear] clearing", window.__fp_points.length, "points");
    window.__fp_points = [];
    window.__fp_dirty = true;
  }

  function attach() {
    const gd = innerPlot();
    if (!gd) {
      console.log("[fp click] waiting for plot...");
      setTimeout(attach, 120);
      return;
    }

    // Remove any existing listeners before adding new ones
    gd.removeEventListener("click", onClick);
    gd.addEventListener("click", onClick);

    const clr = document.getElementById(CLEAR_BTN_ID);
    if (clr) {
      clr.removeEventListener("click", onClearClick);
      clr.addEventListener("click", onClearClick);
    }

    console.log("✓ clicked_points.js attached - mode:", window.__click_mode);
  }

  if (document.readyState === "loading") {
    document.addEventListener("DOMContentLoaded", attach);
  } else {
    attach();
  }
})();