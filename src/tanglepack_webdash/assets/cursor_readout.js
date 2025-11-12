// src/tanglepack_webdash/assets/cursor_readout.js
(function () {
  const READOUT_ID = "cursor-readout";
  const PLOT_WRAPPER = "#plot";
  const RETRY_MS = 120;

  function fmt(v) {
    if (!Number.isFinite(v)) return "—";
    return Math.abs(v) >= 1e4 ? v.toExponential(3) : v.toFixed(3);
  }

  function attach() {
    const wrap = document.querySelector(PLOT_WRAPPER);
    const gd = wrap && wrap.querySelector(".js-plotly-plot");
    const readout = document.getElementById(READOUT_ID);
    if (!gd || !readout || !gd._fullLayout || !gd._fullLayout.xaxis || !gd._fullLayout.yaxis) {
      setTimeout(attach, RETRY_MS);
      return;
    }

    function onMove(evt) {
      const full = gd._fullLayout;
      const rect = gd.getBoundingClientRect();
      const px = evt.clientX - rect.left - full.margin.l;
      const py = evt.clientY - rect.top  - full.margin.t;
      const w = (full._size && full._size.w) || (gd.clientWidth  - full.margin.l - full.margin.r);
      const h = (full._size && full._size.h) || (gd.clientHeight - full.margin.t - full.margin.b);
      if (px < 0 || py < 0 || px > w || py > h) return;
      const x = full.xaxis.p2c(px);
      const y = full.yaxis.p2c(py);
      readout.textContent = `x: ${fmt(x)} | y: ${fmt(y)}`;
    }

    gd.removeEventListener("mousemove", onMove);
    gd.addEventListener("mousemove", onMove);

    // keep stable after relayouts/resizes
    if (gd.on) {
      gd.on("plotly_relayout", () => {
        // handler reads fresh _fullLayout each time, so nothing else needed
      });
    }
  }

  if (document.readyState === "loading") {
    document.addEventListener("DOMContentLoaded", attach);
  } else {
    attach();
  }
})();
