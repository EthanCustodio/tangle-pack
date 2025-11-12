# tanglepack_webdash/callbacks/bridge_ops.py
from __future__ import annotations
import numpy as np
import traceback
from dash import (
    Dash,
    Input,
    Output,
    State as DashState,
    no_update,
    ctx,
    ClientsideFunction,
)
from ..sessions import get_state
from ..utils.figures import (
    plot_tangle,
    plot_all_bridges,
    plot_intersections,
    plot_bridges_with_highlight,
)


def register(app: Dash):

    # Update the button state and info text when selection changes
    @app.callback(
        Output("selected-bridge-info", "children"),
        Output("btn-iterate-bridge", "disabled"),
        Input("selected-bridge-idx", "data"),
        DashState("sid", "data"),
        prevent_initial_call=False,
    )
    def update_selection_ui(selected_idx, sid):
        """Update UI elements when bridge selection changes."""

        if selected_idx is None:
            return "No bridge selected", True

        if not sid:
            return "No session", True

        state = get_state(sid)
        bridges = getattr(state, "bridges", None)

        if not bridges:
            return "No bridges available", True

        bridge_list = list(bridges) if isinstance(bridges, dict) else bridges

        if 0 <= selected_idx < len(bridge_list):
            return f"Bridge {selected_idx + 1} selected", False
        else:
            return "Invalid bridge index", True

    # Update plot to highlight selected bridge
    @app.callback(
        Output("plot", "figure", allow_duplicate=True),
        Input("selected-bridge-idx", "data"),
        DashState("sid", "data"),
        prevent_initial_call=True,
    )
    def highlight_selected_bridge(selected_idx, sid):
        """Redraw the plot with the selected bridge highlighted."""

        if not sid:
            return no_update

        state = get_state(sid)

        if state.workbench is None or state.fp is None:
            return no_update

        bridges = getattr(state, "bridges", None)
        if not bridges:
            return no_update

        # Rebuild the figure
        state.fig = plot_tangle(state, state.workbench, state.fp, "stable")
        state.fig = plot_intersections(state, state.workbench, state.fp)

        bridge_list = list(bridges) if isinstance(bridges, dict) else bridges

        # Plot bridges with highlighting
        state.fig = plot_bridges_with_highlight(
            state, bridge_list, highlight_idx=selected_idx
        )

        return state.fig

    # Clear bridge selection
    @app.callback(
        Output("selected-bridge-idx", "data", allow_duplicate=True),
        Input("btn-clear-bridge-selection", "n_clicks"),
        prevent_initial_call=True,
    )
    def clear_bridge_selection(_n):
        """Clear the current bridge selection."""
        # Also clear the JS variable
        return None

    # Iterate the selected bridge
    @app.callback(
        Output("plot", "figure", allow_duplicate=True),
        Output("status", "children", allow_duplicate=True),
        Input("btn-iterate-bridge", "n_clicks"),
        DashState("sid", "data"),
        DashState("selected-bridge-idx", "data"),
        prevent_initial_call=True,
    )
    def iterate_selected_bridge(_n, sid, selected_idx):
        """Iterate the currently selected bridge forward."""

        if not sid:
            return no_update, "❌ Missing session id."

        state = get_state(sid)

        if state.workbench is None or state.fp is None:
            return no_update, "⚠️ Build a system & fixed point first."

        bridges = getattr(state, "bridges", None)
        if not bridges:
            return no_update, "⚠️ No bridges available."

        if selected_idx is None:
            return no_update, "⚠️ No bridge selected."

        bridge_list = list(bridges) if isinstance(bridges, dict) else bridges

        if selected_idx < 0 or selected_idx >= len(bridge_list):
            return no_update, "❌ Invalid bridge index."

        try:
            # Get the selected bridge
            selected_bridge = bridge_list[selected_idx]

            # Iterate the bridge using ManifoldMachine
            iterated_bridge = state.workbench._man_machine.iterate_bridge(
                selected_bridge
            )

            # Replace the bridge in the list
            # bridge_list[selected_idx] = iterated_bridge
            bridge_list.append(iterated_bridge)
            state.workbench.Tangle.add_manifold(iterated_bridge)

            # Update state
            if isinstance(bridges, dict):
                # If bridges is a dict, update it accordingly
                bridge_keys = list(bridges.keys())
                bridges[bridge_keys[selected_idx]] = iterated_bridge
            else:
                state.bridges = bridge_list

            # Redraw the figure with the iterated bridge highlighted
            state.fig = plot_tangle(state, state.workbench, state.fp, "stable")
            state.fig = plot_intersections(state, state.workbench, state.fp)
            state.fig = plot_bridges_with_highlight(
                state, bridge_list, highlight_idx=selected_idx
            )

            return state.fig, f"✅ Iterated bridge {selected_idx + 1}."

        except Exception as e:
            return no_update, f"❌ Iteration error: {e}\n{traceback.format_exc()}"

    # Push bridge data to JavaScript when bridges are created/updated
    @app.callback(
        Output("bridges-data", "data"),  # Store the data
        Input("plot", "figure"),
        DashState("sid", "data"),
        prevent_initial_call=True,
    )
    def sync_bridges_to_js(fig, sid):
        """Sync bridge data to JavaScript for client-side selection."""
        if not sid:
            return None

        state = get_state(sid)
        bridges = getattr(state, "bridges", None)

        if not bridges:
            return None

        bridge_list = list(bridges) if isinstance(bridges, dict) else bridges

        # Prepare bridge data for JS
        bridges_data = []
        for bridge in bridge_list:
            try:
                arr = bridge.get_point_array()
                if arr.size > 0:
                    bridges_data.append({"points": arr.tolist()})
            except Exception:
                continue

        return bridges_data if bridges_data else None

    # Clientside callback to push bridge data to window.__bridges_data
    app.clientside_callback(
        ClientsideFunction(namespace="bridge", function_name="push_bridges_data"),
        Output("bridges-data-holder", "children"),
        Input("bridges-data", "data"),
        prevent_initial_call=False,
    )
