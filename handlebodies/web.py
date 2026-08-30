"""Plotly rendering, for the archive and for 3-D.

Figures written with include_plotlyjs='directory' share one copy of plotly.js
sitting beside the pages, so each figure costs about 10 KB rather than
embedding a 4 MB library.  The result pans, zooms and rotates offline with no
kernel, which is what the notebooks needed Sage's threejs viewer for.
"""

import math

import numpy as np

try:
    import plotly.graph_objects as go
except ImportError as exc:  # pragma: no cover
    raise ImportError(
        "plotly is needed for the web figures: pip install plotly") from exc

PALETTE = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd',
           '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf']
LAYOUT = dict(template='plotly_white', margin=dict(l=40, r=20, t=50, b=40),
              hovermode='closest')


def _colour(j):
    return PALETTE[j % len(PALETTE)]


def _fade(hex_colour, alpha):
    """The same colour at reduced opacity.

    Baked into the colour rather than set as the trace's opacity, so that the
    connecting line can be faint while the computed points stay solid.
    """
    h = hex_colour.lstrip('#')
    r, g, b = (int(h[i:i + 2], 16) for i in (0, 2, 4))
    return f'rgba({r},{g},{b},{alpha})'


# Each step of a traced path is a separate computation, so the points are drawn
# as points; the line only joins them up and is kept deliberately faint.
POINT_SIZE = 4
POINT_SIZE_3D = 2
LINE_ALPHA = 0.35
AXIS_COLOUR = '#444'
#: height of the t axis relative to the widest of Re and Im
Z_ASPECT = 0.25
GRID_COLOUR = 'rgba(0,0,0,0.13)'


def _equal_axes(fig):
    """Mark a figure as wanting equal scaling for its opening view.

    The page turns the recorded bounds into explicit ranges once it knows the
    pixel size of the plot; deliberately not a scaleanchor, which would also
    constrain zooming, and a zoom rectangle should be free to be any shape.

    The bounds are worked out here rather than in the page because plotly
    serialises numpy arrays as base64 typed arrays, which the browser would
    have to decode before it could measure anything.
    """
    xs, ys = [], []
    for t in fig.data:
        for src, dest in ((getattr(t, 'x', None), xs), (getattr(t, 'y', None), ys)):
            if src is None:
                continue
            a = np.asarray(src, dtype=object)
            for v in a.ravel():
                if isinstance(v, (int, float)) and np.isfinite(v):
                    dest.append(float(v))
    meta = dict(equalScale=True)
    if xs and ys:
        meta['bounds'] = [min(xs), max(xs), min(ys), max(ys)]
    fig.update_layout(meta=meta)
    return fig


def points(values, title=None, name='critical values', origin=None,
           origin_name='reference fibre'):
    """Labelled points in the complex plane.

    `origin` marks the reference fibre in its own colour, so it is obvious at a
    glance whether it sits clear of the critical values.
    """
    v = np.asarray([complex(p) for p in values])
    fig = go.Figure(go.Scatter(
        x=v.real, y=v.imag, mode='markers+text',
        text=[str(i) for i in range(len(v))], textposition='top right',
        marker=dict(size=9, color='#d62728'), name=name,
        hovertemplate='%{text}: %{x:.6f} %{y:+.6f}i<extra></extra>'))
    if origin is not None:
        o = complex(origin)
        fig.add_trace(go.Scatter(
            x=[o.real], y=[o.imag], mode='markers',
            marker=dict(size=13, color='#2ca02c', symbol='diamond',
                        line=dict(width=1, color='white')),
            name=origin_name,
            hovertemplate=f'{origin_name}<br>%{{x:.6f}} %{{y:+.6f}}i<extra></extra>'))
    fig.update_layout(title=title, xaxis_title='Re', yaxis_title='Im', **LAYOUT)
    return _equal_axes(fig)


def _loose_points(frames):
    """Every computed value, tagged with how far along the path it came from.

    Used when the values could not be sorted into branches: each step may hold
    a different number of them, so there is nothing to join up and they are
    drawn as they are.
    """
    xs, ys, ts = [], [], []
    n = max(len(frames) - 1, 1)
    for k, frame in enumerate(frames):
        for value in np.asarray(frame):
            v = complex(value)
            xs.append(v.real); ys.append(v.imag); ts.append(k / n)
    return np.array(xs), np.array(ys), np.array(ts)


def trace(branches, title=None, highlight=None, loose=False, zmax=1.0):
    """Tracked branches projected to the complex plane."""
    if loose:
        x, y, t = _loose_points(branches)
        t = t * zmax
        fig = go.Figure(go.Scatter(
            x=x, y=y, mode='markers',
            marker=dict(size=POINT_SIZE, color=t, colorscale='Viridis',
                        showscale=True,
                        colorbar=dict(title='along the path', thickness=12)),
            name='critical values',
            hovertemplate='%{x:.6f} %{y:+.6f}i<extra></extra>'))
        fig.update_layout(title=title, xaxis_title='Re', yaxis_title='Im',
                          **LAYOUT)
        return _equal_axes(fig)
    b = np.asarray(branches)
    fig = go.Figure()
    for j in range(b.shape[1]):
        strong = highlight is not None and j in highlight
        colour = _colour(j)
        fig.add_trace(go.Scatter(
            x=b[:, j].real, y=b[:, j].imag, mode='lines+markers',
            line=dict(color=_fade(colour, LINE_ALPHA if strong else LINE_ALPHA * 0.7),
                      width=1),
            marker=dict(color=colour, size=POINT_SIZE if strong else POINT_SIZE - 1),
            opacity=1.0 if strong else 0.75, name=f'branch {j}',
            hovertemplate=f'branch {j}<br>%{{x:.6f}} %{{y:+.6f}}i<extra></extra>'))
    fig.add_trace(go.Scatter(
        x=b[0].real, y=b[0].imag, mode='markers+text',
        text=[str(j) for j in range(b.shape[1])], textposition='top right',
        marker=dict(size=8, color='#d62728'), name='start'))
    # Hollow, so a closed path's coincident start and end both stay visible.
    fig.add_trace(go.Scatter(x=b[-1].real, y=b[-1].imag, mode='markers',
                             marker=dict(size=13, symbol='circle-open',
                                         line=dict(width=2.5, color='#008b9a')),
                             name='end'))
    fig.update_layout(title=title, xaxis_title='Re', yaxis_title='Im',
                      # clear of the toolbar, which sits in the top right
                      legend=dict(yanchor='top', y=0.86), **LAYOUT)
    # Every branch stays in view; the collision is found by zooming.
    return _equal_axes(fig)


def _floor(fig, b, origin, z0, z1, divisions=6):
    """A grid on the floor plane and a set of axes through the origin.

    plotly's own 3-D grid is drawn on all three walls, which clutters the view
    with vertical lines of constant t that carry no information; so the floor
    is drawn here instead and the built-in grid is switched off.
    """
    lo_x, hi_x = float(b.real.min()), float(b.real.max())
    lo_y, hi_y = float(b.imag.min()), float(b.imag.max())
    o = complex(origin) if origin is not None else 0j
    lo_x, hi_x = min(lo_x, o.real), max(hi_x, o.real)
    lo_y, hi_y = min(lo_y, o.imag), max(hi_y, o.imag)
    px, py = 0.08 * (hi_x - lo_x or 1), 0.08 * (hi_y - lo_y or 1)
    lo_x, hi_x, lo_y, hi_y = lo_x - px, hi_x + px, lo_y - py, hi_y + py

    xs, ys, zs = [], [], []
    for gx in np.linspace(lo_x, hi_x, divisions + 1):
        xs += [gx, gx, None]; ys += [lo_y, hi_y, None]; zs += [z0, z0, None]
    for gy in np.linspace(lo_y, hi_y, divisions + 1):
        xs += [lo_x, hi_x, None]; ys += [gy, gy, None]; zs += [z0, z0, None]
    fig.add_trace(go.Scatter3d(x=xs, y=ys, z=zs, mode='lines',
                               line=dict(color=GRID_COLOUR, width=1),
                               hoverinfo='skip', showlegend=False))

    # x, y and t axes through the reference fibre, as in the original plots.
    # Each spans only its own direction: giving both the larger of the two
    # widths would push the scene's range out past the data and leave the
    # branches squeezed into a fraction of the frame.
    axes = [([lo_x, hi_x], [o.imag, o.imag], [z0, z0]),
            ([o.real, o.real], [lo_y, hi_y], [z0, z0]),
            ([o.real, o.real], [o.imag, o.imag], [z0, z1])]
    for ax, ay, az in axes:
        fig.add_trace(go.Scatter3d(x=ax, y=ay, z=az, mode='lines',
                                   line=dict(color=AXIS_COLOUR, width=2),
                                   hoverinfo='skip', showlegend=False))

    if origin is not None:
        fig.add_trace(go.Scatter3d(
            x=[o.real], y=[o.imag], z=[z0], mode='markers',
            marker=dict(size=4, color='#2ca02c', symbol='diamond'),
            name='reference fibre',
            hovertemplate='reference fibre<br>%{x:.6f} %{y:+.6f}i<extra></extra>'))

    return (hi_x - lo_x), (hi_y - lo_y)


def _eye(ax, ay, az, pull=1.0):
    """Camera position scaled to the aspect box.

    plotly fits the box's bounding sphere in the viewport, so a fixed camera
    leaves an elongated box marooned in white space; the distance has to follow
    the box's own diagonal.
    """
    diagonal = math.sqrt(ax * ax + ay * ay + az * az)
    direction = (0.62, -0.62, 0.30)
    norm = math.sqrt(sum(v * v for v in direction))
    d = pull * diagonal / norm
    return dict(x=direction[0] * d, y=direction[1] * d, z=direction[2] * d)


def trace3d(branches, title=None, highlight=None, zlabel='t', origin=None,
            loose=False, zmax=1.0):
    """Tracked branches with the path parameter as the vertical axis.

    Each horizontal slice is the set of tracked values over one point of the
    path.  Unlike the flat plots, the axes here are stretched to fill the frame
    rather than held to a common scale: in three dimensions a faithful aspect
    ratio squeezes an elongated region into an unreadable sliver, and what the
    picture is for is seeing which branches meet, not measuring distances.
    """
    if loose:
        x, y, t = _loose_points(branches)
        t = t * zmax
        fig = go.Figure(go.Scatter3d(
            x=x, y=y, z=t, mode='markers',
            marker=dict(size=POINT_SIZE_3D, color=t, colorscale='Viridis'),
            name='critical values',
            hovertemplate=(f'%{{x:.6f}} %{{y:+.6f}}i<br>{zlabel}=%{{z:.4f}}'
                           f'<extra></extra>')))
        bare = dict(showgrid=False, showbackground=False, zeroline=False)
        fig.update_layout(
            title=title, height=620,
            scene=dict(xaxis=dict(title='Re', **bare),
                       yaxis=dict(title='Im', **bare),
                       zaxis=dict(title=zlabel, **bare),
                       aspectmode='manual',
                       aspectratio=dict(x=1.0, y=1.0, z=Z_ASPECT),
                       camera=dict(eye=_eye(1.0, 1.0, Z_ASPECT))),
            margin=dict(l=0, r=0, t=42, b=0), template='plotly_white',
            hovermode='closest')
        return fig

    b = np.asarray(branches)
    steps, n = b.shape
    # The vertical axis is the fraction of the path covered, so a run that
    # stopped early shows as reaching only that far up.
    z = np.linspace(0, zmax, steps)

    fig = go.Figure()
    _floor(fig, b, origin, 0.0, zmax)
    for j in range(n):
        strong = highlight is not None and j in highlight
        colour = _colour(j)
        fig.add_trace(go.Scatter3d(
            x=b[:, j].real, y=b[:, j].imag, z=z, mode='lines+markers',
            line=dict(color=_fade(colour, LINE_ALPHA if strong else LINE_ALPHA * 0.7),
                      width=2),
            marker=dict(color=colour,
                        size=POINT_SIZE_3D if strong else POINT_SIZE_3D * 0.75),
            opacity=1.0 if strong else 0.8, name=f'branch {j}',
            hovertemplate=(f'branch {j}<br>%{{x:.6f}} %{{y:+.6f}}i'
                           f'<br>{zlabel}=%{{z:.4f}}<extra></extra>')))
    fig.add_trace(go.Scatter3d(
        x=b[0].real, y=b[0].imag, z=[z[0]] * n, mode='markers+text',
        text=[str(j) for j in range(n)],
        marker=dict(size=3, color='#d62728'), name='start'))
    fig.add_trace(go.Scatter3d(
        x=b[-1].real, y=b[-1].imag, z=[z[-1]] * n, mode='markers',
        marker=dict(size=3, color='#2ca02c'), name='end'))

    # Re and Im each fill the frame; t is kept short so the slices stay legible.
    ar_x = ar_y = 1.0
    bare = dict(showgrid=False, showbackground=False, zeroline=False)
    fig.update_layout(
        title=title,
        scene=dict(
            xaxis=dict(title='Re', **bare),
            yaxis=dict(title='Im', **bare),
            zaxis=dict(title=zlabel, **bare),
            aspectmode='manual',
            aspectratio=dict(x=ar_x, y=ar_y, z=Z_ASPECT),
            camera=dict(eye=_eye(ar_x, ar_y, Z_ASPECT))),
        height=620,
        legend=dict(yanchor='top', y=0.98, xanchor='left', x=0.99),
        margin=dict(l=0, r=0, t=42, b=0), template='plotly_white',
        hovermode='closest')
    return fig
