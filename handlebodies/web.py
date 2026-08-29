"""Plotly rendering, for the archive and for 3-D.

Figures written with include_plotlyjs='directory' share one copy of plotly.js
sitting beside the pages, so each figure costs about 10 KB rather than
embedding a 4 MB library.  The result pans, zooms and rotates offline with no
kernel, which is what the notebooks needed Sage's threejs viewer for.
"""

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
LINE_ALPHA = 0.35


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
    return fig


def _frame(b, highlight, pad=0.25):
    """Axis ranges that show the collision rather than the whole excursion.

    A branch that swings far away compresses the collision to a point at the
    default scale, which is precisely what the picture exists to show.  So the
    initial view is the box around the colliding branches; plotly's autoscale
    button restores the full extent.  Used for the flat traces only -- the 3-D
    views show everything.
    """
    if not highlight:
        return None, None
    sel = b[:, list(highlight)]
    xs, ys = sel.real, sel.imag
    dx = max(xs.max() - xs.min(), 1e-9)
    dy = max(ys.max() - ys.min(), 1e-9)
    span = max(dx, dy)          # keep it square-ish so shape is not distorted
    cx, cy = (xs.max() + xs.min()) / 2, (ys.max() + ys.min()) / 2
    h = span * (0.5 + pad)
    return [cx - h, cx + h], [cy - h, cy + h]


def trace(branches, title=None, highlight=None):
    """Tracked branches projected to the complex plane."""
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
    fig.update_layout(title=title, xaxis_title='Re', yaxis_title='Im', **LAYOUT)
    xr, yr = _frame(b, highlight)
    if xr:
        fig.update_xaxes(range=xr)
        fig.update_yaxes(range=yr)
    return fig


def trace3d(branches, title=None, highlight=None, zlabel='t'):
    """Tracked branches with the path parameter as the vertical axis.

    This is the picture the notebooks built with Sage's threejs viewer: each
    horizontal slice is the set of tracked values over one point of the path.
    It applies equally to the critical values of rho over a path in the base of
    pi (a matching path) and to the points of a fibre of rho over a path in the
    base of rho.
    """
    b = np.asarray(branches)
    steps, n = b.shape
    z = np.linspace(0, 1, steps)

    fig = go.Figure()
    for j in range(n):
        strong = highlight is not None and j in highlight
        colour = _colour(j)
        fig.add_trace(go.Scatter3d(
            x=b[:, j].real, y=b[:, j].imag, z=z, mode='lines+markers',
            line=dict(color=_fade(colour, LINE_ALPHA if strong else LINE_ALPHA * 0.7),
                      width=2),
            marker=dict(color=colour, size=POINT_SIZE - 1 if strong else POINT_SIZE - 2),
            opacity=1.0 if strong else 0.8, name=f'branch {j}',
            hovertemplate=(f'branch {j}<br>%{{x:.6f}} %{{y:+.6f}}i'
                           f'<br>{zlabel}=%{{z:.4f}}<extra></extra>')))
    fig.add_trace(go.Scatter3d(
        x=b[0].real, y=b[0].imag, z=[z[0]] * n, mode='markers+text',
        text=[str(j) for j in range(n)],
        marker=dict(size=4, color='#d62728'), name='start'))
    fig.add_trace(go.Scatter3d(
        x=b[-1].real, y=b[-1].imag, z=[z[-1]] * n, mode='markers',
        marker=dict(size=4, color='#2ca02c'), name='end'))
    # Deliberately unclamped: every critical value of rho and its whole
    # trajectory has to stay visible here, since the point of the 3-D view is
    # to see them all move together.  Rotating and zooming is how you inspect
    # a particular pair.
    scene = dict(xaxis_title='Re', yaxis_title='Im', zaxis_title=zlabel)
    fig.update_layout(title=title, scene=scene, **LAYOUT)
    return fig

