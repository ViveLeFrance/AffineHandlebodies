from typing import List, Dict
from sage.all import *
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
import numpy as np


    
def NumericalRoots(expr):
    """Returns the numerical roots of the polynomial 'expr'."""
    coeffs = expr.coefficients(sparse=False)
    coeffs = [complex(coefficient) for coefficient in coeffs]
    return np.polynomial.polynomial.polyroots(coeffs)

def set_free_variable_to_one(sol: dict, substitution=1):
    # Identify the free variable
    free_var = None
    for key in sol.keys():
        for variable in sol[key].variables():
            if str(variable).startswith('r'):
                free_var = variable
                # substitute the free variable with 1
                sol[key] = sol[key].subs(free_var==substitution)

    return sol

def set_free_variable_to_one_list(sols: List[dict], substitution=1):
    return [set_free_variable_to_one(sol, substitution=substitution) for sol in sols]


def sort_by_angle(points: List[complex], origin_fibre, anticlockwise: bool = True):
    """Sorts a list of points in the complex plane by their argument,
    starting from the negative real axis."""
    points = [complex(point) for point in points]
    
    points = sorted(points, key=lambda point: (-np.pi+np.angle(point-origin_fibre))%(2*np.pi))
    if not anticlockwise:
        points.reverse()
    return points

def plot_points_ordered(points: List[complex], origin_fibre, title: str = None, anticlockwise=True):
    """Plots a list of points on the complex plane, ordered anticlockwise by argument."""
    # Sage complex type is not compatible with python's, but can be coerced
    points = [complex(point) for point in points]

    # Sort points by argument
    points = sort_by_angle(points, origin_fibre=origin_fibre, anticlockwise=anticlockwise)

  

    real = [point.real for point in points]
    imag = [point.imag for point in points]
    
    
    # subplots returns a tuple (figure, axes), we only need axes
    ax = plt.subplots()[1]

    ax.spines['left'].set_position(('data', 0))
    # Move bottom spine to y=0
    ax.spines['bottom'].set_position(('data', 0))

    # Remove the top and right spines
    ax.spines['right'].set_color('none')
    ax.spines['top'].set_color('none')
    
    ax.set_title(title)
    
    ax.plot(real, imag, 'ro', markersize=5)
    ax.grid(True)

    xlim = ax.get_xlim()
    data_width = xlim[1] - xlim[0]



    for index, point in enumerate(points):
        ax.text(point.real+0.05*data_width, point.imag, str(index), fontsize=12, color='blue')



def plot_path(path: Dict[complex, List[complex]], title: str = None, origin_fibre=0, anticlockwise=True):
    path = {target: [complex(point) for point in path[target]] for target in path.keys()}

    plot_points_origin = path[0]
    plot_points_target = path[1]

    fig, ax = plot_points_ordered(plot_points_origin, title=title, origin_fibre=origin_fibre, anticlockwise=anticlockwise)

    ax.plot([point.real for point in plot_points_target], [point.imag for point in plot_points_target], 'co', markersize=5)

    path_points_x = []
    path_points_y = []
    for step in path.values():

        path_points_x.extend([point.real for point in step])
        
        path_points_y.extend([point.imag for point in step])

    ax.plot(path_points_x, path_points_y, 'bo', markersize=2)

    return fig, ax

def color_generator(n: int):
    cmap = plt.get_cmap('hsv')
    index = 0
    while True:
        yield cmap(index)
        index+=1/n

def plot_path_3d(path: Dict[complex, List[complex]], origin_fibre, title: str = None, anticlockwise=True, trajectory_color = 'blue'):
    path = {target: [complex(point) for point in path[target]] for target in path.keys()}


    points = []

    for index, step in enumerate(path.values()):
        
        for point in step:
            points.append((point.real, point.imag, index/len(path)))

    n = len(path[0])


    color = color_generator(n)

    sorted_points = sort_by_angle(path[0], origin_fibre=origin_fibre, anticlockwise=anticlockwise)
    sorted_points = [complex(pt) for pt in sorted_points]


    # Now, generate labels in order.
    labels = []
    for i, pt in enumerate(sorted_points):
        # Place the label near the point (using its x,y coordinate and z=0).
        label = text3d(str(i), (pt.real, pt.imag, 0.05), color='purple', fontsize=30, fontweight='bold')
        labels.append(label)


    init_points  =[(point.real, point.imag, 0) for point in path[0]]
    final_points = [(point.real, point.imag, 1) for point in path[1]]
    init_plot = point3d(init_points, size=20, color='red')
    target_plot = point3d(final_points, size=20, color='green')
    trajectory_plot = point3d(points, size=5, color=trajectory_color)

    axis_length = 2

    # Create axes (x, y, z) at the origin
    x_axis = line3d([(0,0,0), (axis_length,0,0)], color='black', thickness=1.5)
    y_axis = line3d([(0,0,0), (0,axis_length,0)], color='black', thickness=1.5)
    z_axis = line3d([(0,0,0), (0,0,axis_length)], color='black', thickness=1.5)

    axes = x_axis + y_axis + z_axis

    plot = trajectory_plot + init_plot + target_plot + axes
    for lab in labels:
        plot += lab
   
    if title:
        title_text = text3d(title, (0, 0, axis_length + 0.5), color='black', fontsize=20)
        plot += title_text

    plot.show()

def pl_path(points: List[complex], steps=70):
    path = []

    if len(points) == 1:
        return [points[0]]*steps

    seg_length = steps//(len(points)-1)
    remainder = steps - seg_length * (len(points)-1) 

    for i in range(0, len(points)-1):
        for s in np.linspace(0,1,seg_length):
            path.append((1-s)*points[i] + s*points[(i+1)])
    
    while remainder > 0:
        path.append(points[-1])
        remainder -= 1


    return path

def pl_path_1(origin_fibre, target_fibre, offset = None, steps=70, above=True):
    if offset is None:
        offset = np.abs(target_fibre - origin_fibre)/4

    theta = np.arctan(1/2)
    hyp = offset / np.sin(theta)

    if above is True:
        intermediate_point = origin_fibre + hyp*np.exp(1j*theta)
    else:
        intermediate_point = origin_fibre - hyp*np.exp(1j*theta)

    return pl_path([origin_fibre, intermediate_point, target_fibre], steps=steps)


def plot_pl_path(path: List[complex], steps=70, title=None, spec_points=None):
    points = [complex(point) for point in path]
    real = [point.real for point in points]
    imag = [point.imag for point in points]

    fig, ax = plt.subplots()
    ax.plot(real, imag, 'bo-', markersize=4, label="PL Path")
    ax.set_xlabel('Re')
    ax.set_ylabel('Im')
    ax.set_title(title if title else "Piecewise Linear Path in Complex Plane")
    ax.grid(True)

    # Overlay specified points if provided
    overlay_points = spec_points
    if overlay_points is not None:
        overlay_points = [complex(pt) for pt in overlay_points]
        overlay_real = [pt.real for pt in overlay_points]
        overlay_imag = [pt.imag for pt in overlay_points]
        ax.plot(overlay_real, overlay_imag, 'ro', markersize=7, label="Overlay Points")

    ax.legend()
    plt.show()


def plot_paths(paths: Dict[int, List[complex]], origin_fibre_rho=0):
     # Sage complex type is not compatible with python's, but can be coerced
    points = [complex(point) for point in paths[0]]

    # Sort points by argument
    points = sort_by_angle(points, origin_fibre=origin_fibre_rho, anticlockwise=True)

  

    real = [point.real for point in points]
    imag = [point.imag for point in points]
    
    # if ax is None:
    fig, ax = plt.subplots()

    ax.spines['left'].set_position(('data', 0))
    # Move bottom spine to y=0
    ax.spines['bottom'].set_position(('data', 0))

    # Remove the top and right spines
    ax.spines['right'].set_color('none')
    ax.spines['top'].set_color('none')

    
    
    ax.set_title("trace of critical values for chosen parameter ranges")
    
    ax.plot(real, imag, 'ro', markersize=5)
    ax.grid(True)

    xlim = ax.get_xlim()
    data_width = xlim[1] - xlim[0]



    for index, point in enumerate(points):
        ax.text(point.real+0.05*data_width, point.imag, str(index), fontsize=12, color='blue')

    paths.pop(0)

    cmap = plt.get_cmap('viridis')

    gradient = range(len(paths))
    norm = Normalize(vmin=min(gradient), vmax=max(gradient))

    for index, path in paths.items():
        
        real = [complex(point).real for point in path]
        imag = [complex(point).imag for point in path]

        color = cmap(norm(index))

        ax.scatter(real, imag, color=color, norm=norm, s=1)

    target_points = paths[len(paths)]
    points = [complex(point) for point in target_points]

    real = [point.real for point in points]
    imag = [point.imag for point in points]

    ax.plot(real, imag, 'mo', markersize=5)
        
    plt.show()


