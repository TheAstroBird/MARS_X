import matplotlib.pyplot as plt
import pandas as pd

def plot_2d_models(path: str='dynamic/integral_parameters',
                   x: str='Love_number',
                   x_real: bool=True,
                   y: str='inertia',
                   y_real: bool=True,
                   color_axis: bool=True,
                   c: str='Chandler_period',
                   c_real: bool=True,
                   compare: bool=False,
                   c_comp: str='Chandler_period_elastic',
                   c_comp_real: bool=True,
                   save: bool=False,
                   target: str='png/test.png'):
    """Plots models integrals along 2 axis and 1 'colorbar'-axis

    :param path: The path of file with models integrals to plot, defaults to 'dynamic/integral_parameters'
    :type path: str
    :param x: The name of parameter from the file to draw on X-axis, defaults to 'Love_number'
    :type x: str
    :param x_real: Defines whether to plot real part, if parameter x is complex, defaults to True
    :type x_real: bool
    :param y: The name of parameter from the file to draw on Y-axis, defaults to 'inertia'
    :type y: str
    :param y_real: Defines whether to plot real part, if parameter y is complex, defaults to True
    :type y_real: bool
    :param color_axis: Defines whether to plot a 'colorbar'-axis, defaults to True
    :type color_axis: bool
    :param c: The name of parameter from the file to draw on 'colorbar'-axis, defaults to 'Chandler_period'
    :type c: str
    :param c_real: Defines whether to plot real part, if parameter c is complex, defaults to True
    :type c_real: bool
    :param compare: Defines whether to plot a secondary graph with another parameter at 'colorbar'-axis.
        The parameter is defined by variable c_comp. Defaults to False
    :type compare: bool
    :param c_comp: The name of parameter from the file to draw on secondary graph along 'colorbar'-axis,
        defaults to 'Chandler_period_elastic'
    :type c_comp: str
    :param c_comp_real: Defines whether to plot real part, if parameter c_comp is complex, defaults to True
    :type c_comp_real: bool
    :param save: Defines whether to save the plot, defaults to False
    :type save: bool
    :param target: The path of the file of save the plot, defaults to 'png/test.png'
    :type target: str"""
    path = '../data/' + path + '.xlsx'
    df = pd.read_excel(path)
    for k, r in zip([x, y, c, c_comp], [x_real, y_real, c_real, c_comp_real]):
        if k in ['Love_number', 'Love_number_CW']:
            df[k] = df[k].apply(lambda z: complex(z).real) if r else df[k].apply(lambda z: complex(z).imag)
    x_axis = df[x]
    y_axis = df[y]
    c_axis = df[c]
    c_comp_axis = df[c_comp]

    cols = 1
    x_figsize = 7
    if compare and color_axis:
        cols = 2
        x_figsize = 14
    fig, ax = plt.subplots(1, cols, figsize=(x_figsize, 5))
    if not (compare and color_axis):
        ax = [ax]

    sc = []
    if color_axis:
        sc.append(ax[0].scatter(x_axis, y_axis, c=c_axis, cmap='plasma'))
    else:
        sc.append(ax[0].scatter(x_axis, y_axis, c='k'))
    if color_axis and compare:
        sc.append(ax[1].scatter(x_axis, y_axis, c=c_comp_axis, cmap='plasma'))
    for a, c in zip(ax, sc):
        if x == 'Love_number':
            a.axvline(0.166, c='k')
            a.axvline(0.182, c='k')
            a.set_xlabel('k$_2$')
        elif x == 'inertia':
            a.axvline(0.3634, c='k')
            a.axvline(0.3646, c='k')
            a.set_xlabel('I/(MR$^2$)')
        if y == 'Love_number':
            a.axhline(0.166, c='k')
            a.axhline(0.182, c='k')
            a.set_ylabel('k$_2$')
        elif y == 'inertia':
            a.axhline(0.3634, c='k')
            a.axhline(0.3646, c='k')
            a.set_ylabel('I/(MR$^2$)')
        a.tick_params(axis='x', labelrotation=45)
        a.grid()
        if color_axis:
            fig.colorbar(c, ax=a)
    plt.show()
    if save:
        target = '../data/' + target
        fig.savefig(target, dpi=600, bbox_inches='tight')
    pass