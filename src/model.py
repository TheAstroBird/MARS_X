import math as m
import numpy as np
import scipy.integrate as ing
import scipy.interpolate as inp
import scipy.optimize as opt
import scipy.linalg as lin
import pandas as pd
import matplotlib.pyplot as plt

class Mars:
    """A class of models of internal structure of Mars. Contains all parameters to calculate for model"""

    def __init__(self,
                 composition: str = 'BF97',
                 areoterm: str = "ATM",
                 density_crust: float = 2.9,
                 depth_crust: float = 50,
                 density_core: float = 6.1,
                 hydro_core: float = 0,
                 sulfur_core: float = 0.8,
                 viscosity: float = 10**21,
                 viscosity_melt: float = 10**11,
                 andrade: float = 0.3,
                 bconstcore: bool = False,
                 bviscosity: bool = True,
                 bmeltlayer: bool = True,
                 bcreepfunction: bool = False
                ):
        """Creates an instance of internal structure model that about to be calculated

            :param composition: The name of the composition used. Defaults to "BF97"
            :type composition: str
            :param areoterm: The name of the areoterm used. Defaults to "ATM"
            :type areoterm: str
            :param density_crust: Defines crust density of the model, expressed in g/cm3. Defaults to 2.9
            :type density_crust: float
            :param depth_crust: Defines crust depth of the model, expressed in km. Defaults to 50
            :type depth_crust: float
            :param density_core: Is used only when bconstcore = True!
                Defines core density of the model, expressed in g/cm3. Defaults to 6.1
            :type density_core: float
            :param hydro_core: Is used only when bconstcore = False (by default)!
                Defines the molar fraction of FeH in the core. Defaults to 0
            :type hydro_core: float
            :param sulfur_core: Is used only when bconstcore = False (by default)!
                Defines the molar fraction of FeS in the core. Defaults to 0.8
            :type sulfur_core: float
            :param viscosity: Is used only when bviscosity = True (by default)!
                Defines the viscosity parameter of the model, expressed in Pa*s. Defaults to 10**21
            :type viscosity: float
            :param viscosity_melt: Is used only when bviscosity = True and bmeltlayer = True (by default)!
                Defines the viscosity of the melt layer, expressed in Pa*s. Defaults to 10**11
            :type viscosity_melt: float
            :param andrade: Is used only when bviscosity = True (by default)!
                Defines the Andrade parameter of the model. Defaults to 0.3
            :type andrade: float
            :param bconstcore: Defines whether the model is with the constant density in the core.
                Density of the core is defined by parameter "density_core". Defaults to False
            :type bconstcore: bool
            :param bviscosity: Defines whether the model accounts for inelasticity.
                Parameters of inelasticity are defined through "viscosity" and "andrade". Defaults to True
            :type bviscosity: bool
            :param bmeltlayer: Defines whether the model containts the melt layer at the bottom of the mantle.
                Viscosity of the melt layer is defined through parameter "viscosity_melt". Defaults to True
            :type bmeltlayer: bool
            :param bcreepfunction: Defines whether there is a need in calculating inelastic value of Chandler period through
                creep function. Defaults to False
            :type bcreepfunction: bool
            """
        if sulfur_core + hydro_core > 1:
            raise ValueError("Total amount of FeH and FeS cannot exceed 1!")
        self.composition = composition.lower()
        self.areoterm = areoterm.lower()
        self.density_crust = density_crust
        self.depth_crust = depth_crust
        self.density_core = density_core
        self.hydro_core = hydro_core
        self.sulfur_core = sulfur_core
        self.ferrum_core = 1 - sulfur_core - hydro_core
        self.viscosity = viscosity
        self.viscosity_melt = viscosity_melt
        self.andrade = andrade
        self.bconstcore = bconstcore
        self.bviscosity = bviscosity
        self.bmeltlayer = bmeltlayer
        self.bcreepfunction = bcreepfunction

        self.R_Mars: float = 3389.92
        self.M_Mars = 6.4185e23  # kg
        self.den_av = 3 * self.M_Mars / (4 * m.pi * self.R_Mars ** 3) * 1e-12  # g/cm3
        self.Gravity_const = 6.6743e-11  # in SI
        self.grav_av = 3.7279  # м/с^2x
        self.P_const = 3 * self.Gravity_const * self.M_Mars ** 2 / (4 * m.pi * self.R_Mars ** 4) * 1e-21
        self.T_exp_const = 3e-5 * self.Gravity_const * self.M_Mars / 5e2 / self.R_Mars / 1e3
        self.H_melt = 200  # km - thickness of the melt layer
        self.freq_cw = 2 * m.pi / 206.9 / 86400  # tidal frequency - Chandler frequency, 1/s
        self.freq_sun = 2 * m.pi / 44340  # Sun tidal frequency, 1/s
        self.n_grid = 10000 # number of steps of the grid
        self.grid_step = 1 / self.n_grid
        self.t_core_bound = 0.02  # stop point in the core (defines the limit of the integration, then the density considered to be constant)
        self.creep_exp = 0.4

        self.density_core_boundary: float | None = None
        self.radius_core: float | None = None
        self._radius_core_dl: float | None = None
        self.mass_model: float | None = None
        self.MOI: float | None = None
        self.MOI_core: float | None = None
        self.Molodensky_number: float | None = None
        self.Love_elast: float | None = None
        self.Love_cw: float | None = None
        self.Love_sun: float | None = None
        self.Love_sec: float | None = None
        self.period_cw_elast: float | None = None
        self.period_cw: float | None = None
        self.period_cw_creep: float | None = None
        self.period_cw_creep_delta: float | None = None
        self.l_crust: int | None = None
        self.l_to_core: int | None = None
        self.l_planet: int | None = None

        self.radius: list[float] | None = None
        self._radius_dl: list[float] | None = None
        self.density: list[float] | None = None
        self._density_dl: list[float] | None = None
        self.gravity: list[float] | None = None
        self._gravity_dl: list[float] | None = None
        self.mass: list[float] | None = None
        self._mass_dl: list[float] | None = None
        self.pressure: list[float] | None = None
        self.temperature: list[float] | None = None
        self.vp: list[float] | None = None
        self.vs: list[float] | None = None
        self.visc_distr: list[float] | None = None
        self.lame1_elast: list[float] | None = None
        self._lame1_elast_dl: list[float] | None = None
        self.lame2_elast: list[float] | None = None
        self._lame2_elast_dl: list[float] | None = None
        self.lame1: list[float] | None = None
        self._lame1_dl: list[float] | None = None
        self.lame2: list[float] | None = None
        self._lame2_dl: list[float] | None = None
        self.lame1_cw: list[float] | None = None
        self._lame1_cw_dl: list[float] | None = None
        self.lame2_cw: list[float] | None = None
        self._lame2_cw_dl: list[float] | None = None
        self.Q: list[float] | None = None

    def ced(mars):
        """'Calculate Elastic Distributions'

            Calculates the distributions of main parameters such as density, pressure, temperature, mass and gravity"""

        # Initialization of the main parameters of the method
        print('\n')
        print("-------------------------------------------------------")
        print("Initialization the calculation of elastic distributions")
        print("-------------------------------------------------------")
        rad_core_const = 1650  # km

        mars._radius_dl = [1.0]
        mars._density_dl = []
        mars.gravity = []
        mars._mass_dl = [1.0]
        mars.pressure = [0.0]
        mars.temperature = [300.0]
        mars.vp = []
        mars.vs = []

        den_crust_const = mars.density_crust / mars.den_av
        den_core_const = mars.density_core / mars.den_av
        depth_crust_const = mars.depth_crust / mars.R_Mars
        rad_core_const /= mars.R_Mars

        DATAs_perplex = pd.read_excel('../data/dynamic/mantle_distributions_' + mars.composition + '_' + mars.areoterm
                                      + '.xlsx')
        DATAs_perplex['density'] /= mars.den_av * 1e3
        DATAs_perplex['pressure'] /= 1e4

        mix = mars.ferrum_core * 55.85 + mars.sulfur_core * 87.92 + mars.hydro_core * 56.86
        x_fe = mars.ferrum_core * 55.85 / mix
        x_fes = mars.sulfur_core * 87.92 / mix
        x_feh = mars.hydro_core * 56.86 / mix

        # Calculating distributions through depth
        class system15:
            def __init__(self, dens, temp):
                self.dens = dens
                self.temp = temp

            def __call__(self, t, y):
                der1 = 3 * t ** 2 * self.dens(t, y)
                der2 = - mars.P_const * self.dens(t, y) * y[0] / t ** 2
                der3 = self.temp(t, y)
                return [der1, der2, der3]

        # 1 - Crust
        mars._density_dl.append(den_crust_const)

        def den_crust_func(t, y):
            return den_crust_const

        def temp_crust_func(t, y):
            fun0 = inp.interp1d(DATAs_perplex['pressure'], DATAs_perplex['temperature'],
                                fill_value='extrapolate')
            fun = inp.interp1d([1 - depth_crust_const, 1], [float(fun0(y[1])), 300])
            return float(fun(t))

        sys15_crust = system15(den_crust_func,
                               lambda t, y: 0)  # temperature is calculated by interpolation and not by dif. equation
        crust_distr = ing.RK45(sys15_crust, 1, [1, 0, 300], 1 - depth_crust_const, max_step=mars.grid_step)
        while crust_distr.status == 'running':
            crust_distr.step()
            mars._radius_dl.append(crust_distr.t)
            mars._density_dl.append(den_crust_const)
            mars._mass_dl.append(crust_distr.y[0])
            mars.pressure.append(crust_distr.y[1])
            mars.temperature.append(temp_crust_func(crust_distr.t, crust_distr.y))
        mars.l_crust = len(mars._radius_dl)

        # 2 - Mantle
        print("0% - calculating of the mantle")
        def den_mantle_func(t, y):
            fun = inp.interp1d(DATAs_perplex['pressure'], DATAs_perplex['density'], fill_value='extrapolate')
            return float(fun(y[1]))

        def temp_mantle_func(t, y):
            fun = inp.interp1d(DATAs_perplex['pressure'], DATAs_perplex['temperature'], fill_value='extrapolate')
            return float(fun(y[1]))

        mars._radius_dl.append(mars._radius_dl[-1] - mars.grid_step ** 2)
        mars._density_dl.append(den_mantle_func(mars._radius_dl[-1],
                                                        [mars._mass_dl[-1],
                                                         mars.pressure[-1] * (1 + mars.grid_step ** 2),
                                                         mars.temperature[-1]]))
        mars._mass_dl.append(mars._mass_dl[-1])
        mars.pressure.append(mars.pressure[-1] * (1 + mars.grid_step ** 2))
        mars.temperature.append(mars.temperature[-1])

        sys15_mantle = system15(den_mantle_func,
                                lambda t, y: 0)  # temperature is calculated by interpolation and not by dif. equation
        mantle_distr = ing.RK45(sys15_mantle, mars._radius_dl[-1],
                                [mars._mass_dl[-1], mars.pressure[-1],
                                 mars.temperature[-1]],
                                rad_core_const, max_step=mars.grid_step)

        while mantle_distr.status == 'running':
            mantle_distr.step()
            mars._radius_dl.append(mantle_distr.t)
            mars._density_dl.append(den_mantle_func(mantle_distr.t, mantle_distr.y))
            mars._mass_dl.append(mantle_distr.y[0])
            mars.pressure.append(mantle_distr.y[1])
            mars.temperature.append(temp_mantle_func(mantle_distr.t, mantle_distr.y))

        # 3 - Core (first calculation)
        class Core_state:
            def __init__(self, temp, mode):
                self.temp = temp
                self.mode = mode
                # parameter lists [den_0, alpha_0, K_T_0, K'_T_0, Kdot_T_0]
                if mode == 'Fe':
                    self.param = [8 / mars.den_av, 75e-6, 183, 4, 0]
                elif mode == 'FeS':
                    self.param = [4.94 / mars.den_av, 68.52e-6, 54, 4, 0.02]
                elif mode == 'FeH':
                    self.param = [6.7 / mars.den_av, 0, 121, 5.31, 0]

            def __call__(self, den: np.ndarray):
                den = den[0]
                if self.mode == 'Fe' or self.mode == 'FeS':
                    f = 1 / 2 * ((den / self.param[0]) ** (2 / 3) - 1)
                    p0 = 3 * f * (1 + 2 * f) ** (5 / 2) * self.param[2] * (1 + 3 / 2 * (self.param[3] - 4) * f)
                    K = self.param[2] + self.param[4] * (self.temp - 800)
                    p = p0 + self.param[1] * K * (self.temp - 800)
                    return p
                elif self.mode == 'FeH':
                    f: float = 1 - (self.param[0] / den) ** (1 / 3)
                    p = 3 * self.param[2] * (den / self.param[0]) ** (2 / 3) * f * m.exp(3 / 2 * (self.param[3] - 1) * f)
                    return p

        def den_core_func(t, y: np.ndarray) -> float:
            if mars.bconstcore:
                return den_core_const
            else:
                pfe = Core_state(y[2], 'Fe')
                pfes = Core_state(y[2], 'FeS')
                pfeh = Core_state(y[2], 'FeH')
                den_fe: float = opt.fsolve(lambda k: pfe(k) - y[1], np.array(8.5 / mars.den_av))[0]
                den_fes: float = opt.fsolve(lambda k: pfes(k) - y[1], np.array(6 / mars.den_av))[0]
                den_feh: float = opt.fsolve(lambda k: pfeh(k) - y[1], np.array(7 / mars.den_av))[0]
                return 1 / (x_fe / den_fe + x_fes / den_fes + x_feh / den_feh)

        def temp_core_func(t, y):
            return -mars.T_exp_const * y[0] * y[2] / t ** 2

        distr_core = {'density': [den_core_func(mars._radius_dl[-1] - mars.grid_step ** 2,
                                                np.array([mars._mass_dl[-1],
                                                 mars.pressure[-1] * (1 + mars.grid_step ** 2),
                                                 mars.temperature[-1]]))],
                      'radius': [mars._radius_dl[-1] - mars.grid_step ** 2],
                      'mass': [mars._mass_dl[-1]],
                      'pressure': [mars.pressure[-1] * (1 + mars.grid_step ** 2)],
                      'temperature': [mars.temperature[-1]]}
        sys15_core = system15(den_core_func, temp_core_func)
        core_distr = ing.RK45(sys15_core, distr_core['radius'][-1],
                              [distr_core['mass'][-1], distr_core['pressure'][-1], distr_core['temperature'][-1]],
                              mars.t_core_bound)

        while core_distr.status == 'running':
            core_distr.step()
            distr_core['radius'].append(core_distr.t)
            distr_core['density'].append(den_core_func(core_distr.t, core_distr.y))
            distr_core['mass'].append(core_distr.y[0])
            distr_core['pressure'].append(core_distr.y[1])
            distr_core['temperature'].append(core_distr.y[2])

        distr_core['radius'].append(0)
        distr_core['density'].append(distr_core['density'][-1])
        distr_core['mass'].append(distr_core['mass'][-1] - distr_core['density'][-1] * mars.t_core_bound ** 3)
        distr_core['pressure'].append(
            distr_core['pressure'][-1] + mars.P_const * (distr_core['density'][-1] * mars.t_core_bound) ** 2 / 2)
        distr_core['temperature'].append(
            distr_core['temperature'][-1] * m.exp(mars.T_exp_const * distr_core['density'][-1] * mars.t_core_bound ** 2 / 2))

        mars.l_to_core = len(mars._radius_dl)

        # 4 - "jump" - approximate shift of the boundary according to the mass difference in the center of the planet
        r3_c = mars._radius_dl[-1] ** 3 + distr_core['mass'][-1] / (
                    distr_core['density'][0] - mars._density_dl[-1])
        r_c = r3_c ** (1 / 3) if r3_c >= 0 else -(-r3_c) ** (1 / 3)

        while mars._radius_dl[-1] < r_c:
            for l in [mars._radius_dl, mars._density_dl, mars._mass_dl, mars.pressure, mars.temperature]:
                l.pop()

        mantle_distr = ing.RK45(sys15_mantle, mars._radius_dl[-1],
                                [mars._mass_dl[-1], mars.pressure[-1],
                                 mars.temperature[-1]],
                                r_c, max_step=mars.grid_step)
        while mantle_distr.status == 'running':
            mantle_distr.step()
            mars._radius_dl.append(mantle_distr.t)
            mars._density_dl.append(den_mantle_func(mantle_distr.t, mantle_distr.y))
            mars._mass_dl.append(mantle_distr.y[0])
            mars.pressure.append(mantle_distr.y[1])
            mars.temperature.append(temp_mantle_func(mantle_distr.t, mantle_distr.y))

        # 5 - shift the boundary while the mass at the center is to be zero
        # step 1 - integrating the large grid
        # step 2 - integrating the fine grid
        steps_in_core = [np.inf, mars.grid_step]
        print("10% - calculating of the core")
        for i in range(len(steps_in_core)):
            distr_core = {'density': [den_core_func(mars._radius_dl[-1] - mars.grid_step ** 2,
                                                    np.array([mars._mass_dl[-1],
                                                     mars.pressure[-1] * (1 + mars.grid_step ** 2),
                                                     mars.temperature[-1]]))],
                          'radius': [mars._radius_dl[-1] - mars.grid_step ** 2],
                          'mass': [mars._mass_dl[-1]],
                          'pressure': [mars.pressure[-1] * (1 + mars.grid_step ** 2)],
                          'temperature': [mars.temperature[-1]]}
            core_distr = ing.RK45(sys15_core, distr_core['radius'][-1],
                                  [distr_core['mass'][-1], distr_core['pressure'][-1], distr_core['temperature'][-1]],
                                  mars.t_core_bound, max_step=steps_in_core[i])

            if i == 1:
                step = mars._radius_dl[-1]/4
                now = 3
            while core_distr.status == 'running':
                if i == 1 and core_distr.t < step*now and now > 0:
                    print(f"{50-10*now}% - calculating of the core")
                    now -= 1
                core_distr.step()
                distr_core['radius'].append(core_distr.t)
                distr_core['density'].append(den_core_func(core_distr.t, core_distr.y))
                distr_core['mass'].append(core_distr.y[0])
                distr_core['pressure'].append(core_distr.y[1])
                distr_core['temperature'].append(core_distr.y[2])

            distr_core['radius'].append(0)
            distr_core['density'].append(distr_core['density'][-1])
            distr_core['mass'].append(distr_core['mass'][-1] - distr_core['density'][-1] * mars.t_core_bound ** 3)
            distr_core['pressure'].append(
                distr_core['pressure'][-1] + mars.P_const * (distr_core['density'][-1] * mars.t_core_bound) ** 2 / 2)
            distr_core['temperature'].append(
                distr_core['temperature'][-1] * m.exp(mars.T_exp_const * distr_core['density'][-1] * mars.t_core_bound ** 2 / 2))

            mars.l_to_core = len(mars._radius_dl)
            if i == 1:
                print("50% - calculating of the core")

            mass_center_count = 0
            if distr_core['mass'][-1] > 0:
                mass_center_count = 1
            elif distr_core['mass'][-1] < 0:
                mass_center_count = -1

            distr_numeric_last = {}
            while mass_center_count:
                if mass_center_count > 0:
                    distr_numeric_last['radius'] = mars._radius_dl[-1]
                    distr_numeric_last['density'] = mars._density_dl[-1]
                    distr_numeric_last['mass'] = mars._mass_dl[-1]
                    distr_numeric_last['pressure'] = mars.pressure[-1]
                    distr_numeric_last['temperature'] = mars.temperature[-1]
                    for l in [mars._radius_dl, mars._density_dl, mars._mass_dl, mars.pressure, mars.temperature]:
                        l.pop()
                else:
                    mantle_distr.t_bound -= mars.grid_step
                    mantle_distr.status = 'running'
                    while mantle_distr.status == 'running':
                        mantle_distr.step()
                        mars._radius_dl.append(mantle_distr.t)
                        mars._density_dl.append(den_mantle_func(mantle_distr.t, mantle_distr.y))
                        mars._mass_dl.append(mantle_distr.y[0])
                        mars.pressure.append(mantle_distr.y[1])
                        mars.temperature.append(temp_mantle_func(mantle_distr.t, mantle_distr.y))

                distr_core_prev = distr_core
                distr_core = {'density': [den_core_func(mars._radius_dl[-1] - mars.grid_step ** 2,
                                                        np.array([mars._mass_dl[-1],
                                                         mars.pressure[-1] * (1 + mars.grid_step ** 2),
                                                         mars.temperature[-1]]))],
                              'radius': [mars._radius_dl[-1] - mars.grid_step ** 2],
                              'mass': [mars._mass_dl[-1]],
                              'pressure': [mars.pressure[-1] * (1 + mars.grid_step ** 2)],
                              'temperature': [mars.temperature[-1]]}
                core_distr = ing.RK45(sys15_core, distr_core['radius'][-1],
                                      [distr_core['mass'][-1], distr_core['pressure'][-1],
                                       distr_core['temperature'][-1]],
                                      mars.t_core_bound, max_step=steps_in_core[i])

                if i == 1:
                    step = mars._radius_dl[-1] / 4
                    now = 3
                while core_distr.status == 'running':
                    if i == 1 and core_distr.t < step * now and now > 0:
                        print(f"{90 - 10 * now}% - calculating of the core")
                        now -= 1
                        if now == -1:
                            print(f"{90}% - need some more time")
                    core_distr.step()
                    distr_core['radius'].append(core_distr.t)
                    distr_core['density'].append(den_core_func(core_distr.t, core_distr.y))
                    distr_core['mass'].append(core_distr.y[0])
                    distr_core['pressure'].append(core_distr.y[1])
                    distr_core['temperature'].append(core_distr.y[2])

                distr_core['radius'].append(0)
                distr_core['density'].append(distr_core['density'][-1])
                distr_core['mass'].append(distr_core['mass'][-1] - distr_core['density'][-1] * mars.t_core_bound ** 3)
                distr_core['pressure'].append(distr_core['pressure'][-1] +
                                              mars.P_const * (distr_core['density'][-1] * mars.t_core_bound) ** 2 / 2)
                distr_core['temperature'].append(distr_core['temperature'][-1] *
                                                 m.exp(mars.T_exp_const * distr_core['density'][-1] * mars.t_core_bound ** 2 / 2))

                mars.l_to_core = len(mars._radius_dl)

                if mass_center_count * distr_core['mass'][-1] <= 0:
                    if abs(distr_core_prev['mass'][-1]) < abs(distr_core['mass'][-1]):
                        distr_core = distr_core_prev
                        if mass_center_count > 0:
                            mars._radius_dl.append(distr_numeric_last['radius'])
                            mars._density_dl.append(distr_numeric_last['density'])
                            mars._mass_dl.append(distr_numeric_last['mass'])
                            mars.pressure.append(distr_numeric_last['pressure'])
                            mars.temperature.append(distr_numeric_last['temperature'])
                        else:
                            while mars._radius_dl[-1] < mantle_distr.t_bound + mars.grid_step:
                                for l in [mars._radius_dl, mars._density_dl, mars._mass_dl, mars.pressure, mars.temperature]:
                                    l.pop()
                        mars.l_to_core = len(mars._radius_dl)
                    mass_center_count = 0

        for i in range(len(distr_core['mass'])):
            distr_core['mass'][i] -= distr_core['mass'][-1]

        # Calculating for the seismic wave speed

        # 1 - Crust
        mars.vp = [7.0] * mars.l_crust
        mars.vs = [4.0] * mars.l_crust

        # 2 - Mantle
        vel_p_mantle_func = inp.interp1d(DATAs_perplex['pressure'], DATAs_perplex['compressional velocity'],
                                         fill_value='extrapolate')
        vel_s_mantle_func = inp.interp1d(DATAs_perplex['pressure'], DATAs_perplex['shear velocity'],
                                         fill_value='extrapolate')
        for i in range(mars.l_crust, mars.l_to_core):
            mars.vp.append(vel_p_mantle_func(mars.pressure[i]))
            mars.vs.append(vel_s_mantle_func(mars.pressure[i]))

        # 3 - Core
        distr_core['compr velocity'] = []
        distr_core['shear velocity'] = []

        class Core_elasticity:
            def __init__(self, mode):
                self.mode = mode
                # parameter lists [den_0, alpha_0, K_T_0, K'_T_0, Kdot_T_0, T_0]
                if mode == 'Fe':
                    self.param = [7.03 / mars.den_av, 75e-6, 105, 4.5, 0.025, 2100]
                elif mode == 'FeS':
                    self.param = [4.94 / mars.den_av, 68.52e-6, 54, 4, 0.02, 1100]
                elif mode == 'FeH':
                    self.param = [6.7 / mars.den_av, 1, 121, 5.31, 0, 0]

            def __call__(self, T):
                if self.mode == 'FeH':
                    self.param[5] = T
                den_star = self.param[0] * m.exp(-self.param[1] * (T - self.param[5]))
                expon = self.param[4] / self.param[1] / self.param[2]
                K = self.param[2] * (den_star / self.param[0]) ** expon
                K_der = self.param[3] * m.exp(self.param[1] * (T - self.param[5]))
                return den_star, K, K_der

        Fe_vel = Core_elasticity('Fe')
        FeS_vel = Core_elasticity('FeS')
        FeH_vel = Core_elasticity('FeH')

        for i in range(len(distr_core['radius'])):
            T = distr_core['temperature'][i]

            den_fe, K_fe, K_der_fe = Fe_vel(T)
            den_fes, K_fes, K_der_fes = FeS_vel(T)
            den_feh, K_feh, K_der_feh = FeH_vel(T)

            denv = x_fe * den_fe + x_fes * den_fes + x_feh * den_feh
            denr = x_fe / den_fe + x_fes / den_fes + x_feh / den_feh
            Kv = x_fe * K_fe + x_fes * K_fes + x_feh * K_feh
            Kr = x_fe / K_fe + x_fes / K_fes + x_feh / K_feh
            Kv_der = x_fe * K_der_fe + x_fes * K_der_fes + x_feh * K_der_feh
            Kr_der = x_fe / K_der_fe + x_fes / K_der_fes + x_feh / K_der_feh

            den_star = (denv + 1 / denr) / 2
            K = (Kv + 1 / Kr) / 2
            K_der = (Kv_der + 1 / Kr_der) / 2

            L1 = K
            L2 = 5 * K - 3 * K * K_der
            epsil = (1 - (distr_core['density'][i] / den_star) ** (2 / 3)) / 2

            distr_core['compr velocity'].append(
                m.sqrt((1 - 2 * epsil) ** (5 / 2) * (L1 + L2 * epsil) / distr_core['density'][i] / mars.den_av))
            distr_core['shear velocity'].append(0)

        mars._radius_dl += distr_core['radius']
        mars._density_dl += distr_core['density']
        mars._mass_dl += distr_core['mass']
        mars.pressure += distr_core['pressure']
        mars.temperature += distr_core['temperature']
        mars.vp += distr_core['compr velocity']
        mars.vs += distr_core['shear velocity']

        mars.l_planet = len(mars._radius_dl)
        mars.radius = [mars._radius_dl[i]*mars.R_Mars for i in range(mars.l_planet)]
        mars.density = [mars._density_dl[i]*mars.den_av for i in range(mars.l_planet)]
        mars.mass = [mars._mass_dl[i]*mars.M_Mars for i in range(mars.l_planet)]

        mars.gravity = [mars.Gravity_const * mars.mass[i] / (mars.radius[i] * 1e3) ** 2 for i in range(mars.l_planet - 1)]
        mars.gravity.append(0)

        print("95% - finishing the calculation")

        mars.radius_core = mars.radius[mars.l_to_core]
        mars.density_core_boundary = mars.density[mars.l_to_core]
        mars.density_core = mars.density[-1]
        print('100% - succesfully completed')
        print('----------------------------')

    def cid(mars):
        """'Calculate Inelastic Distributions'

        Calculates the distributions of inelastic parameters such as lame parameters"""

        # Initialization of the main parameters of the method
        print('\n')
        print('---------------------------------------------------------')
        print('Initialization the calculation of inelastic distributions')
        print('---------------------------------------------------------')
        if mars.bviscosity:
            distr_mineral = pd.read_excel(
                '../data/dynamic/mineral_distributions_' + mars.composition + '_' + mars.areoterm + '.xlsx')

        mars.lame1_elast = []
        mars.lame2_elast = []
        mars.lame1 = []
        mars.lame2 = []
        mars.lame1_cw = []
        mars.lame2_cw = []
        mars.Q = []
        mars.visc_distr = []

        # Introducing inelasticity
        if mars.bviscosity:
            # 1 - Inelastic parameters
            visc_0_exp = m.log10(mars.viscosity) - 9  # viscosity of the crust, exponent in GPa*s
            visc_melt_exp = m.log10(mars.viscosity_melt) - 9

            # 2 - Calculation of the layers
            rad_trans = [mars.R_Mars, mars.R_Mars - mars.depth_crust, mars.R_Mars - mars.depth_crust -
                         mars.grid_step ** 2]
            visc_exp = [visc_0_exp, visc_0_exp, visc_0_exp - 2]
            pres_trans = []
            for i in range(1, distr_mineral.shape[0]):
                if len(pres_trans) == 0 and not m.isnan(distr_mineral['Wad'][i]) and m.isnan(
                        distr_mineral['Wad'][i - 1]):
                    pres_trans += [distr_mineral['pressure'][i] / 1e4]
                if len(pres_trans) == 1 and m.isnan(distr_mineral['O'][i]) and not m.isnan(distr_mineral['O'][i - 1]):
                    pres_trans += [distr_mineral['pressure'][i] * (1 + mars.grid_step ** 2) / 1e4]
                if len(pres_trans) == 2 and not m.isnan(distr_mineral['Ring'][i]) and m.isnan(
                        distr_mineral['Ring'][i - 1]):
                    pres_trans += [distr_mineral['pressure'][i] / 1e4]
                if len(pres_trans) == 3 and m.isnan(distr_mineral['Wad'][i]) and not m.isnan(
                        distr_mineral['Wad'][i - 1]):
                    pres_trans += [distr_mineral['pressure'][i] * (1 + mars.grid_step ** 2) / 1e4]
                    break
            pres_to_rad = inp.interp1d(mars.pressure, mars.radius)
            visc_exp += [visc_0_exp - 2, visc_0_exp - 1, visc_0_exp - 1, visc_0_exp]
            for p in pres_trans:
                r = pres_to_rad(p)
                if r < mars.radius_core + mars.H_melt * mars.bmeltlayer:
                    break
                else:
                    rad_trans += [r]
            visc_exp = visc_exp[:len(rad_trans)]

            visc_exp += [visc_exp[-1]]
            if mars.bmeltlayer:
                rad_trans += [mars.radius_core + mars.H_melt, mars.radius_core + mars.H_melt - mars.grid_step ** 2]
                visc_exp += [visc_melt_exp] * 2
            rad_trans += [mars.radius_core]

            def visc(r):
                fun = inp.interp1d(rad_trans, visc_exp)
                return 10 ** (fun(r))

            # 3 - Calculating of the Lame parameters
            for i in range(mars.l_to_core):
                lame2_elast = (mars.vs[i] ** 2) * mars.density[i]
                mars.lame2_elast.append(lame2_elast)
                J = ((1 + (1j * visc(mars.radius[i]) / lame2_elast * mars.freq_sun) ** (-mars.andrade) * m.gamma(
                    1 + mars.andrade))
                     / lame2_elast
                     - 1j / (visc(mars.radius[i]) * mars.freq_sun))
                lame2_inelast = 1 / J
                mars.lame2.append(lame2_inelast)
                bulk_elast = (mars.vp[i] ** 2) * mars.density[i] - 4 / 3 * \
                             mars.lame2_elast[i]
                bulk_inelast = (mars.vp[i] ** 2) * mars.density[i] - 4 / 3 * \
                               mars.lame2[i]
                mars.lame1_elast.append(bulk_elast - 2 / 3 * lame2_elast)
                mars.lame1.append(bulk_inelast - 2 / 3 * lame2_inelast)

                J_cw = ((1 + (1j * visc(mars.radius[i]) / lame2_elast * mars.freq_cw) ** (-mars.andrade) * m.gamma(
                    1 + mars.andrade))
                        / lame2_elast
                        - 1j / (visc(mars.radius[i]) * mars.freq_cw))
                mars.lame2_cw.append(1 / J_cw)
                bulk_cw = ((mars.vp[i] ** 2) * mars.density[i] -
                           4 / 3 * mars.lame2_cw[i])
                mars.lame1_cw.append(bulk_cw - 2 / 3 / J_cw)

                mars.Q.append(lame2_inelast.real / lame2_inelast.imag)
                mars.visc_distr.append(visc(mars.radius[i]) * 1e9)

        # Calculating of the Lame parameters for elastic model
        else:
            for i in range(mars.l_to_core):
                mars.lame2.append((mars.vs[i] ** 2) * mars.density[i])
                mars.lame2_elast.append(mars.lame2[i])
                mars.lame1.append((mars.vp[i] ** 2) * mars.density[i] -
                                            4 / 3 * mars.lame2[i])
                mars.lame1_elast.append(mars.lame1[i])

                mars.lame2_cw.append(mars.lame2[i])
                mars.lame1_cw.append(mars.lame1[i])

                mars.Q.append(m.inf)
                mars.visc_distr.append(0)

        for i in range(mars.l_to_core, mars.l_planet):
            mars.lame2.append(0)
            mars.lame2_elast.append(mars.lame2[i])
            mars.lame1.append((mars.vp[i] ** 2) * mars.density[i])
            mars.lame1_elast.append(mars.lame1[i])

            mars.lame2_cw.append(0)
            mars.lame1_cw.append(mars.lame1[i])

            mars.Q.append(m.inf)
            mars.visc_distr.append(0)
        print('Finishing the calculation')
        print('-------------------------')
        pass

    def cip(mars):
        """'Calculate Integral Parameters'

        Calculates integral parameters of the model"""

        n = 2  # order of the defining Love number
        print('\n')
        print('-----------------------------------------------------')
        print('Initialization the calculation of integral parameters')
        print('-----------------------------------------------------')

        # Making dimensionless distributions
        mars._gravity_dl = [mars.gravity[i] / mars.grav_av for i in range(mars.l_planet)]
        mars._lame1_elast_dl = [mars.lame1_elast[i] / (mars.den_av * mars.R_Mars * mars.grav_av / 1e3) for i in range(mars.l_planet)]
        mars._lame2_elast_dl = [mars.lame2_elast[i] / (mars.den_av * mars.R_Mars * mars.grav_av / 1e3) for i in range(mars.l_planet)]
        mars._lame1_dl = [mars.lame1[i] / (mars.den_av * mars.R_Mars * mars.grav_av / 1e3) for i in range(mars.l_planet)]
        mars._lame2_dl = [mars.lame2[i] / (mars.den_av * mars.R_Mars * mars.grav_av / 1e3) for i in range(mars.l_planet)]
        mars._lame1_cw_dl = [mars.lame1_cw[i] / (mars.den_av * mars.R_Mars * mars.grav_av / 1e3) for i in range(mars.l_planet)]
        mars._lame2_cw_dl = [mars.lame2_cw[i] / (mars.den_av * mars.R_Mars * mars.grav_av / 1e3) for i in range(mars.l_planet)]
        mars._radius_core_dl = mars.radius_core / mars.R_Mars

        # Calculation of the model values of mass and moment of inertia with the use of Simpson formula
        integral_dict = {'mass': [], 'inertia': []}
        for i in range(mars.l_planet):
            integral_dict['mass'] += [4 * m.pi * mars._density_dl[i] * mars._radius_dl[i] ** 2]
            integral_dict['inertia'] += [8 / 3 * m.pi * mars._density_dl[i] * mars._radius_dl[i] ** 4]
        mars.mass_model = - ing.trapezoid(integral_dict['mass'], x=mars._radius_dl)
        mars.mass_model *= 3 / 4 / m.pi
        mars.MOI = - ing.trapezoid(integral_dict['inertia'], x=mars._radius_dl)
        mars.MOI *= 3 / 4 / m.pi
        mars.MOI_core = - ing.trapezoid(integral_dict['inertia'][mars.l_to_core:], x=mars._radius_dl[mars.l_to_core:])
        mars.MOI_core *= 3 / 4 / m.pi

        # Calculation of the Molodensky number
        print('0% - calculation of the Molodensky number')
        den_func = inp.interp1d(mars._radius_dl, mars._density_dl)
        grav_func = inp.interp1d(mars._radius_dl, mars._gravity_dl)

        def Molodensky_func(t, y):
            return - y[0] ** 2 - 2 * (n + 1) / t * y[0] - 3 * (
                    den_func(t + 0.0001) - den_func(t - 0.0001)) / 0.0002 / grav_func(t)

        Molodensky_solution = ing.RK45(Molodensky_func, mars.t_core_bound, [0], mars._radius_dl[mars.l_to_core] - 0.00011,
                                       max_step=mars.grid_step)
        while Molodensky_solution.status == 'running':
            Molodensky_solution.step()
        mars.Molodensky_number = Molodensky_solution.y[0]

        # Calculation of the Love number k2 for 3 different frequency
        lame1_func = inp.interp1d(mars._radius_dl, mars._lame1_dl)
        lame2_func = inp.interp1d(mars._radius_dl, mars._lame2_dl)
        lame1_cw_func = inp.interp1d(mars._radius_dl, mars._lame1_cw_dl)
        lame2_cw_func = inp.interp1d(mars._radius_dl, mars._lame2_cw_dl)
        lame1_elst_func = inp.interp1d(mars._radius_dl, mars._lame1_elast_dl)
        lame2_elst_func = inp.interp1d(mars._radius_dl, mars._lame2_elast_dl)

        class Love_RHS:
            def __init__(self, freq_mode):
                if freq_mode == 'sun':
                    self.lame1 = lame1_func
                    self.lame2 = lame2_func
                elif freq_mode == 'cw':
                    self.lame1 = lame1_cw_func
                    self.lame2 = lame2_cw_func
                elif freq_mode == 'elast':
                    self.lame1 = lame1_elst_func
                    self.lame2 = lame2_elst_func

            def __call__(self, x, y):
                lame1 = self.lame1(x)
                lame2 = self.lame2(x)
                den = den_func(x)
                grav = grav_func(x)
                eq = [0] * 6
                eq[0] = -2 * lame1 / (lame1 + 2 * lame2) * y[0] / x + y[1] / (lame1 + 2 * lame2) + lame1 * n * (
                            n + 1) / (lame1 + 2 * lame2) * y[4] / x
                eq[1] = ((-4 * den * grav * x + 4 * lame2 * (3 * lame1 + 2 * lame2) / (lame1 + 2 * lame2)) * y[
                    0] / x ** 2 - 4 * lame2 / (lame1 + 2 * lame2) * y[1] / x +
                         (n * (n + 1) * den * grav * x - 2 * lame2 * (3 * lame1 + 2 * lame2) * n * (n + 1) / (
                                     lame1 + 2 * lame2)) * y[4] / x ** 2 + n * (n + 1) / x * y[5] -
                         den * y[3])
                eq[2] = 3 * den * y[0] + y[3]
                eq[3] = -3 * den * (n + 1) * n / x * y[4] + n * (n + 1) / x ** 2 * y[2] - 2 / x * y[3]
                eq[4] = -y[0] / x + y[4] / x + y[5] / lame2
                eq[5] = ((den * grav / x - 2 * lame2 * (3 * lame1 + 2 * lame2) / (lame1 + 2 * lame2) / x ** 2) * y[0] -
                         lame1 / (lame1 + 2 * lame2) * y[1] / x +
                         2 * lame2 / (lame1 + 2 * lame2) * (
                                     lame1 * (2 * n ** 2 + 2 * n - 1) + 2 * lame2 * (n ** 2 + n - 1)) * y[
                             4] / x ** 2 - 3 * y[5] / x -
                         den * y[2] / x)
                return eq

        # 1 - Elastic Love number
        Love_RHS_elast = Love_RHS('elast')
        Love_elast_solution = [ing.RK45(Love_RHS_elast, mars._radius_dl[mars.l_to_core - 1],
                                        np.array([1,
                                                  mars._density_dl[mars.l_to_core - 1] * mars._gravity_dl[
                                                      mars.l_to_core - 1],
                                                  0, -3 * mars._density_dl[mars.l_to_core - 1], 0, 0],
                                                 dtype="complex_"),
                                        1, max_step=mars.grid_step),
                               ing.RK45(Love_RHS_elast, mars._radius_dl[mars.l_to_core - 1],
                                        np.array([0, -mars._density_dl[mars.l_to_core - 1], 1,
                                                  n / mars._radius_dl[mars.l_to_core - 1] + mars.Molodensky_number, 0, 0],
                                                 dtype="complex_"),
                                        1, max_step=mars.grid_step),
                               ing.RK45(Love_RHS_elast, mars._radius_dl[mars.l_to_core - 1],
                                        np.array([0, 0, 0, 0, 1, 0], dtype="complex_"), 1, max_step=mars.grid_step)]
        now = 0
        for sol in Love_elast_solution:
            print(f'{10+10*now}% - calculation of the elastic Love number')
            while sol.status == 'running':
                sol.step()
            now += 1

        bound_cond_matrix = [[Love_elast_solution[0].y[1], Love_elast_solution[1].y[1], Love_elast_solution[2].y[1]],
                             [Love_elast_solution[0].y[3] + (n + 1) * Love_elast_solution[0].y[2],
                              Love_elast_solution[1].y[3] + (n + 1) * Love_elast_solution[1].y[2],
                              Love_elast_solution[2].y[3] + (n + 1) * Love_elast_solution[2].y[2]],
                             [Love_elast_solution[0].y[5], Love_elast_solution[1].y[5], Love_elast_solution[2].y[5]]]
        bound_cond_RHS = [0, 2 * n + 1, 0]
        solutions_coeffs = lin.solve(bound_cond_matrix, bound_cond_RHS)
        mars.Love_elast = (
                    solutions_coeffs[0] * Love_elast_solution[0].y[2] + solutions_coeffs[1] * Love_elast_solution[1].y[
                2] +
                    solutions_coeffs[2] * Love_elast_solution[2].y[2] - 1)

        # 2 - Love number on the Chandler frequency
        Love_RHS_CW = Love_RHS('cw')
        Love_CW_solution = [ing.RK45(Love_RHS_CW, mars._radius_dl[mars.l_to_core - 1],
                                     np.array([1,
                                               mars._density_dl[mars.l_to_core - 1] * mars._gravity_dl[
                                                   mars.l_to_core - 1],
                                               0, -3 * mars._density_dl[mars.l_to_core - 1], 0, 0],
                                              dtype="complex_"),
                                     1, max_step=mars.grid_step),
                            ing.RK45(Love_RHS_CW, mars._radius_dl[mars.l_to_core - 1],
                                     np.array([0, -mars._density_dl[mars.l_to_core - 1], 1,
                                               n / mars._radius_dl[mars.l_to_core - 1] + mars.Molodensky_number, 0, 0],
                                              dtype="complex_"), 1, max_step=mars.grid_step),
                            ing.RK45(Love_RHS_CW, mars._radius_dl[mars.l_to_core - 1], np.array([0, 0, 0, 0, 1, 0],
                                                                                                   dtype="complex_"),
                                     1, max_step=mars.grid_step)]
        now = 0
        for sol in Love_CW_solution:
            print(f'{40+10*now}% - calculation of the Love number on the Chandler frequency')
            while sol.status == 'running':
                sol.step()
            now += 1

        bound_cond_matrix = [[Love_CW_solution[0].y[1], Love_CW_solution[1].y[1], Love_CW_solution[2].y[1]],
                             [Love_CW_solution[0].y[3] + (n + 1) * Love_CW_solution[0].y[2],
                              Love_CW_solution[1].y[3] + (n + 1) * Love_CW_solution[1].y[2],
                              Love_CW_solution[2].y[3] + (n + 1) * Love_CW_solution[2].y[2]],
                             [Love_CW_solution[0].y[5], Love_CW_solution[1].y[5], Love_CW_solution[2].y[5]]]
        bound_cond_RHS = [0, 2 * n + 1, 0]
        solutions_coeffs = lin.solve(bound_cond_matrix, bound_cond_RHS)
        mars.Love_cw = (solutions_coeffs[0] * Love_CW_solution[0].y[2] + solutions_coeffs[1] * Love_CW_solution[1].y[2] +
                   solutions_coeffs[2] * Love_CW_solution[2].y[2] - 1)

        # 3 - Love number on the Sun tides frequency
        Love_RHS_Sun = Love_RHS('sun')
        Love_Sun_solution = [ing.RK45(Love_RHS_Sun, mars._radius_dl[mars.l_to_core - 1],
                                      np.array([1,
                                                mars._density_dl[mars.l_to_core - 1] * mars._gravity_dl[
                                                    mars.l_to_core - 1],
                                                0, -3 * mars._density_dl[mars.l_to_core - 1], 0, 0],
                                               dtype="complex_"),
                                      1, max_step=mars.grid_step),
                             ing.RK45(Love_RHS_Sun, mars._radius_dl[mars.l_to_core - 1],
                                      np.array([0, -mars._density_dl[mars.l_to_core - 1], 1,
                                                n / mars._radius_dl[mars.l_to_core - 1] + mars.Molodensky_number, 0, 0],
                                               dtype="complex_"),
                                      1, max_step=mars.grid_step),
                             ing.RK45(Love_RHS_Sun, mars._radius_dl[mars.l_to_core - 1],
                                      np.array([0, 0, 0, 0, 1, 0], dtype="complex_"), 1, max_step=mars.grid_step)]
        now = 0
        for sol in Love_Sun_solution:
            print(f'{70 + 10 * now}% - calculation of the Love number on the Sun tides frequency')
            while sol.status == 'running':
                sol.step()
            now += 1

        bound_cond_matrix = [[Love_Sun_solution[0].y[1], Love_Sun_solution[1].y[1], Love_Sun_solution[2].y[1]],
                             [Love_Sun_solution[0].y[3] + (n + 1) * Love_Sun_solution[0].y[2],
                              Love_Sun_solution[1].y[3] + (n + 1) * Love_Sun_solution[1].y[2],
                              Love_Sun_solution[2].y[3] + (n + 1) * Love_Sun_solution[2].y[2]],
                             [Love_Sun_solution[0].y[5], Love_Sun_solution[1].y[5], Love_Sun_solution[2].y[5]]]
        bound_cond_RHS = [0, 2 * n + 1, 0]
        solutions_coeffs = lin.solve(bound_cond_matrix, bound_cond_RHS)
        mars.Love_sun = (solutions_coeffs[0] * Love_Sun_solution[0].y[2] + solutions_coeffs[1] * Love_Sun_solution[1].y[2] +
                    solutions_coeffs[2] * Love_Sun_solution[2].y[2] - 1)

        # Calculation of the Chandler period
        A = 0.362976
        B = 0.363229
        C = 0.365067
        period_rot = 24.6229 * 3600

        A_red = (C - B) / A
        B_red = (C - A) / B
        A_aver = (A + B) / 2
        freq_rot = 2 * m.pi / period_rot
        period_Euler = period_rot / m.sqrt(A_red * B_red)
        mars.Love_sec = 3 * mars.Gravity_const * (C - A_aver) * mars.M_Mars / freq_rot ** 2 / (mars.R_Mars * 1e3) ** 3

        mars.period_cw = period_Euler * (1 - mars.MOI_core / m.sqrt(A * B)) / (1 - mars.Love_cw.real / mars.Love_sec) / 86400

        mars.period_cw_elast = period_Euler * (1 - mars.MOI_core / m.sqrt(A * B)) / (1 - mars.Love_elast.real / mars.Love_sec) / 86400

        # Calculation of the Chandler period with the use of creep function
        if mars.bviscosity and mars.bcreepfunction:
            mars.period_cw_creep = period_Euler * (1 - mars.MOI_core / m.sqrt(A * B)) / (1 - mars.Love_sun.real / mars.Love_sec) / 86400
            mars.period_cw_creep_delta = mars.period_cw_creep / 1118 * ((mars.freq_sun / mars.freq_cw) ** mars.creep_exp - 1) / (mars.Love_sec - mars.Love_cw.real) / m.tan(
                mars.creep_exp * m.pi / 2)

        print('100% - successfully completed')
        print('---------------------------')
        print('Results of the calculation:')
        print(f'Mass,                                    M = {mars.mass_model:.4f}')
        print(f'Moment of inertia,                       I = {mars.MOI:.4f}')
        print(f'Molodensky number,                   gamma = {mars.Molodensky_number:.4f}')
        print(f'Love number, elastic,                 k2_e = {mars.Love_elast:.4f}')
        print(f'Love number on Chandler frequency,   k2_cw = {mars.Love_cw:.4f}')
        print(f'Love number on Sun tides frequency, k2_sun = {mars.Love_sun:.4f}')
        print(f'Core radius,                        R_core = {mars.radius_core:.0f}')
        print(f'Chandler period, elastic,             Tw0e = {mars.period_cw_elast:.1f}')
        print(f'Chandler period,                       Tw0 = {mars.period_cw:.1f}')
        if mars.bviscosity and mars.bcreepfunction:
            print(f'Chandler period, by creep function,  Tw_cr = {mars.period_cw_creep:.1f} + '
                  f'{mars.period_cw_creep_delta:.1f}')

        print('-------------------------')
        print('Finishing the calculation')
        print('-------------------------')
        pass

    def calculate(mars):
        mars.ced()
        mars.cid()
        mars.cip()
        pass

    def plot_d(mars,
               show: bool=True,
               save: bool=False,
               name_fig: str="test.png"):
        """Plots calculated distributions

        :param show: Defines whether to show the graph, defaults to True
        :type show: bool
        :param save: Defines whether to save the graph, defaults to False
        :type save: bool
        :param name_fig: Is used only when save = True!
            The name of graph to save with extension, defaults to "test.png"
        :type name_fig: str"""
        fig, ax = plt.subplots()
        fig.subplots_adjust(left=0.25, right=0.75)

        pres = ax.twinx()
        temp = ax.twinx()
        velo = ax.twinx()

        velo.yaxis.tick_left()
        velo.yaxis.set_label_position('left')
        temp.spines.right.set_position(("axes", 1.2))
        velo.spines.left.set_position(("axes", -0.2))

        pr, = ax.plot(mars.radius, mars.density, "C0")
        ax.text(700, mars.density[mars.l_planet - mars.l_planet // 10] + 0.1, r"$\mathrm{\rho}$")
        pp, = pres.plot(mars.radius, mars.pressure, "C1")
        pres.text(1000, mars.pressure[mars.l_planet - 3 * mars.l_planet // 10], "$P$",
                  horizontalalignment='right', verticalalignment='top')
        pt, = temp.plot(mars.radius, mars.temperature, "C2")
        temp.text(1500, mars.temperature[mars.l_planet - 4 * mars.l_planet // 10] + 30, '$T$')
        pvp, = velo.plot(mars.radius, mars.vp, "C3")
        velo.text(700, mars.vp[mars.l_planet - mars.l_planet // 10] + 0.3, r"$v_\mathrm{p}$")
        velo.text(3000, mars.vp[mars.l_planet - 9 * mars.l_planet // 10] - 0.1, r"$v_\mathrm{p}$",
                  horizontalalignment='right', verticalalignment='top')
        pvs, = velo.plot(mars.radius, mars.vs, "C4")
        velo.text(2200, mars.vs[mars.l_planet - 7 * mars.l_planet // 10], r"$v_\mathrm{s}$",
                  horizontalalignment='right', verticalalignment='top')

        ax.set(xlim=(0, mars.R_Mars), ylim=(0, 10), xlabel="Радиус, $км$", ylabel="Плотность, $г/см^3$")
        pres.set(ylim=(0, 50), ylabel="Давление, $ГПа$")
        temp.set(ylim=(0, 2500), ylabel="Температура, $К$")
        velo.set(ylim=(0, 25), ylabel="Скорость, $км/с$")

        if save:
            name_fig = '../data/png/' + name_fig
            plt.savefig(fname=name_fig, dpi=600, bbox_inches='tight', pad_inches=0.1)
        if show:
            plt.show()

    def save_distributions(mars,
                           name_file: str='model_distributions'):
        """Saves the data of distribution calculated

        :param name_file: The name of the file to save without extension, defaults to "model_distributions"
        :type name_file: str"""
        DATAs = pd.DataFrame({'radius': mars.radius,
                              'density': mars.density,
                              'gravity': mars.gravity,
                              'mass': mars.mass,
                              'pressure': mars.pressure,
                              'temperature': mars.temperature,
                              'compr velocity': mars.vp,
                              'shear velocity': mars.vs,
                              'viscosity': mars.visc_distr,
                              'lame1_elast': mars.lame1_elast,
                              'lame2_elast': mars.lame2_elast,
                              'lame1_Sun': mars.lame1,
                              'lame2_Sun': mars.lame2,
                              'lame1_CW': mars.lame1_cw,
                              'lame2_CW': mars.lame2_cw,
                              'Q': mars.Q})
        DATAs.to_excel("../data/dynamic/" + name_file + ".xlsx", index=False)

    def save_integrals(mars,
                       target: str="dynamic/integral_parameters",
                       rewrite: bool=False,
                       add_composition_name: bool=False,
                       suffix: str=""):
        """Saves the data of integrals calculated

        :param target: The path and name of the file to save, defaults to "dynamic/integral_parameters"
        :type target: str
        :param rewrite: Defines whether to rewrite the file, defaults to False
        :type rewrite: bool
        :param add_composition_name: Defines whether to add composition name (and areotherm name) to the target file
            name, defaults to False
        :type add_composition_name: bool
        :param suffix: The name of the suffix to add to the name, defaults to ""
        :type suffix: str"""
        integrals = {
            'crust_density': [mars.density_crust],
            'crust_depth': [mars.depth_crust],
            'core_density_center': [mars.density_core],
            'core_density_boundary': [mars.density_core_boundary],
            'core_radius': [mars.radius_core],
            'FeH': [mars.hydro_core],
            'FeS': [mars.sulfur_core],
            'has_viscosity': [mars.bviscosity],
            'viscosity': [mars.viscosity],
            'andrade': [mars.andrade],
            'has_melt_layer': [mars.bmeltlayer],
            'viscosity_melt': [mars.viscosity_melt],
            'inertia': [mars.MOI],
            'Molodensky_number': [mars.Molodensky_number],
            'Love_number_elastic': [mars.Love_elast.real],
            'Love_number': [mars.Love_sun],
            'Love_number_CW': [mars.Love_cw],
            'Love_number_secular': [mars.Love_sec],
            'Chandler_period': [mars.period_cw],
            'Chandler_period_elastic': [mars.period_cw_elast]
        }
        df_integrals = pd.DataFrame(integrals)
        s_output = '../data/' + target
        # PhD_2nd_year/comparing_chem_models
        if add_composition_name:
            s_output += '_' + mars.composition + '_' + mars.areoterm
        if suffix:
            s_output += '_' + suffix
        s_output += '.xlsx'
        if not rewrite:
            df_integrals_old = pd.read_excel(s_output)
            df_integrals = pd.concat([df_integrals_old, df_integrals], ignore_index=True)
        df_integrals.to_excel(s_output, index=False)
        pass